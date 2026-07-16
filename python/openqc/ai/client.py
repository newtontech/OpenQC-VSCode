"""Dependency-free OpenAI Responses API and Ollama client for OpenQC."""

from __future__ import annotations

import argparse
import json
import os
import re
import socket
import sys
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from typing import Any, Dict, Optional

from .parser import ResponseParser
from .prompts import PromptTemplates

DEFAULT_TIMEOUT_SECONDS = 120.0
DEFAULT_MAX_OUTPUT_CHARS = 100_000


@dataclass
class AIRequest:
    """AI request data structure."""

    type: str
    content: str
    software: Optional[str] = None
    context: Optional[Dict[str, Any]] = None


@dataclass
class AIResponse:
    """AI response returned to the extension process."""

    success: bool
    content: Optional[str] = None
    suggestions: Optional[list] = None
    generatedInput: Optional[str] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class ProviderError(RuntimeError):
    """A provider response that is safe to surface after sanitization."""


class AIClient:
    """Client for OpenAI Responses API and Ollama operations."""

    def __init__(self) -> None:
        self.provider = os.environ.get("OPENQC_AI_PROVIDER", "ollama").strip().lower()
        self.model = os.environ.get("OPENQC_AI_MODEL", "llama2").strip()
        self.api_key = os.environ.get("OPENQC_AI_API_KEY", "")
        self.openai_url = os.environ.get(
            "OPENQC_AI_OPENAI_URL", "https://api.openai.com/v1"
        ).rstrip("/")
        self.ollama_url = os.environ.get(
            "OPENQC_AI_OLLAMA_URL", "http://localhost:11434"
        ).rstrip("/")
        self.max_tokens = _bounded_int("OPENQC_AI_MAX_TOKENS", 2048, 1, 131_072)
        self.max_output_chars = _bounded_int(
            "OPENQC_AI_MAX_OUTPUT_CHARS", DEFAULT_MAX_OUTPUT_CHARS, 256, 1_000_000
        )
        self.temperature = _bounded_float("OPENQC_AI_TEMPERATURE", 0.7, 0.0, 2.0)
        self.timeout = _bounded_float(
            "OPENQC_AI_TIMEOUT_SECONDS", DEFAULT_TIMEOUT_SECONDS, 0.1, 600.0
        )
        self.parser = ResponseParser()

    def check(self) -> AIResponse:
        """Check whether the configured provider is reachable."""
        if self.provider == "openai":
            if not self.api_key:
                return AIResponse(
                    success=False, error="OpenAI API key is not configured"
                )
            return self._request_status(
                f"{self.openai_url}/models",
                headers={"Authorization": f"Bearer {self.api_key}"},
                provider="OpenAI",
            )
        if self.provider == "ollama":
            return self._request_status(
                f"{self.ollama_url}/api/tags", provider="Ollama"
            )
        return AIResponse(success=False, error=f"Unknown AI provider: {self.provider}")

    def optimize(self, request: AIRequest) -> AIResponse:
        """Optimize an input file."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        response = self._call_llm(
            PromptTemplates.get_optimize_prompt(
                request.content, request.software, request.context
            )
        )
        if not response.success:
            return response
        parsed = self.parser.parse_optimization(response.content or "")
        return AIResponse(
            success=True,
            content=parsed.get("content"),
            suggestions=parsed.get("suggestions", []),
            metadata=response.metadata,
        )

    def generate(self, request: AIRequest) -> AIResponse:
        """Generate an input file."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        response = self._call_llm(
            PromptTemplates.get_generate_prompt(
                request.content, request.software, request.context
            )
        )
        if not response.success:
            return response
        return AIResponse(
            success=True,
            generatedInput=self.parser.parse_generated_input(response.content or ""),
            metadata=response.metadata,
        )

    def explain(self, request: AIRequest) -> AIResponse:
        """Explain input parameters."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        return self._call_llm(
            PromptTemplates.get_explain_prompt(request.content, request.software)
        )

    def debug(self, request: AIRequest) -> AIResponse:
        """Debug a failed calculation."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        output = request.context.get("output", "") if request.context else ""
        return self._call_llm(
            PromptTemplates.get_debug_prompt(request.content, output, request.software)
        )

    def _request_status(
        self, url: str, provider: str, headers: Optional[Dict[str, str]] = None
    ) -> AIResponse:
        try:
            self._http_json(url, headers=headers)
            return AIResponse(success=True, content=f"{provider} connection OK")
        except (
            Exception
        ) as exc:  # Boundary converts all network failures to protocol JSON.
            return AIResponse(success=False, error=self._safe_error(provider, exc))

    def _call_llm(self, prompt: str) -> AIResponse:
        if self.provider == "openai":
            return self._call_openai(prompt)
        if self.provider == "ollama":
            return self._call_ollama(prompt)
        return AIResponse(success=False, error=f"Unknown AI provider: {self.provider}")

    def _call_openai(self, prompt: str) -> AIResponse:
        if not self.api_key:
            return AIResponse(success=False, error="OpenAI API key is not configured")
        payload = {
            "model": self.model,
            "instructions": "You are a quantum chemistry expert.",
            "input": prompt,
            "max_output_tokens": self.max_tokens,
            "temperature": self.temperature,
        }
        try:
            data = self._http_json(
                f"{self.openai_url}/responses",
                payload,
                {"Authorization": f"Bearer {self.api_key}"},
            )
            content = _extract_openai_text(data)
            return self._bounded_success(content)
        except (
            Exception
        ) as exc:  # Boundary converts all provider failures to protocol JSON.
            return AIResponse(success=False, error=self._safe_error("OpenAI", exc))

    def _call_ollama(self, prompt: str) -> AIResponse:
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {
                "temperature": self.temperature,
                "num_predict": self.max_tokens,
            },
        }
        try:
            data = self._http_json(f"{self.ollama_url}/api/generate", payload)
            content = data.get("response")
            if not isinstance(content, str):
                raise ProviderError("provider response did not contain text")
            return self._bounded_success(content)
        except (
            Exception
        ) as exc:  # Boundary converts all provider failures to protocol JSON.
            return AIResponse(success=False, error=self._safe_error("Ollama", exc))

    def _bounded_success(self, content: str) -> AIResponse:
        truncated = len(content) > self.max_output_chars
        bounded = content[: self.max_output_chars]
        metadata = {"truncated": True} if truncated else None
        return AIResponse(success=True, content=bounded, metadata=metadata)

    def _http_json(
        self,
        url: str,
        payload: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        request_headers = {"Accept": "application/json", "User-Agent": "OpenQC/0.0.1"}
        request_headers.update(headers or {})
        body = None
        if payload is not None:
            body = json.dumps(payload).encode("utf-8")
            request_headers["Content-Type"] = "application/json"
        request = urllib.request.Request(url, data=body, headers=request_headers)
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                raw = response.read(self.max_output_chars * 4 + 1)
                if len(raw) > self.max_output_chars * 4:
                    raise ProviderError(
                        "provider response exceeded the safe size limit"
                    )
                result = json.loads(raw.decode("utf-8"))
                if not isinstance(result, dict):
                    raise ProviderError("provider returned invalid JSON")
                return result
        except urllib.error.HTTPError as exc:
            detail = ""
            try:
                body_data = json.loads(exc.read(4096).decode("utf-8", errors="replace"))
                detail = _provider_error_message(body_data)
            except (json.JSONDecodeError, UnicodeDecodeError):
                pass
            finally:
                exc.close()
            suffix = f": {detail}" if detail else ""
            raise ProviderError(f"HTTP {exc.code}{suffix}") from None
        except (urllib.error.URLError, socket.timeout, TimeoutError) as exc:
            reason = getattr(exc, "reason", exc)
            if (
                isinstance(reason, (socket.timeout, TimeoutError))
                or "timed out" in str(reason).lower()
            ):
                raise ProviderError(
                    f"request timed out after {self.timeout:g} seconds"
                ) from None
            raise ProviderError("provider is unavailable") from None
        except json.JSONDecodeError:
            raise ProviderError("provider returned invalid JSON") from None

    def _safe_error(self, provider: str, error: Exception) -> str:
        message = str(error) or error.__class__.__name__
        if self.api_key:
            message = message.replace(self.api_key, "[REDACTED]")
        message = re.sub(r"(?i)bearer\s+[^\s,;]+", "Bearer [REDACTED]", message)
        message = re.sub(r"\bsk-[A-Za-z0-9_-]{8,}\b", "[REDACTED]", message)
        return f"{provider} request failed: {message[:500]}"


def _extract_openai_text(data: Dict[str, Any]) -> str:
    output_text = data.get("output_text")
    if isinstance(output_text, str):
        return output_text
    parts = []
    for item in data.get("output", []):
        if not isinstance(item, dict):
            continue
        for content in item.get("content", []):
            if isinstance(content, dict) and isinstance(content.get("text"), str):
                parts.append(content["text"])
    if not parts:
        raise ProviderError("provider response did not contain output text")
    return "".join(parts)


def _provider_error_message(data: Any) -> str:
    if not isinstance(data, dict):
        return ""
    error = data.get("error", data)
    if isinstance(error, dict) and isinstance(error.get("message"), str):
        return error["message"][:300]
    if isinstance(error, str):
        return error[:300]
    return ""


def _bounded_int(name: str, default: int, minimum: int, maximum: int) -> int:
    try:
        value = int(os.environ.get(name, str(default)))
    except ValueError:
        value = default
    return max(minimum, min(maximum, value))


def _bounded_float(name: str, default: float, minimum: float, maximum: float) -> float:
    try:
        value = float(os.environ.get(name, str(default)))
    except ValueError:
        value = default
    return max(minimum, min(maximum, value))


def main() -> None:
    """CLI/module entry point."""
    parser = argparse.ArgumentParser(description="OpenQC AI Client")
    parser.add_argument(
        "command", choices=["check", "optimize", "generate", "explain", "debug"]
    )
    args = parser.parse_args()
    client = AIClient()

    if args.command == "check":
        print(json.dumps(asdict(client.check())))
        return

    try:
        data = json.loads(sys.stdin.read())
        if not isinstance(data, dict):
            raise ValueError("request must be a JSON object")
        request = AIRequest(
            type=str(data.get("type", "")),
            content=str(data.get("content", "")),
            software=data.get("software"),
            context=data.get("context"),
        )
        operation = getattr(client, args.command)
        response = operation(request)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        response = AIResponse(success=False, error=f"Invalid request: {str(exc)[:300]}")

    print(json.dumps(asdict(response)))


if __name__ == "__main__":
    main()
