"""
AI Client for LLM Integration

Supports OpenAI API and Ollama (local models).
"""

import os
import sys
import json
import argparse
from typing import Dict, Optional, Any
from dataclasses import dataclass
from .prompts import PromptTemplates
from .parser import ResponseParser


@dataclass
class AIRequest:
    """AI request data structure."""
    type: str
    content: str
    software: Optional[str] = None
    context: Optional[Dict[str, Any]] = None


@dataclass
class AIResponse:
    """AI response data structure."""
    success: bool
    content: Optional[str] = None
    suggestions: Optional[list] = None
    generatedInput: Optional[str] = None
    error: Optional[str] = None
    metadata: Optional[Dict[str, Any]] = None


class AIClient:
    """Client for AI/LLM operations."""

    def __init__(self):
        self.provider = os.environ.get('OPENQC_AI_PROVIDER', 'ollama')
        self.model = os.environ.get('OPENQC_AI_MODEL', 'llama2')
        self.api_key = os.environ.get('OPENQC_AI_API_KEY', '')
        self.ollama_url = os.environ.get('OPENQC_AI_OLLAMA_URL', 'http://localhost:11434')
        self.max_tokens = int(os.environ.get('OPENQC_AI_MAX_TOKENS', '2048'))
        self.temperature = float(os.environ.get('OPENQC_AI_TEMPERATURE', '0.7'))
        self.parser = ResponseParser()

    def check(self) -> AIResponse:
        """Check if AI backend is available."""
        try:
            if self.provider == 'openai':
                return self._check_openai()
            elif self.provider == 'ollama':
                return self._check_ollama()
            else:
                return AIResponse(success=False, error=f"Unknown provider: {self.provider}")
        except Exception as e:
            return AIResponse(success=False, error=str(e))

    def _check_openai(self) -> AIResponse:
        """Check OpenAI availability."""
        if not self.api_key:
            return AIResponse(success=False, error="OpenAI API key not configured")
        
        try:
            import openai
            openai.api_key = self.api_key
            # Test with a simple request
            models = openai.Model.list()
            return AIResponse(
                success=True,
                content=f"OpenAI connection OK. Available models: {len(models.data)}"
            )
        except ImportError:
            return AIResponse(success=False, error="OpenAI package not installed. Run: pip install openai")
        except Exception as e:
            return AIResponse(success=False, error=f"OpenAI connection failed: {e}")

    def _check_ollama(self) -> AIResponse:
        """Check Ollama availability."""
        try:
            import requests
            response = requests.get(f"{self.ollama_url}/api/tags", timeout=5)
            if response.status_code == 200:
                data = response.json()
                models = data.get('models', [])
                return AIResponse(
                    success=True,
                    content=f"Ollama connection OK. Available models: {len(models)}"
                )
            else:
                return AIResponse(success=False, error=f"Ollama returned status {response.status_code}")
        except ImportError:
            return AIResponse(success=False, error="requests package not installed")
        except Exception as e:
            return AIResponse(success=False, error=f"Ollama connection failed: {e}")

    def optimize(self, request: AIRequest) -> AIResponse:
        """Optimize input file."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        
        prompt = PromptTemplates.get_optimize_prompt(
            request.content,
            request.software,
            request.context
        )
        
        response = self._call_llm(prompt)
        if not response.success:
            return response
        
        # Parse suggestions from response
        parsed = self.parser.parse_optimization(response.content or '')
        return AIResponse(
            success=True,
            content=parsed.get('content'),
            suggestions=parsed.get('suggestions', []),
            metadata={'raw_response': response.content}
        )

    def generate(self, request: AIRequest) -> AIResponse:
        """Generate input file."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        
        prompt = PromptTemplates.get_generate_prompt(
            request.content,
            request.software,
            request.context
        )
        
        response = self._call_llm(prompt)
        if not response.success:
            return response
        
        # Extract input file from response
        generated = self.parser.parse_generated_input(response.content or '')
        return AIResponse(
            success=True,
            generatedInput=generated,
            metadata={'raw_response': response.content}
        )

    def explain(self, request: AIRequest) -> AIResponse:
        """Explain parameters."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        
        prompt = PromptTemplates.get_explain_prompt(
            request.content,
            request.software
        )
        
        return self._call_llm(prompt)

    def debug(self, request: AIRequest) -> AIResponse:
        """Debug calculation."""
        if not request.software:
            return AIResponse(success=False, error="Software type not specified")
        
        output = request.context.get('output', '') if request.context else ''
        
        prompt = PromptTemplates.get_debug_prompt(
            request.content,
            output,
            request.software
        )
        
        return self._call_llm(prompt)

    def _call_llm(self, prompt: str) -> AIResponse:
        """Call the LLM with the given prompt."""
        try:
            if self.provider == 'openai':
                return self._call_openai(prompt)
            elif self.provider == 'ollama':
                return self._call_ollama(prompt)
            else:
                return AIResponse(success=False, error=f"Unknown provider: {self.provider}")
        except Exception as e:
            return AIResponse(success=False, error=f"LLM call failed: {e}")

    def _call_openai(self, prompt: str) -> AIResponse:
        """Call OpenAI API."""
        try:
            import openai
            openai.api_key = self.api_key
            
            response = openai.ChatCompletion.create(
                model=self.model,
                messages=[
                    {"role": "system", "content": "You are a quantum chemistry expert."},
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )
            
            content = response.choices[0].message.content
            return AIResponse(success=True, content=content)
        except ImportError:
            return AIResponse(
                success=False,
                error="OpenAI package not installed. Run: pip install openai"
            )
        except Exception as e:
            return AIResponse(success=False, error=f"OpenAI API error: {e}")

    def _call_ollama(self, prompt: str) -> AIResponse:
        """Call Ollama API."""
        try:
            import requests
            
            response = requests.post(
                f"{self.ollama_url}/api/generate",
                json={
                    "model": self.model,
                    "prompt": prompt,
                    "stream": False,
                    "options": {
                        "temperature": self.temperature,
                        "num_predict": self.max_tokens
                    }
                },
                timeout=120
            )
            
            if response.status_code == 200:
                data = response.json()
                content = data.get('response', '')
                return AIResponse(success=True, content=content)
            else:
                return AIResponse(
                    success=False,
                    error=f"Ollama returned status {response.status_code}: {response.text}"
                )
        except ImportError:
            return AIResponse(
                success=False,
                error="requests package not installed. Run: pip install requests"
            )
        except Exception as e:
            return AIResponse(success=False, error=f"Ollama API error: {e}")


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description='OpenQC AI Client')
    parser.add_argument('command', choices=['check', 'optimize', 'generate', 'explain', 'debug'])
    args = parser.parse_args()
    
    client = AIClient()
    
    if args.command == 'check':
        response = client.check()
        print(json.dumps({
            'success': response.success,
            'content': response.content,
            'error': response.error
        }))
        return
    
    # Read request from stdin
    request_data = sys.stdin.read()
    if not request_data:
        print(json.dumps({'success': False, 'error': 'No input provided'}))
        return
    
    try:
        data = json.loads(request_data)
        request = AIRequest(
            type=data.get('type', ''),
            content=data.get('content', ''),
            software=data.get('software'),
            context=data.get('context')
        )
    except json.JSONDecodeError as e:
        print(json.dumps({'success': False, 'error': f'Invalid JSON: {e}'}))
        return
    
    # Execute command
    if args.command == 'optimize':
        response = client.optimize(request)
    elif args.command == 'generate':
        response = client.generate(request)
    elif args.command == 'explain':
        response = client.explain(request)
    elif args.command == 'debug':
        response = client.debug(request)
    else:
        response = AIResponse(success=False, error=f"Unknown command: {args.command}")
    
    # Output response
    print(json.dumps({
        'success': response.success,
        'content': response.content,
        'suggestions': response.suggestions,
        'generatedInput': response.generatedInput,
        'error': response.error,
        'metadata': response.metadata
    }))


if __name__ == '__main__':
    main()
