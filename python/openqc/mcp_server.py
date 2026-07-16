"""Dependency-free JSON-RPC 2.0 MCP server for OpenQC scientific bridges."""

from __future__ import annotations

import json
import re
import sys
from typing import Any, Callable, Dict, Optional

TOOLS = [
    {
        "name": "openqc.parseStructure",
        "description": "Parse a molecular or crystal structure into OpenQCStructure JSON",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "File path"},
                "format": {"type": "string", "description": "Optional format hint"},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.parseCalculationOutput",
        "description": "Parse calculation output with cclib and a native fallback",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Output file path"},
                "software": {"type": "string", "description": "Optional software hint"},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.checkPythonBackend",
        "description": "Check the Scientific Python backend and optional dependencies",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "openqc.generateSupercell",
        "description": "Generate a bounded supercell preview",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string"},
                "na": {"type": "integer", "minimum": 1, "maximum": 100, "default": 2},
                "nb": {"type": "integer", "minimum": 1, "maximum": 100, "default": 2},
                "nc": {"type": "integer", "minimum": 1, "maximum": 100, "default": 2},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.summarizeDataset",
        "description": "Summarize an atomistic dataset when dpdata is installed",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string"},
                "format": {"type": "string", "default": "auto"},
            },
            "required": ["path"],
        },
    },
]

_TOOL_SCHEMAS = {tool["name"]: tool["inputSchema"] for tool in TOOLS}


class RpcError(Exception):
    """JSON-RPC error with a stable code and safe message."""

    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


def _parse_structure(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from openqc.bridge.structure_bridge import parse_file

    return parse_file(arguments["path"], arguments.get("format"))


def _parse_output(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from openqc.bridge.output_bridge import _parse_output_native, _parse_with_cclib

    path = arguments["path"]
    software = arguments.get("software", "auto")
    try:
        result = _parse_with_cclib(path, software)
        if result is None:
            raise ValueError("cclib returned no output")
        result["cclibAvailable"] = True
        return result
    except Exception as cclib_error:
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            result = _parse_output_native(handle.read(), software)
        result["sourceFile"] = path
        result["cclibAvailable"] = False
        result.setdefault("warnings", []).append(
            f"cclib parser unavailable or failed: {_safe_message(cclib_error)}"
        )
        return result


def _check_backend(_arguments: Dict[str, Any]) -> Dict[str, Any]:
    from openqc.bridge.check_backend import check_backend

    return check_backend()


def _generate_supercell(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from openqc.bridge.structure_bridge import _generate_supercell, parse_file

    structure = parse_file(arguments["path"])
    return _generate_supercell(
        structure,
        arguments.get("na", 2),
        arguments.get("nb", 2),
        arguments.get("nc", 2),
    )


def _summarize_dataset(arguments: Dict[str, Any]) -> Dict[str, Any]:
    from openqc.bridge.dpdata_bridge import _summarize_dataset

    return _summarize_dataset(arguments["path"], arguments.get("format", "auto"))


_HANDLERS: Dict[str, Callable[[Dict[str, Any]], Dict[str, Any]]] = {
    "openqc.parseStructure": _parse_structure,
    "openqc.parseCalculationOutput": _parse_output,
    "openqc.checkPythonBackend": _check_backend,
    "openqc.generateSupercell": _generate_supercell,
    "openqc.summarizeDataset": _summarize_dataset,
}


def _validate_request(request: Any) -> Dict[str, Any]:
    if not isinstance(request, dict) or request.get("jsonrpc") != "2.0":
        raise RpcError(-32600, "Invalid Request")
    if not isinstance(request.get("method"), str):
        raise RpcError(-32600, "Invalid Request")
    params = request.get("params", {})
    if not isinstance(params, dict):
        raise RpcError(-32602, "Invalid params: params must be an object")
    return request


def _validate_arguments(tool_name: str, arguments: Any) -> Dict[str, Any]:
    if not isinstance(arguments, dict):
        raise RpcError(-32602, "Invalid params: arguments must be an object")
    schema = _TOOL_SCHEMAS[tool_name]
    result = dict(arguments)
    for name in schema.get("required", []):
        if name not in result:
            raise RpcError(
                -32602, f"Invalid params: missing required argument '{name}'"
            )
    for name, definition in schema.get("properties", {}).items():
        if name not in result and "default" in definition:
            result[name] = definition["default"]
        if name not in result:
            continue
        value = result[name]
        expected = definition.get("type")
        if expected == "string" and (not isinstance(value, str) or not value.strip()):
            raise RpcError(
                -32602, f"Invalid params: '{name}' must be a non-empty string"
            )
        if expected == "integer" and (
            not isinstance(value, int) or isinstance(value, bool)
        ):
            raise RpcError(-32602, f"Invalid params: '{name}' must be an integer")
        minimum = definition.get("minimum")
        maximum = definition.get("maximum")
        if minimum is not None and value < minimum:
            raise RpcError(
                -32602, f"Invalid params: '{name}' is below minimum {minimum}"
            )
        if maximum is not None and value > maximum:
            raise RpcError(
                -32602, f"Invalid params: '{name}' exceeds maximum {maximum}"
            )
    return result


def _result(request_id: Any, payload: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "jsonrpc": "2.0",
        "id": request_id,
        "result": {
            "content": [
                {"type": "text", "text": json.dumps(payload, ensure_ascii=False)}
            ]
        },
    }


def _error(request_id: Any, code: int, message: str) -> Dict[str, Any]:
    return {
        "jsonrpc": "2.0",
        "id": request_id,
        "error": {"code": code, "message": message},
    }


def _safe_message(error: Exception) -> str:
    message = str(error) or error.__class__.__name__
    message = re.sub(r"(?i)bearer\s+[^\s,;]+", "Bearer [REDACTED]", message)
    message = re.sub(r"\bsk-[A-Za-z0-9_-]{8,}\b", "[REDACTED]", message)
    return message[:500]


def handle_request(request: Any) -> Optional[Dict[str, Any]]:
    """Validate and dispatch one JSON-RPC request without terminating the server."""
    request_id = request.get("id") if isinstance(request, dict) else None
    is_notification = isinstance(request, dict) and "id" not in request
    try:
        validated = _validate_request(request)
        method = validated["method"]
        params = validated.get("params", {})

        if method == "notifications/initialized":
            return None
        if method == "initialize":
            response = {
                "jsonrpc": "2.0",
                "id": request_id,
                "result": {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {"tools": {}},
                    "serverInfo": {"name": "openqc-mcp", "version": "0.1.0"},
                },
            }
        elif method == "tools/list":
            response = {
                "jsonrpc": "2.0",
                "id": request_id,
                "result": {"tools": TOOLS},
            }
        elif method == "tools/call":
            tool_name = params.get("name")
            if not isinstance(tool_name, str):
                raise RpcError(-32602, "Invalid params: tool name must be a string")
            handler = _HANDLERS.get(tool_name)
            if handler is None:
                raise RpcError(-32601, f"Unknown tool: {tool_name}")
            arguments = _validate_arguments(tool_name, params.get("arguments", {}))
            try:
                response = _result(request_id, handler(arguments))
            except ImportError as exc:
                raise RpcError(
                    -32001, f"Optional dependency unavailable: {_safe_message(exc)}"
                ) from None
            except (
                FileNotFoundError,
                IsADirectoryError,
                PermissionError,
                ValueError,
            ) as exc:
                raise RpcError(
                    -32002, f"Tool execution failed: {_safe_message(exc)}"
                ) from None
            except Exception as exc:
                raise RpcError(
                    -32603, f"Internal error: {_safe_message(exc)}"
                ) from None
        else:
            raise RpcError(-32601, f"Unknown method: {method}")

        return None if is_notification else response
    except RpcError as exc:
        return None if is_notification else _error(request_id, exc.code, exc.message)


def main() -> None:
    """Run newline-delimited JSON-RPC over stdin/stdout."""
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
        except json.JSONDecodeError:
            response = _error(None, -32700, "Parse error")
        else:
            response = handle_request(request)
        if response is not None:
            sys.stdout.write(json.dumps(response, ensure_ascii=False) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
