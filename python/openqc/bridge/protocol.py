"""
JSON stdin/stdout protocol for OpenQC Python bridge.

Protocol:
  - Read JSON request from stdin.
  - Write JSON response to stdout.
  - Errors go to stderr as JSON.

Request:
  {
    "command": "check_backend" | "parse_structure" | ...,
    "args": { ... }
  }

Response:
  {
    "success": true/false,
    "data": { ... },       # on success
    "error": { ... }       # on failure
  }
"""

import json
import sys
from typing import Any, Dict, Optional


def read_request() -> Dict[str, Any]:
    """Read a JSON request from stdin."""
    raw = sys.stdin.read()
    try:
        request = json.loads(raw)
    except json.JSONDecodeError as e:
        write_error(f"Invalid JSON request: {e}")
        sys.exit(1)

    if not isinstance(request, dict):
        write_error("Request must be a JSON object")
        sys.exit(1)

    return request


def write_response(data: Any) -> None:
    """Write a success response to stdout."""
    response = {"success": True, "data": data}
    sys.stdout.write(json.dumps(response) + "\n")
    sys.stdout.flush()


def write_error(message: str, details: Optional[Dict[str, Any]] = None) -> None:
    """Write an error response to stderr."""
    error = {"message": message}
    if details:
        error.update(details)
    response = {"success": False, "error": error}
    sys.stderr.write(json.dumps(response) + "\n")
    sys.stderr.flush()


def get_command(request: Dict[str, Any]) -> str:
    """Extract the command from a request."""
    command = request.get("command")
    if not isinstance(command, str) or not command:
        write_error("Missing or invalid 'command' field")
        sys.exit(1)
    return command


def get_args(request: Dict[str, Any]) -> Dict[str, Any]:
    """Extract args from a request, defaulting to empty dict."""
    args = request.get("args", {})
    if not isinstance(args, dict):
        write_error("'args' must be a JSON object")
        sys.exit(1)
    return args
