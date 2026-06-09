"""
OpenQC MCP Server - Exposes OpenQC operations as MCP tools.

This is a minimal MCP server that lists available tools and can be
launched by MCP-compatible clients (Claude, Codex, Cursor, etc.).

Usage:
  python -m openqc.mcp_server

Protocol: JSON-RPC 2.0 over stdin/stdout (MCP standard transport)
"""

import json
import sys


TOOLS = [
    {
        "name": "openqc.parseStructure",
        "description": "Parse a molecular/crystal structure file into OpenQCStructure JSON",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "File path"},
                "format": {"type": "string", "description": "Format hint (auto, xyz, poscar, cif)"},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.parseCalculationOutput",
        "description": "Parse a quantum chemistry output file using cclib",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Output file path"},
                "software": {"type": "string", "description": "Software hint"},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.checkPythonBackend",
        "description": "Check Scientific Python backend status",
        "inputSchema": {"type": "object", "properties": {}},
    },
    {
        "name": "openqc.generateSupercell",
        "description": "Generate a supercell preview",
        "inputSchema": {
            "type": "object",
            "properties": {
                "path": {"type": "string"},
                "na": {"type": "number", "default": 2},
                "nb": {"type": "number", "default": 2},
                "nc": {"type": "number", "default": 2},
            },
            "required": ["path"],
        },
    },
    {
        "name": "openqc.summarizeDataset",
        "description": "Summarize an atomistic dataset (dpdata)",
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


def handle_request(request: dict) -> dict:
    """Handle a JSON-RPC request."""
    method = request.get("method", "")
    req_id = request.get("id")
    params = request.get("params", {})

    if method == "initialize":
        return {
            "jsonrpc": "2.0",
            "id": req_id,
            "result": {
                "protocolVersion": "2024-11-05",
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "openqc-mcp", "version": "0.1.0"},
            },
        }

    elif method == "tools/list":
        return {"jsonrpc": "2.0", "id": req_id, "result": {"tools": TOOLS}}

    elif method == "tools/call":
        tool_name = params.get("name", "")
        arguments = params.get("arguments", {})

        # Route to appropriate bridge
        if tool_name == "openqc.checkPythonBackend":
            from openqc.bridge.check_backend import check_backend
            result = check_backend()
            return {"jsonrpc": "2.0", "id": req_id, "result": {"content": [{"type": "text", "text": json.dumps(result)}]}}

        elif tool_name == "openqc.parseStructure":
            from openqc.bridge.structure_bridge import parse_file
            result = parse_file(arguments["path"], arguments.get("format"))
            return {"jsonrpc": "2.0", "id": req_id, "result": {"content": [{"type": "text", "text": json.dumps(result)}]}}

        elif tool_name == "openqc.parseCalculationOutput":
            from openqc.bridge.output_bridge import _parse_with_cclib
            try:
                result = _parse_with_cclib(arguments["path"], arguments.get("software"))
                return {"jsonrpc": "2.0", "id": req_id, "result": {"content": [{"type": "text", "text": json.dumps(result)}]}}
            except ImportError as e:
                return {"jsonrpc": "2.0", "id": req_id, "result": {"content": [{"type": "text", "text": json.dumps({"error": str(e)})}]}}

        else:
            return {
                "jsonrpc": "2.0",
                "id": req_id,
                "error": {"code": -32601, "message": f"Unknown tool: {tool_name}"},
            }

    elif method == "notifications/initialized":
        return None  # No response for notifications

    else:
        return {
            "jsonrpc": "2.0",
            "id": req_id,
            "error": {"code": -32601, "message": f"Unknown method: {method}"},
        }


def main():
    """Run the MCP server."""
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue

        try:
            request = json.loads(line)
            response = handle_request(request)
            if response is not None:
                sys.stdout.write(json.dumps(response) + "\n")
                sys.stdout.flush()
        except json.JSONDecodeError as e:
            error_response = {
                "jsonrpc": "2.0",
                "id": None,
                "error": {"code": -32700, "message": f"Parse error: {e}"},
            }
            sys.stdout.write(json.dumps(error_response) + "\n")
            sys.stdout.flush()
        except Exception as e:
            error_response = {
                "jsonrpc": "2.0",
                "id": None,
                "error": {"code": -32603, "message": f"Internal error: {e}"},
            }
            sys.stdout.write(json.dumps(error_response) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()
