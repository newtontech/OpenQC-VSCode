import json
import unittest
from unittest.mock import patch

from openqc.mcp_server import _HANDLERS, TOOLS, handle_request


def request(request_id, method, params=None):
    value = {"jsonrpc": "2.0", "id": request_id, "method": method}
    if params is not None:
        value["params"] = params
    return value


class McpContractTest(unittest.TestCase):
    def test_tools_list_and_every_declared_tool_are_callable(self):
        listed = handle_request(request("list-id", "tools/list"))
        tools = listed["result"]["tools"]
        self.assertEqual(tools, TOOLS)
        self.assertEqual(len(tools), 5)

        fake_handlers = {
            tool["name"]: (lambda args, name=tool["name"]: {"tool": name, "args": args})
            for tool in tools
        }
        arguments = {
            "openqc.parseStructure": {"path": "structure.xyz"},
            "openqc.parseCalculationOutput": {"path": "calculation.out"},
            "openqc.checkPythonBackend": {},
            "openqc.generateSupercell": {"path": "POSCAR"},
            "openqc.summarizeDataset": {"path": "dataset"},
        }
        with patch.dict(_HANDLERS, fake_handlers, clear=True):
            for index, tool in enumerate(tools):
                request_id = f"tool-{index}"
                response = handle_request(
                    request(
                        request_id,
                        "tools/call",
                        {"name": tool["name"], "arguments": arguments[tool["name"]]},
                    )
                )
                self.assertEqual(response["id"], request_id)
                payload = json.loads(response["result"]["content"][0]["text"])
                self.assertEqual(payload["tool"], tool["name"])

        self.assertEqual(
            json.loads(response["result"]["content"][0]["text"])["args"]["format"],
            "auto",
        )

    def test_invalid_jsonrpc_and_params_have_stable_codes_and_ids(self):
        invalid = handle_request({"jsonrpc": "1.0", "id": 0, "method": "tools/list"})
        self.assertEqual(invalid["id"], 0)
        self.assertEqual(invalid["error"]["code"], -32600)

        missing = handle_request(
            request(
                42,
                "tools/call",
                {"name": "openqc.parseStructure", "arguments": {}},
            )
        )
        self.assertEqual(missing["id"], 42)
        self.assertEqual(missing["error"]["code"], -32602)

        wrong_type = handle_request(
            request(
                "same-id",
                "tools/call",
                {
                    "name": "openqc.generateSupercell",
                    "arguments": {"path": "POSCAR", "na": 1.5},
                },
            )
        )
        self.assertEqual(wrong_type["id"], "same-id")
        self.assertEqual(wrong_type["error"]["code"], -32602)

    def test_optional_dependency_and_internal_errors_are_normalized(self):
        def missing(_arguments):
            raise ImportError("install dpdata")

        def broken(_arguments):
            raise RuntimeError("Bearer example-credential")

        call = lambda request_id: request(
            request_id,
            "tools/call",
            {"name": "openqc.summarizeDataset", "arguments": {"path": "dataset"}},
        )

        with patch.dict(_HANDLERS, {"openqc.summarizeDataset": missing}):
            dependency = handle_request(call("dependency-id"))
        self.assertEqual(dependency["id"], "dependency-id")
        self.assertEqual(dependency["error"]["code"], -32001)

        with patch.dict(_HANDLERS, {"openqc.summarizeDataset": broken}):
            internal = handle_request(call("internal-id"))
        self.assertEqual(internal["id"], "internal-id")
        self.assertEqual(internal["error"]["code"], -32603)
        self.assertNotIn("example-credential", internal["error"]["message"])

    def test_notifications_do_not_receive_responses(self):
        self.assertIsNone(
            handle_request({"jsonrpc": "2.0", "method": "notifications/initialized"})
        )
        self.assertIsNone(handle_request({"jsonrpc": "2.0", "method": "unknown"}))


if __name__ == "__main__":
    unittest.main()
