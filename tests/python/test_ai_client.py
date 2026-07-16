import json
import os
import socketserver
import subprocess
import sys
import threading
import time
import unittest
from http.server import BaseHTTPRequestHandler
from pathlib import Path
from unittest.mock import patch

from openqc.ai.client import AIClient, AIRequest

PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


class _ThreadingServer(socketserver.ThreadingMixIn, socketserver.TCPServer):
    allow_reuse_address = True
    daemon_threads = True


class _ProviderHandler(BaseHTTPRequestHandler):
    requests = []
    delay = 0.0
    status = 200
    response = {"output_text": "fake response"}

    def do_POST(self):
        if self.delay:
            time.sleep(self.delay)
        length = int(self.headers.get("Content-Length", "0"))
        body = json.loads(self.rfile.read(length))
        self.__class__.requests.append((self.path, dict(self.headers), body))
        payload = json.dumps(self.response).encode()
        self.send_response(self.status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        try:
            self.wfile.write(payload)
        except BrokenPipeError:
            pass

    def log_message(self, _format, *_args):
        return


class AIClientFakeHttpTest(unittest.TestCase):
    def setUp(self):
        _ProviderHandler.requests = []
        _ProviderHandler.delay = 0.0
        _ProviderHandler.status = 200
        _ProviderHandler.response = {"output_text": "fake response"}
        self.server = _ThreadingServer(("127.0.0.1", 0), _ProviderHandler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.base_url = f"http://127.0.0.1:{self.server.server_address[1]}"

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join(timeout=2)

    def _environment(self, **overrides):
        values = {
            "OPENQC_AI_PROVIDER": "openai",
            "OPENQC_AI_MODEL": "gpt-test",
            "OPENQC_AI_API_KEY": "test-only-secret",
            "OPENQC_AI_OPENAI_URL": self.base_url,
            "OPENQC_AI_TIMEOUT_SECONDS": "1",
            "OPENQC_AI_MAX_OUTPUT_CHARS": "1000",
        }
        values.update(overrides)
        return patch.dict(os.environ, values, clear=False)

    def test_openai_uses_responses_api_and_bearer_auth(self):
        with self._environment():
            response = AIClient().explain(AIRequest("explain", "ENCUT=500", "vasp"))

        self.assertTrue(response.success)
        self.assertEqual(response.content, "fake response")
        path, headers, body = _ProviderHandler.requests[0]
        self.assertEqual(path, "/responses")
        self.assertEqual(headers["Authorization"], "Bearer test-only-secret")
        self.assertEqual(body["model"], "gpt-test")
        self.assertEqual(body["max_output_tokens"], 2048)
        self.assertIn("ENCUT=500", body["input"])

    def test_module_entrypoint_succeeds_against_fake_http_server(self):
        env = os.environ.copy()
        env.update(
            {
                "PYTHONPATH": str(PYTHON_ROOT),
                "OPENQC_AI_PROVIDER": "openai",
                "OPENQC_AI_MODEL": "gpt-test",
                "OPENQC_AI_API_KEY": "test-only-secret",
                "OPENQC_AI_OPENAI_URL": self.base_url,
            }
        )
        completed = subprocess.run(
            [sys.executable, "-m", "openqc.ai.client", "explain"],
            input=json.dumps(
                {"type": "explain", "content": "input", "software": "orca"}
            ),
            text=True,
            capture_output=True,
            env=env,
            cwd=PYTHON_ROOT,
            timeout=5,
            check=True,
        )

        payload = json.loads(completed.stdout)
        self.assertTrue(payload["success"])
        self.assertEqual(payload["content"], "fake response")
        self.assertNotIn("test-only-secret", completed.stdout + completed.stderr)
        self.assertNotIn("RuntimeWarning", completed.stderr)

    def test_package_exports_remain_available_with_lazy_client_import(self):
        from openqc.ai import AIClient as ExportedClient

        self.assertIs(ExportedClient, AIClient)

    def test_missing_openai_credential_fails_without_http(self):
        with self._environment(OPENQC_AI_API_KEY=""):
            response = AIClient().explain(AIRequest("explain", "input", "vasp"))

        self.assertFalse(response.success)
        self.assertEqual(response.error, "OpenAI API key is not configured")
        self.assertEqual(_ProviderHandler.requests, [])

    def test_provider_error_is_redacted(self):
        _ProviderHandler.status = 401
        _ProviderHandler.response = {
            "error": {"message": "bad Bearer test-only-secret"}
        }
        with self._environment():
            response = AIClient().explain(AIRequest("explain", "input", "vasp"))

        self.assertFalse(response.success)
        self.assertIn("HTTP 401", response.error)
        self.assertNotIn("test-only-secret", response.error)

    def test_timeout_is_bounded(self):
        _ProviderHandler.delay = 0.2
        with self._environment(OPENQC_AI_TIMEOUT_SECONDS="0.05"):
            response = AIClient().explain(AIRequest("explain", "input", "vasp"))

        self.assertFalse(response.success)
        self.assertIn("timed out", response.error)

    def test_output_is_truncated_to_configured_limit(self):
        _ProviderHandler.response = {"output_text": "x" * 1000}
        with self._environment(OPENQC_AI_MAX_OUTPUT_CHARS="256"):
            response = AIClient().explain(AIRequest("explain", "input", "vasp"))

        self.assertTrue(response.success)
        self.assertEqual(len(response.content), 256)
        self.assertEqual(response.metadata, {"truncated": True})

    def test_ollama_generate_api(self):
        _ProviderHandler.response = {"response": "local response"}
        with self._environment(
            OPENQC_AI_PROVIDER="ollama", OPENQC_AI_OLLAMA_URL=self.base_url
        ):
            response = AIClient().explain(AIRequest("explain", "input", "vasp"))

        self.assertTrue(response.success)
        self.assertEqual(response.content, "local response")
        self.assertEqual(_ProviderHandler.requests[0][0], "/api/generate")


if __name__ == "__main__":
    unittest.main()
