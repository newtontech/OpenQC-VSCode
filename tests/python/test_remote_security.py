import importlib
import sys
import types
import unittest


class FakeAutoAddPolicy:
    pass


class FakeRejectPolicy:
    pass


class FakeSSHClient:
    def __init__(self):
        self.loaded_system_host_keys = False
        self.loaded_host_keys = []
        self.policy = None
        self.connected_with = None
        self.exec_calls = []

    def load_system_host_keys(self):
        self.loaded_system_host_keys = True

    def load_host_keys(self, path):
        self.loaded_host_keys.append(path)

    def set_missing_host_key_policy(self, policy):
        self.policy = policy

    def connect(self, **kwargs):
        self.connected_with = kwargs

    def exec_command(self, command, timeout=None):
        self.exec_calls.append((command, timeout))
        return None, FakeStream(b"ok"), FakeStream(b"")


class FakeChannel:
    def recv_exit_status(self):
        return 0


class FakeStream:
    def __init__(self, data):
        self.data = data
        self.channel = FakeChannel()

    def read(self, size=-1):
        if size < 0:
            return self.data
        return self.data[:size]


class FakeSSH:
    def __init__(self):
        self.calls = []

    def execute_args(self, args):
        self.calls.append(args)
        return 0, "Submitted batch job 12345", ""


def install_fake_paramiko():
    fake_paramiko = types.SimpleNamespace(
        SSHClient=FakeSSHClient,
        AutoAddPolicy=FakeAutoAddPolicy,
        RejectPolicy=FakeRejectPolicy,
    )
    sys.modules["paramiko"] = fake_paramiko


class RemoteSecurityTests(unittest.TestCase):
    def setUp(self):
        install_fake_paramiko()
        for module_name in [
            "server.ssh_handler",
            "server.slurm_interface",
            "server.file_transfer",
        ]:
            sys.modules.pop(module_name, None)

    def test_ssh_rejects_unknown_hosts_by_default(self):
        ssh_handler = importlib.import_module("server.ssh_handler")

        handler = ssh_handler.SSHHandler(
            "cluster.example",
            "alice",
            known_hosts_path="/tmp/known_hosts",
        )

        self.assertTrue(handler.connect())
        self.assertTrue(handler._client.loaded_system_host_keys)
        self.assertEqual(handler._client.loaded_host_keys, ["/tmp/known_hosts"])
        self.assertIsInstance(handler._client.policy, FakeRejectPolicy)
        self.assertEqual(handler._client.connected_with["timeout"], 10)
        self.assertEqual(handler._client.connected_with["banner_timeout"], 10)
        self.assertEqual(handler._client.connected_with["auth_timeout"], 10)

    def test_unknown_host_trust_must_be_explicit(self):
        ssh_handler = importlib.import_module("server.ssh_handler")

        handler = ssh_handler.SSHHandler(
            "cluster.example",
            "alice",
            allow_unknown_host=True,
        )

        self.assertTrue(handler.connect())
        self.assertIsInstance(handler._client.policy, FakeAutoAddPolicy)

    def test_execute_args_quotes_shell_arguments_and_sets_timeout(self):
        ssh_handler = importlib.import_module("server.ssh_handler")
        handler = ssh_handler.SSHHandler("cluster.example", "alice")
        self.assertTrue(handler.connect())

        result = handler.execute_args(["mkdir", "-p", "/scratch/alice/a dir; rm -rf /"])

        self.assertEqual(result, (0, "ok", ""))
        self.assertEqual(
            handler._client.exec_calls,
            [("mkdir -p '/scratch/alice/a dir; rm -rf /'", 300)]
        )

    def test_execute_enforces_output_limit(self):
        ssh_handler = importlib.import_module("server.ssh_handler")
        handler = ssh_handler.SSHHandler("cluster.example", "alice", max_output_bytes=1)

        with self.assertRaises(ssh_handler.RemoteCommandError):
            handler._read_limited(FakeStream(b"too much"), "stdout")

    def test_slurm_uses_structured_commands_and_rejects_injected_identifiers(self):
        slurm_interface = importlib.import_module("server.slurm_interface")
        ssh = FakeSSH()
        slurm = slurm_interface.SlurmInterface(ssh)

        job_id = slurm.submit_job("/scratch/alice/job script.sh", job_name="qc_job-1")
        cancelled = slurm.cancel_job("12345")

        self.assertEqual(job_id, "12345")
        self.assertTrue(cancelled)
        self.assertEqual(
            ssh.calls,
            [
                ["sbatch", "--job-name=qc_job-1", "/scratch/alice/job script.sh"],
                ["scancel", "12345"],
            ],
        )
        with self.assertRaises(ValueError):
            slurm.cancel_job("12345; rm -rf /")

    def test_file_transfer_creates_remote_directory_with_structured_command(self):
        file_transfer = importlib.import_module("server.file_transfer")

        class FakeSftp:
            def stat(self, path):
                raise FileNotFoundError()

            def put(self, local_path, remote_path, callback=None):
                self.put_call = (local_path, remote_path, callback)

        class FakeClient:
            def __init__(self):
                self.sftp = FakeSftp()

            def open_sftp(self):
                return self.sftp

        ssh = FakeSSH()
        ssh._client = FakeClient()
        transfer = file_transfer.FileTransfer(ssh)

        self.assertTrue(transfer.upload("local.inp", "/scratch/alice/a dir/job.inp"))
        self.assertEqual(ssh.calls, [["mkdir", "-p", "/scratch/alice/a dir"]])


if __name__ == "__main__":
    unittest.main()
