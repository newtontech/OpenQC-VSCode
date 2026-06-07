"""
SSH connection management for remote servers.
"""

import paramiko
import shlex
from typing import Optional, Sequence


class RemoteCommandError(RuntimeError):
    """Raised when a remote command cannot be executed safely."""


class SSHHandler:
    """
    SSH connection handler for remote compute servers.
    """
    
    def __init__(
        self,
        host: str,
        username: str,
        port: int = 22,
        key_path: Optional[str] = None,
        password: Optional[str] = None,
        known_hosts_path: Optional[str] = None,
        allow_unknown_host: bool = False,
        connect_timeout: int = 10,
        command_timeout: int = 300,
        max_output_bytes: int = 1024 * 1024
    ):
        """
        Initialize SSH handler.
        
        Args:
            host: Server hostname or IP
            username: SSH username
            port: SSH port (default 22)
            key_path: Path to SSH private key
            password: Password (if not using key)
            known_hosts_path: Optional known_hosts file to trust
            allow_unknown_host: Explicitly trust and save unknown host keys
            connect_timeout: Connection timeout in seconds
            command_timeout: Remote command timeout in seconds
            max_output_bytes: Maximum bytes read from stdout or stderr
        """
        self.host = host
        self.username = username
        self.port = port
        self.key_path = key_path
        self.password = password
        self.known_hosts_path = known_hosts_path
        self.allow_unknown_host = allow_unknown_host
        self.connect_timeout = connect_timeout
        self.command_timeout = command_timeout
        self.max_output_bytes = max_output_bytes
        
        self._client: Optional[paramiko.SSHClient] = None
    
    def connect(self) -> bool:
        """Establish SSH connection."""
        try:
            self._client = paramiko.SSHClient()
            self._client.load_system_host_keys()
            if self.known_hosts_path:
                self._client.load_host_keys(self.known_hosts_path)
            if self.allow_unknown_host:
                self._client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            else:
                self._client.set_missing_host_key_policy(paramiko.RejectPolicy())
            
            connect_kwargs = {
                'hostname': self.host,
                'port': self.port,
                'username': self.username,
                'timeout': self.connect_timeout,
                'banner_timeout': self.connect_timeout,
                'auth_timeout': self.connect_timeout,
                'look_for_keys': self.key_path is None and self.password is None,
            }
            
            if self.key_path:
                connect_kwargs['key_filename'] = self.key_path
            if self.password:
                connect_kwargs['password'] = self.password
            
            self._client.connect(**connect_kwargs)
            return True
            
        except Exception as e:
            print(f"SSH connection failed: {e}")
            return False
    
    def disconnect(self) -> None:
        """Close SSH connection."""
        if self._client:
            self._client.close()
            self._client = None
    
    def execute(self, command: str, *, timeout: Optional[int] = None) -> tuple[int, str, str]:
        """
        Execute command on remote server.
        
        Prefer execute_args() for command execution with untrusted arguments.
        
        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        if not self._client:
            raise RuntimeError("Not connected to server")
        if not isinstance(command, str) or not command.strip():
            raise ValueError("command must be a non-empty string")
        
        stdin, stdout, stderr = self._client.exec_command(
            command,
            timeout=timeout or self.command_timeout
        )
        exit_code = stdout.channel.recv_exit_status()
        
        return (
            exit_code,
            self._read_limited(stdout, "stdout"),
            self._read_limited(stderr, "stderr")
        )

    def execute_args(self, argv: Sequence[str], *, timeout: Optional[int] = None) -> tuple[int, str, str]:
        """
        Execute a command from an argv-style sequence.

        Arguments are shell-quoted before being sent to Paramiko because SSH
        servers execute command strings through the user's login shell.
        """
        return self.execute(self.build_command(argv), timeout=timeout)

    @staticmethod
    def build_command(argv: Sequence[str]) -> str:
        """Build a safely quoted remote shell command from argv."""
        if not argv:
            raise ValueError("argv must contain at least one argument")
        if any(not isinstance(arg, str) or arg == "" for arg in argv):
            raise ValueError("argv entries must be non-empty strings")
        return " ".join(shlex.quote(arg) for arg in argv)

    def _read_limited(self, stream, stream_name: str) -> str:
        data = stream.read(self.max_output_bytes + 1)
        if len(data) > self.max_output_bytes:
            raise RemoteCommandError(
                f"Remote command {stream_name} exceeded {self.max_output_bytes} bytes"
            )
        return data.decode('utf-8', errors='replace')
    
    @property
    def is_connected(self) -> bool:
        """Check if connection is active."""
        return self._client is not None and self._client.get_transport() is not None
    
    def __enter__(self):
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()
