"""
SSH connection management for remote servers.
"""

import paramiko
from typing import Optional, Dict, Any
from pathlib import Path


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
        password: Optional[str] = None
    ):
        """
        Initialize SSH handler.
        
        Args:
            host: Server hostname or IP
            username: SSH username
            port: SSH port (default 22)
            key_path: Path to SSH private key
            password: Password (if not using key)
        """
        self.host = host
        self.username = username
        self.port = port
        self.key_path = key_path
        self.password = password
        
        self._client: Optional[paramiko.SSHClient] = None
    
    def connect(self) -> bool:
        """Establish SSH connection."""
        try:
            self._client = paramiko.SSHClient()
            self._client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            
            connect_kwargs = {
                'hostname': self.host,
                'port': self.port,
                'username': self.username,
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
    
    def execute(self, command: str) -> tuple[int, str, str]:
        """
        Execute command on remote server.
        
        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        if not self._client:
            raise RuntimeError("Not connected to server")
        
        stdin, stdout, stderr = self._client.exec_command(command)
        exit_code = stdout.channel.recv_exit_status()
        
        return (
            exit_code,
            stdout.read().decode('utf-8'),
            stderr.read().decode('utf-8')
        )
    
    @property
    def is_connected(self) -> bool:
        """Check if connection is active."""
        return self._client is not None and self._client.get_transport() is not None
    
    def __enter__(self):
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.disconnect()
