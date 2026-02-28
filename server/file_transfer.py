"""
Remote file transfer operations.
"""

from pathlib import Path
from typing import Optional, List, Callable
import os


class FileTransfer:
    """
    File transfer operations using SFTP.
    """
    
    def __init__(self, ssh_handler):
        """
        Initialize file transfer.
        
        Args:
            ssh_handler: SSHHandler instance
        """
        self.ssh = ssh_handler
        self._sftp = None
    
    def _get_sftp(self):
        """Get SFTP client."""
        if self._sftp is None:
            self._sftp = self.ssh._client.open_sftp()
        return self._sftp
    
    def upload(
        self,
        local_path: str,
        remote_path: str,
        callback: Optional[Callable] = None
    ) -> bool:
        """
        Upload file to remote server.
        
        Args:
            local_path: Local file path
            remote_path: Remote destination path
            callback: Optional progress callback
            
        Returns:
            True if successful
        """
        try:
            sftp = self._get_sftp()
            
            # Create remote directory if needed
            remote_dir = os.path.dirname(remote_path)
            try:
                sftp.stat(remote_dir)
            except FileNotFoundError:
                self.ssh.execute(f"mkdir -p {remote_dir}")
            
            sftp.put(local_path, remote_path, callback=callback)
            return True
            
        except Exception as e:
            print(f"Upload failed: {e}")
            return False
    
    def download(
        self,
        remote_path: str,
        local_path: str,
        callback: Optional[Callable] = None
    ) -> bool:
        """
        Download file from remote server.
        
        Args:
            remote_path: Remote file path
            local_path: Local destination path
            callback: Optional progress callback
            
        Returns:
            True if successful
        """
        try:
            sftp = self._get_sftp()
            
            # Create local directory if needed
            Path(local_path).parent.mkdir(parents=True, exist_ok=True)
            
            sftp.get(remote_path, local_path, callback=callback)
            return True
            
        except Exception as e:
            print(f"Download failed: {e}")
            return False
    
    def list_remote(self, remote_dir: str) -> List[str]:
        """List files in remote directory."""
        try:
            sftp = self._get_sftp()
            return sftp.listdir(remote_dir)
        except Exception as e:
            print(f"List failed: {e}")
            return []
    
    def delete_remote(self, remote_path: str) -> bool:
        """Delete file on remote server."""
        try:
            sftp = self._get_sftp()
            sftp.remove(remote_path)
            return True
        except Exception as e:
            print(f"Delete failed: {e}")
            return False
    
    def mkdir_remote(self, remote_path: str) -> bool:
        """Create directory on remote server."""
        try:
            sftp = self._get_sftp()
            sftp.mkdir(remote_path)
            return True
        except Exception as e:
            print(f"mkdir failed: {e}")
            return False
    
    def close(self):
        """Close SFTP connection."""
        if self._sftp:
            self._sftp.close()
            self._sftp = None
