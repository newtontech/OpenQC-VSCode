"""
Slurm job management interface.
"""

from typing import Optional, List, Dict, Any
from dataclasses import dataclass
from datetime import datetime


@dataclass
class SlurmJob:
    """Represents a Slurm job."""
    job_id: str
    name: str
    user: str
    state: str
    partition: str
    num_nodes: int
    num_cpus: int
    time_limit: str
    time_used: str
    submit_time: Optional[datetime] = None
    start_time: Optional[datetime] = None


class SlurmInterface:
    """
    Interface for Slurm workload manager.
    """
    
    def __init__(self, ssh_handler):
        """
        Initialize Slurm interface.
        
        Args:
            ssh_handler: SSHHandler instance for remote execution
        """
        self.ssh = ssh_handler
    
    def submit_job(
        self,
        script_path: str,
        job_name: Optional[str] = None,
        **kwargs
    ) -> str:
        """
        Submit a job script to Slurm.
        
        Args:
            script_path: Path to job script
            job_name: Optional job name override
            
        Returns:
            Job ID
        """
        cmd = f"sbatch {script_path}"
        if job_name:
            cmd += f" --job-name={job_name}"
        
        exit_code, stdout, stderr = self.ssh.execute(cmd)
        
        if exit_code != 0:
            raise RuntimeError(f"Job submission failed: {stderr}")
        
        # Parse job ID from output
        # Output format: "Submitted batch job 12345"
        job_id = stdout.strip().split()[-1]
        return job_id
    
    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running job."""
        exit_code, stdout, stderr = self.ssh.execute(f"scancel {job_id}")
        return exit_code == 0
    
    def get_job_status(self, job_id: str) -> Optional[SlurmJob]:
        """Get status of a specific job."""
        cmd = f"squeue -j {job_id} --format='%i %j %u %T %P %D %C %l %M' --noheader"
        exit_code, stdout, stderr = self.ssh.execute(cmd)
        
        if exit_code != 0 or not stdout.strip():
            return None
        
        parts = stdout.strip().split()
        if len(parts) >= 9:
            return SlurmJob(
                job_id=parts[0],
                name=parts[1],
                user=parts[2],
                state=parts[3],
                partition=parts[4],
                num_nodes=int(parts[5]),
                num_cpus=int(parts[6]),
                time_limit=parts[7],
                time_used=parts[8]
            )
        return None
    
    def list_jobs(self, user: Optional[str] = None) -> List[SlurmJob]:
        """List all jobs (optionally filtered by user)."""
        cmd = "squeue --format='%i %j %u %T %P %D %C %l %M' --noheader"
        if user:
            cmd += f" --user={user}"
        
        exit_code, stdout, stderr = self.ssh.execute(cmd)
        
        if exit_code != 0:
            return []
        
        jobs = []
        for line in stdout.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 9:
                jobs.append(SlurmJob(
                    job_id=parts[0],
                    name=parts[1],
                    user=parts[2],
                    state=parts[3],
                    partition=parts[4],
                    num_nodes=int(parts[5]),
                    num_cpus=int(parts[6]),
                    time_limit=parts[7],
                    time_used=parts[8]
                ))
        
        return jobs
    
    def get_queue_info(self) -> Dict[str, Any]:
        """Get queue/partition information."""
        cmd = "sinfo --format='%P %a %l %D %T' --noheader"
        exit_code, stdout, stderr = self.ssh.execute(cmd)
        
        if exit_code != 0:
            return {}
        
        partitions = {}
        for line in stdout.strip().split('\n'):
            parts = line.split()
            if len(parts) >= 5:
                partitions[parts[0]] = {
                    'availability': parts[1],
                    'timelimit': parts[2],
                    'nodes': parts[3],
                    'state': parts[4]
                }
        
        return partitions
