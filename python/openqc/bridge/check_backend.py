"""
Check scientific Python backend status.

Reports Python version, available scientific packages,
and external tool availability.
"""

import importlib
import json
import shutil
import subprocess
import sys
import platform
from typing import Any, Dict, List, Optional


def _get_python_info() -> Dict[str, str]:
    """Get Python interpreter information."""
    return {
        "executable": sys.executable,
        "version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "platform": platform.platform(),
    }


def _check_package(name: str, install_hint: str = "") -> Dict[str, Any]:
    """Check if a Python package is available and return its version."""
    result: Dict[str, Any] = {"available": False}
    if install_hint:
        result["installHint"] = install_hint

    try:
        mod = importlib.import_module(name)
        version = getattr(mod, "__version__", "unknown")
        result["available"] = True
        result["version"] = version
    except ImportError:
        pass

    return result


def _check_external_tool(name: str) -> Dict[str, Any]:
    """Check if an external tool is available on PATH."""
    path = shutil.which(name)
    return {
        "available": path is not None,
        "path": path,
    }


def check_backend() -> Dict[str, Any]:
    """Run full backend check and return structured result."""
    # Core Python
    python_info = _get_python_info()

    # Scientific packages
    packages = {
        "pymatgen": _check_package("pymatgen", "pip install pymatgen"),
        "ase": _check_package("ase", "pip install ase"),
        "cclib": _check_package("cclib", "pip install cclib"),
        "dpdata": _check_package("dpdata", "pip install dpdata"),
        "spglib": _check_package("spglib", "pip install spglib"),
        "rdkit": _check_package("rdkit", "conda install -c conda-forge rdkit"),
        "numpy": _check_package("numpy", "pip install numpy"),
    }

    # External tools
    external_tools = {
        "multiwfn": _check_external_tool("Multiwfn"),
        "openbabel": _check_external_tool("obabel"),
        "c2x": _check_external_tool("c2x"),
    }

    return {
        "success": True,
        "python": python_info,
        "packages": packages,
        "externalTools": external_tools,
    }


def main():
    """Entry point for check_backend command."""
    result = check_backend()
    sys.stdout.write(json.dumps(result) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
