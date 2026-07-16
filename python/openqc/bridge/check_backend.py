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

    capabilities = _derive_capabilities(packages, external_tools)
    required_packages = ["numpy", "ase", "pymatgen", "cclib", "dpdata"]
    missing_packages = [
        name for name in required_packages
        if not packages.get(name, {}).get("available", False)
    ]

    if len(missing_packages) == len(required_packages):
        status = "missing"
    elif missing_packages:
        status = "degraded"
    else:
        status = "installed"

    status_detail = _status_detail(status, missing_packages)

    return {
        "success": True,
        "status": status,
        "statusDetail": status_detail,
        "python": python_info,
        "packages": packages,
        "externalTools": external_tools,
        "capabilities": capabilities,
        "missingPackages": missing_packages,
    }


def _status_detail(status: str, missing_packages: List[str]) -> str:
    """Return a concise human-facing backend health detail."""
    if status == "installed":
        return "Core scientific Python packages are installed."

    if status == "degraded":
        missing = ", ".join(missing_packages)
        return (
            "Python is available, but these core packages are missing: "
            f"{missing}. Related OpenQC features are degraded or disabled."
        )

    return (
        "Python is available, but no core scientific backend packages were found. "
        "Install numpy plus at least one OpenQC backend package such as ase, pymatgen, cclib, or dpdata."
    )


def _capability(
    label: str,
    status: str,
    detail: str,
    requires: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Build a structured capability status entry."""
    result: Dict[str, Any] = {"label": label, "status": status, "detail": detail}
    if requires:
        result["requires"] = requires
    return result


def _status_for(available: List[bool]) -> str:
    """Map dependency availability to available/degraded/missing."""
    if all(available):
        return "available"
    if any(available):
        return "degraded"
    return "missing"


def _missing_requirements(requirements: Dict[str, bool]) -> List[str]:
    """Return dependency names that are not available."""
    return [name for name, available in requirements.items() if not available]


def _derive_capabilities(
    packages: Dict[str, Dict[str, Any]],
    external_tools: Dict[str, Dict[str, Any]],
) -> Dict[str, Dict[str, Any]]:
    """Summarize feature-level backend readiness."""
    has_ase = packages["ase"]["available"]
    has_pymatgen = packages["pymatgen"]["available"]
    has_cclib = packages["cclib"]["available"]
    has_dpdata = packages["dpdata"]["available"]
    has_spglib = packages["spglib"]["available"]
    has_multiwfn = external_tools["multiwfn"]["available"]
    has_c2x = external_tools["c2x"]["available"]
    has_openbabel = external_tools["openbabel"]["available"]

    structure_requirements = {"native": True, "ase": has_ase, "pymatgen": has_pymatgen}
    conversion_requirements = {"ase": has_ase, "dpdata": has_dpdata}
    analysis_requirements = {
        "Multiwfn": has_multiwfn,
        "c2x": has_c2x,
        "Open Babel": has_openbabel,
    }

    return {
        "structureParsing": _capability(
            "Structure parsing",
            _status_for(list(structure_requirements.values())),
            "Native parsers cover XYZ, POSCAR/CONTCAR, Gaussian, ORCA, CP2K, and GAMESS coordinate inputs; ASE and pymatgen broaden periodic and molecular structure support.",
            _missing_requirements({"ase": has_ase, "pymatgen": has_pymatgen}),
        ),
        "formatConversion": _capability(
            "Format conversion",
            _status_for(list(conversion_requirements.values())),
            "ASE handles broad structure-file conversion; dpdata handles molecular dynamics dataset conversion.",
            _missing_requirements(conversion_requirements),
        ),
        "outputParsing": _capability(
            "Calculation output parsing",
            "available" if has_cclib else "missing",
            "cclib parses calculation outputs for energies, orbitals, frequencies, charges, and structures.",
            [] if has_cclib else ["cclib"],
        ),
        "supercell": _capability(
            "Supercell generation",
            "available" if (has_ase or has_pymatgen) else "missing",
            "ASE or pymatgen enables backend supercell generation for parsed periodic structures.",
            [] if (has_ase or has_pymatgen) else ["ase", "pymatgen"],
        ),
        "standardization": _capability(
            "Symmetry standardization",
            "available" if has_spglib else "missing",
            "spglib enables symmetry-aware standardization workflows.",
            [] if has_spglib else ["spglib"],
        ),
        "externalAnalysis": _capability(
            "External analysis tools",
            _status_for(list(analysis_requirements.values())),
            "Multiwfn, c2x, and Open Babel are optional command-line analyzers for wavefunction, density, and structure conversion workflows.",
            _missing_requirements(analysis_requirements),
        ),
    }


def main():
    """Entry point for check_backend command."""
    result = check_backend()
    sys.stdout.write(json.dumps(result) + "\n")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
