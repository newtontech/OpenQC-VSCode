"""
cclib-powered calculation output parser bridge.

Parses output files from Gaussian, ORCA, NWChem, GAMESS, Psi4, Molpro,
OpenMolcas, DALTON, Turbomole, Q-Chem, and others via cclib.

Provides CLI commands:
  parse     - Parse an output file into OpenQCResults JSON
  summarize - Return a brief summary of calculation results
  trajectory - Extract optimization trajectory frames
"""

import json
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

RY_TO_EV = 13.605693122994
_FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"

_ELEMENT_MAP = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
    9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
    16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 26: "Fe", 29: "Cu",
    30: "Zn", 47: "Ag", 79: "Au",
}


# ---------------------------------------------------------------------------
# Results DTO
# ---------------------------------------------------------------------------

def _make_results(**kwargs) -> Dict[str, Any]:
    """Create an OpenQCResults dict."""
    results = {"schemaVersion": "openqc.results.v1"}
    for k, v in kwargs.items():
        if v is not None:
            results[k] = v
    return results


def _element_symbol(atomic_number: Any, fallback: str) -> str:
    """Convert an atomic number to an element symbol when known."""
    try:
        number = int(atomic_number)
    except (TypeError, ValueError):
        return fallback
    return _ELEMENT_MAP.get(number, f"X{number}")


# ---------------------------------------------------------------------------
# cclib-based parser
# ---------------------------------------------------------------------------

def _parse_with_cclib(filepath: str, software: str = None) -> Dict[str, Any]:
    """Parse output file using cclib."""
    try:
        import cclib
    except ImportError:
        raise ImportError(
            "cclib is required for output parsing. Install: pip install cclib"
        )

    data = cclib.io.ccread(filepath)
    if data is None:
        raise ValueError(f"cclib could not parse output file: {filepath}")

    results = _make_results(
        sourceFile=filepath,
        software=software or getattr(data, "metadata", {}).get("package", "unknown"),
        success=True,
    )

    # Final energy
    if hasattr(data, "scfenergies") and data.scfenergies is not None and len(data.scfenergies) > 0:
        final_energy = float(data.scfenergies[-1])
        results["finalEnergy"] = {
            "value": final_energy,
            "unit": "eV",
        }
        results["scfEnergies"] = [float(e) for e in data.scfenergies]

    # Optimization
    if hasattr(data, "geovalues") and data.geovalues is not None and len(data.geovalues) > 0:
        opt_energies = []
        if data.scfenergies is not None and len(data.scfenergies) > 0:
            opt_energies = [float(e) for e in data.scfenergies]
        results["optimization"] = {
            "energies": opt_energies,
            "converged": len(data.geovalues) > 0,
            "stepCount": len(data.geovalues),
        }

    # Frequencies
    if hasattr(data, "vibfreqs") and data.vibfreqs is not None and len(data.vibfreqs) > 0:
        results["frequencies"] = [float(f) for f in data.vibfreqs]

    # Molecular orbitals
    if hasattr(data, "moenergies") and data.moenergies is not None and len(data.moenergies) > 0:
        mo_e = data.moenergies[0]
        results["orbitals"] = {
            "energies": [float(e) for e in mo_e],
        }
        # Try to find HOMO
        if hasattr(data, "homos") and data.homos is not None and len(data.homos) > 0:
            homo_idx = int(data.homos[0])
            results["orbitals"]["homo"] = homo_idx
            if homo_idx + 1 < len(mo_e):
                results["orbitals"]["lumo"] = homo_idx + 1

    # Charges
    if hasattr(data, "atomcharges") and data.atomcharges is not None:
        charges = {}
        for method, values in data.atomcharges.items():
            charges[method] = [float(v) for v in values]
        results["charges"] = charges

    # Dipole moment
    if hasattr(data, "moments") and data.moments is not None and len(data.moments) > 0:
        dipole = data.moments[0]
        if len(dipole) == 3:
            results["dipole"] = [float(dipole[0]), float(dipole[1]), float(dipole[2])]

    # Atom coordinates from final geometry
    if hasattr(data, "atomcoords") and data.atomcoords is not None and len(data.atomcoords) > 0:
        coords = data.atomcoords[-1]
        elements = []
        if hasattr(data, "atomnos") and data.atomnos is not None:
            # Convert atomic numbers to element symbols
            elements = [_element_symbol(n, f"X{n}") for n in data.atomnos]
        atoms = []
        for i, coord in enumerate(coords):
            elem = elements[i] if i < len(elements) else f"X{i+1}"
            atoms.append({"element": elem, "x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2])})
        results["finalStructure"] = {
            "schemaVersion": "openqc.structure.v1",
            "kind": "molecule",
            "atoms": atoms,
        }

    return results


def _extract_trajectory_with_cclib(filepath: str, software: str = None) -> Dict[str, Any]:
    """Extract geometry optimization frames as OpenQCStructure DTOs."""
    try:
        import cclib
    except ImportError as exc:
        return {
            "sourceFile": filepath,
            "software": software or "auto",
            "supported": False,
            "frames": [],
            "frameCount": 0,
            "cclibAvailable": False,
            "warnings": [f"cclib is required for trajectory extraction: {exc}"],
        }

    try:
        data = cclib.io.ccread(filepath)
    except Exception as exc:
        return {
            "sourceFile": filepath,
            "software": software or "auto",
            "supported": False,
            "frames": [],
            "frameCount": 0,
            "cclibAvailable": False,
            "warnings": [f"cclib parser failed during trajectory extraction: {exc}"],
        }

    if data is None:
        return {
            "sourceFile": filepath,
            "software": software or "auto",
            "supported": False,
            "frames": [],
            "frameCount": 0,
            "cclibAvailable": True,
            "warnings": ["cclib could not parse output file for trajectory extraction"],
        }

    software_name = software or getattr(data, "metadata", {}).get("package", "unknown")
    atomcoords = getattr(data, "atomcoords", None)
    atomnos = list(getattr(data, "atomnos", []) or [])

    if atomcoords is None or len(atomcoords) == 0:
        return {
            "sourceFile": filepath,
            "software": software_name,
            "supported": False,
            "frames": [],
            "frameCount": 0,
            "cclibAvailable": True,
            "warnings": ["No atom coordinate trajectory was found in this output file"],
        }

    frames = []
    for frame_index, coords in enumerate(atomcoords):
        atoms = []
        for atom_index, coord in enumerate(coords):
            fallback = f"X{atom_index + 1}"
            element = (
                _element_symbol(atomnos[atom_index], fallback)
                if atom_index < len(atomnos)
                else fallback
            )
            atoms.append(
                {
                    "element": element,
                    "x": float(coord[0]),
                    "y": float(coord[1]),
                    "z": float(coord[2]),
                }
            )

        frames.append(
            {
                "schemaVersion": "openqc.structure.v1",
                "name": f"Frame {frame_index + 1}",
                "kind": "molecule",
                "atoms": atoms,
            }
        )

    energies = []
    scf_energies = getattr(data, "scfenergies", None)
    if scf_energies is not None:
        energies = [float(energy) for energy in scf_energies]

    return {
        "sourceFile": filepath,
        "software": software_name,
        "supported": True,
        "frames": frames,
        "frameCount": len(frames),
        "energies": energies,
        "cclibAvailable": True,
    }


# ---------------------------------------------------------------------------
# Shared parse/summarize entry points
# ---------------------------------------------------------------------------

def parse_output_file(filepath: str, software: str = "auto") -> Dict[str, Any]:
    """Parse an output file with cclib and native fallback."""
    try:
        result = _parse_with_cclib(filepath, software)
        result["cclibAvailable"] = True
        return result
    except Exception as cclib_error:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as handle:
            content = handle.read()
        result = _parse_output_native(content, software)
        result["sourceFile"] = filepath
        result["cclibAvailable"] = False
        warnings = result.get("warnings", [])
        warnings.append(f"cclib parser unavailable or failed: {cclib_error}")
        result["warnings"] = warnings
        return result


def summarize_results(result: Dict[str, Any]) -> Dict[str, Any]:
    """Create a compact output summary from OpenQCResults-like data."""
    return {
        "sourceFile": result.get("sourceFile"),
        "software": result.get("software"),
        "success": result.get("success"),
        "finalEnergy": result.get("finalEnergy"),
        "scfStepCount": len(result.get("scfEnergies", [])),
        "optimizationConverged": result.get("optimization", {}).get("converged"),
        "frequencyCount": len(result.get("frequencies", [])),
        "cclibAvailable": result.get("cclibAvailable"),
        "warnings": result.get("warnings", []),
    }


# ---------------------------------------------------------------------------
# Native fallback for common output formats
# ---------------------------------------------------------------------------

def _parse_output_native(content: str, software: str = None) -> Dict[str, Any]:
    """Try to extract basic info from output without cclib."""
    detected_software = _detect_native_software(content, software)
    results = _make_results(software=detected_software, success=True)
    scf_energies, energy_unit = _extract_native_energies(content, detected_software)

    if scf_energies:
        results["finalEnergy"] = {"value": scf_energies[-1], "unit": energy_unit}
        results["scfEnergies"] = scf_energies
    else:
        results["success"] = False
        results["warnings"] = [
            "No calculation output data could be extracted from this file."
        ]
        results["errors"] = [
            "Unsupported or unrecognized calculation output."
        ]

    return results


def _detect_native_software(content: str, software: Optional[str]) -> str:
    """Detect common output software families for native fallback parsing."""
    hint = (software or "auto").lower()
    if hint and hint != "auto":
        return hint

    signatures = [
        ("gaussian", [r"Entering Gaussian System", r"SCF Done:"]),
        ("orca", [r"^\s*\* O\s+R\s+C\s+A \*", r"FINAL SINGLE POINT ENERGY", r"ORCA TERMINATED"]),
        ("cp2k", [r"^CP2K\|", r"^ENERGY\|"]),
        ("qe", [r"Program PWSCF", r"!\s+total energy\s+="]),
        ("gamess", [r"GAMESS VERSION", r"^\s*TOTAL ENERGY\s*="]),
    ]

    for software_name, patterns in signatures:
        if any(re.search(pattern, content, re.IGNORECASE | re.MULTILINE) for pattern in patterns):
            return software_name

    return software or "auto"


def _extract_native_energies(content: str, software: str) -> Tuple[List[float], str]:
    """Extract final/scf energies from common text outputs."""
    software_lower = (software or "auto").lower()

    if software_lower == "gaussian":
        return _extract_gaussian_energies(content), "hartree"
    if software_lower == "orca":
        return _extract_regex_energies(
            content,
            [
                rf"FINAL SINGLE POINT ENERGY\s+({_FLOAT_RE})",
                rf"Total Energy\s*:\s*({_FLOAT_RE})\s*Eh",
            ],
        ), "hartree"
    if software_lower == "cp2k":
        return _extract_regex_energies(
            content,
            [
                rf"^ENERGY\|.*?({_FLOAT_RE})\s*$",
                rf"Total energy:\s*({_FLOAT_RE})",
            ],
        ), "hartree"
    if software_lower in {"qe", "quantum_espresso"}:
        energies_ry = _extract_regex_energies(
            content,
            [
                rf"!\s+total energy\s+=\s+({_FLOAT_RE})\s+Ry",
                rf"total energy\s+=\s+({_FLOAT_RE})\s+Ry",
            ],
        )
        return [energy * RY_TO_EV for energy in energies_ry], "eV"
    if software_lower == "gamess":
        return _extract_regex_energies(
            content,
            [
                rf"^\s*TOTAL ENERGY\s*=\s*({_FLOAT_RE})",
                rf"FINAL\s+\S+\s+ENERGY\s+IS\s+({_FLOAT_RE})",
            ],
        ), "hartree"

    gaussian_energies = _extract_gaussian_energies(content)
    if gaussian_energies:
        return gaussian_energies, "hartree"

    return _extract_regex_energies(content, [rf"converged to\s+({_FLOAT_RE})"]), "hartree"


def _extract_gaussian_energies(content: str) -> List[float]:
    energies: List[float] = []
    for match in re.finditer(rf"SCF Done:\s+.*?=\s+({_FLOAT_RE})", content, re.IGNORECASE):
        energies.append(_parse_float(match.group(1)))
    if energies:
        return energies

    return _extract_regex_energies(content, [rf"converged to\s+({_FLOAT_RE})"])


def _extract_regex_energies(content: str, patterns: List[str]) -> List[float]:
    energies: List[float] = []
    for pattern in patterns:
        for match in re.finditer(pattern, content, re.IGNORECASE | re.MULTILINE):
            try:
                energies.append(_parse_float(match.group(1)))
            except (IndexError, ValueError):
                continue
        if energies:
            return energies
    return energies


def _parse_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    """CLI entry point for output bridge."""
    from openqc.bridge.protocol import read_request, write_response, write_error, get_command, get_args

    request = read_request()
    command = get_command(request)
    args = get_args(request)

    try:
        if command == "parse":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            software = args.get("software", "auto")
            result = parse_output_file(filepath, software)
            write_response(result)

        elif command == "summarize":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            result = parse_output_file(filepath, args.get("software", "auto"))
            summary = summarize_results(result)
            write_response(summary)

        elif command == "trajectory":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            result = _extract_trajectory_with_cclib(filepath, args.get("software"))
            write_response(result)

        else:
            write_error(f"Unknown command: {command}")

    except Exception as e:
        write_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
