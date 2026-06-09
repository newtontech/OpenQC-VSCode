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
import sys
from typing import Any, Dict, List, Optional


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
            _ELEMENT_MAP = {
                1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O",
                9: "F", 10: "Ne", 11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P",
                16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 26: "Fe", 29: "Cu",
                30: "Zn", 47: "Ag", 79: "Au",
            }
            elements = [_ELEMENT_MAP.get(n, f"X{n}") for n in data.atomnos]
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


# ---------------------------------------------------------------------------
# Native fallback for common output formats
# ---------------------------------------------------------------------------

def _parse_output_native(content: str, software: str = None) -> Dict[str, Any]:
    """Try to extract basic info from output without cclib."""
    results = _make_results(sourceSoftware=software)

    # Try to find final energy
    for line in reversed(content.split("\n")):
        line_lower = line.lower().strip()
        if "scf done" in line_lower or "converged to" in line_lower:
            parts = line.split()
            for i, p in enumerate(parts):
                try:
                    val = float(p)
                    if abs(val) > 0.1 and abs(val) < 1e6:
                        results["finalEnergy"] = {"value": val, "unit": "hartree"}
                        break
                except ValueError:
                    continue
            break

    return results


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
            try:
                result = _parse_with_cclib(filepath, software)
            except ImportError:
                # Fallback: read file and try native extraction
                with open(filepath, "r") as f:
                    content = f.read()
                result = _parse_output_native(content, software)
                result["cclibAvailable"] = False
            write_response(result)

        elif command == "summarize":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            result = _parse_with_cclib(filepath, args.get("software"))
            summary = {
                "sourceFile": result.get("sourceFile"),
                "software": result.get("software"),
                "success": result.get("success"),
                "finalEnergy": result.get("finalEnergy"),
                "scfStepCount": len(result.get("scfEnergies", [])),
                "optimizationConverged": result.get("optimization", {}).get("converged"),
                "frequencyCount": len(result.get("frequencies", [])),
            }
            write_response(summary)

        else:
            write_error(f"Unknown command: {command}")

    except Exception as e:
        write_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
