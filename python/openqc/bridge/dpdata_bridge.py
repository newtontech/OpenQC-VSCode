"""
dpdata bridge for trajectories, MD data, and ML potential datasets.

Provides CLI commands:
  summarize  - Summarize a dataset (frame count, fields, energy range)
  trajectory - Extract frames as OpenQCStructure list
  check      - Quality check on dataset
"""

import json
import sys
from typing import Any, Dict, List


def _summarize_dataset(path: str, fmt: str = "auto") -> Dict[str, Any]:
    """Summarize an atomistic dataset."""
    try:
        import dpdata
    except ImportError:
        raise ImportError("dpdata is required. Install: pip install dpdata")

    if fmt == "auto":
        ds = dpdata.System(path, fmt="auto")
    else:
        ds = dpdata.System(path, fmt=fmt)

    summary = {
        "schemaVersion": "openqc.dataset.v1",
        "sourcePath": path,
        "sourceFormat": fmt,
        "frameCount": len(ds),
        "atomCount": ds.get_natoms() if hasattr(ds, "get_natoms") else None,
        "elements": list(ds.get_atom_names()) if hasattr(ds, "get_atom_names") else [],
        "hasEnergy": "energies" in ds.data,
        "hasForce": "forces" in ds.data,
        "hasVirial": "virials" in ds.data,
        "hasCell": "cells" in ds.data,
    }

    if "energies" in ds.data:
        energies = ds.data["energies"]
        summary["energyRange"] = [float(energies.min()), float(energies.max())]

    if "forces" in ds.data:
        import numpy as np
        forces = ds.data["forces"]
        norms = np.linalg.norm(forces.reshape(-1, 3), axis=1)
        summary["maxForceNorm"] = float(norms.max())

    return summary


def _extract_trajectory(path: str, fmt: str = "auto", max_frames: int = 100) -> Dict[str, Any]:
    """Extract trajectory frames as OpenQCStructure list."""
    try:
        import dpdata
    except ImportError:
        raise ImportError("dpdata is required. Install: pip install dpdata")

    if fmt == "auto":
        ds = dpdata.System(path, fmt="auto")
    else:
        ds = dpdata.System(path, fmt=fmt)

    frames = []
    for idx in range(min(len(ds), max_frames)):
        frame = ds[idx]
        atoms = []
        atom_names = frame.get_atom_names() if hasattr(frame, "get_atom_names") else []
        coords = frame.get_positions() if hasattr(frame, "get_positions") else []

        for i, coord in enumerate(coords):
            elem = atom_names[i] if i < len(atom_names) else f"X{i+1}"
            atoms.append({"element": elem, "x": float(coord[0]), "y": float(coord[1]), "z": float(coord[2])})

        frames.append({
            "schemaVersion": "openqc.structure.v1",
            "kind": "molecule",
            "atoms": atoms,
            "name": f"Frame {idx}",
        })

    return {
        "frameCount": len(frames),
        "totalFrames": len(ds),
        "frames": frames,
    }


def _check_quality(path: str, fmt: str = "auto") -> Dict[str, Any]:
    """Run quality checks on a dataset."""
    try:
        import dpdata
        import numpy as np
    except ImportError:
        raise ImportError("dpdata and numpy are required. Install: pip install dpdata numpy")

    if fmt == "auto":
        ds = dpdata.System(path, fmt="auto")
    else:
        ds = dpdata.System(path, fmt=fmt)

    warnings = []
    errors = []

    if len(ds) == 0:
        errors.append("Dataset is empty")
        return {"schemaVersion": "openqc.dataset.v1", "errors": errors, "warnings": warnings}

    natoms = ds.get_natoms() if hasattr(ds, "get_natoms") else None
    if natoms and natoms == 0:
        errors.append("Dataset has zero atoms")

    if "energies" in ds.data:
        energies = ds.data["energies"]
        if np.any(np.isnan(energies)):
            warnings.append("NaN energies detected")
        if np.any(np.isinf(energies)):
            warnings.append("Inf energies detected")

    if "forces" in ds.data:
        forces = ds.data["forces"]
        if np.any(np.isnan(forces)):
            warnings.append("NaN forces detected")
        norms = np.linalg.norm(forces.reshape(-1, 3), axis=1)
        if norms.max() > 100:
            warnings.append(f"Suspiciously large force: {norms.max():.2f}")

    return {
        "schemaVersion": "openqc.dataset.v1",
        "frameCount": len(ds),
        "warnings": warnings,
        "errors": errors,
        "passed": len(errors) == 0,
    }


def main():
    from openqc.bridge.protocol import read_request, write_response, write_error, get_command, get_args

    request = read_request()
    command = get_command(request)
    args = get_args(request)

    try:
        if command == "summarize":
            path = args.get("path")
            if not path:
                raise ValueError("Missing 'path'")
            result = _summarize_dataset(path, args.get("format", "auto"))
            write_response(result)

        elif command == "trajectory":
            path = args.get("path")
            if not path:
                raise ValueError("Missing 'path'")
            result = _extract_trajectory(path, args.get("format", "auto"), args.get("maxFrames", 100))
            write_response(result)

        elif command == "check":
            path = args.get("path")
            if not path:
                raise ValueError("Missing 'path'")
            result = _check_quality(path, args.get("format", "auto"))
            write_response(result)

        else:
            write_error(f"Unknown command: {command}")

    except Exception as e:
        write_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
