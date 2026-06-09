"""
Structure bridge using pymatgen/ASE for parsing, conversion, and supercell workflows.

Provides CLI commands:
  parse    - Parse a structure file into OpenQCStructure JSON
  convert  - Convert a structure to another format
  supercell - Generate a supercell preview
  standardize - Standardize a structure using spglib if available
"""

import json
import sys
import os
from typing import Any, Dict, List, Optional, Tuple


def _make_atom(element: str, x: float, y: float, z: float, **extra) -> Dict[str, Any]:
    """Create an atom dict for OpenQCStructure."""
    atom: Dict[str, Any] = {"element": element, "x": x, "y": y, "z": z}
    for k in ("fx", "fy", "fz", "charge", "magmom", "label"):
        if k in extra:
            atom[k] = extra[k]
    return atom


def _make_cell(a, b, c, pbc=None, coordinate_mode=None) -> Dict[str, Any]:
    """Create a cell dict for OpenQCStructure."""
    cell: Dict[str, Any] = {
        "a": list(a),
        "b": list(b),
        "c": list(c),
        "pbc": pbc or [True, True, True],
    }
    if coordinate_mode:
        cell["coordinateMode"] = coordinate_mode
    return cell


def _make_structure(atoms, cell=None, name=None, kind=None, source_format=None,
                    source_software=None, warnings=None) -> Dict[str, Any]:
    """Create an OpenQCStructure dict."""
    if kind is None:
        kind = "periodic" if cell else "molecule"
    structure = {
        "schemaVersion": "openqc.structure.v1",
        "kind": kind,
        "atoms": atoms,
    }
    if name:
        structure["name"] = name
    if cell:
        structure["cell"] = cell
    structure["metadata"] = {
        "source": {
            "format": source_format,
            "software": source_software,
            "parser": "python-structure-bridge",
        },
        "provenance": {
            "createdAt": __import__("datetime").datetime.utcnow().isoformat() + "Z",
            "warnings": warnings or [],
        },
    }
    return structure


# ---------------------------------------------------------------------------
# Native parsers (no dependency required)
# ---------------------------------------------------------------------------

def _parse_xyz(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse XYZ format natively."""
    lines = content.strip().split("\n")
    if len(lines) < 3:
        raise ValueError("Invalid XYZ: too few lines")
    n_atoms = int(lines[0].strip())
    comment = lines[1].strip()
    atoms = []
    for i in range(2, min(len(lines), n_atoms + 2)):
        parts = lines[i].strip().split()
        if len(parts) >= 4:
            atoms.append(_make_atom(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return _make_structure(atoms, name=filename or comment, source_format="xyz")


def _parse_poscar(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse VASP POSCAR/CONTCAR natively."""
    lines = content.strip().split("\n")
    if len(lines) < 8:
        raise ValueError("Invalid POSCAR: too few lines")

    comment = lines[0].strip()
    scale = float(lines[1].strip())
    if scale == 0:
        scale = 1.0

    def parse_vec(line):
        parts = line.strip().split()
        return [float(parts[0]) * scale, float(parts[1]) * scale, float(parts[2]) * scale]

    a = parse_vec(lines[2])
    b = parse_vec(lines[3])
    c = parse_vec(lines[4])

    type_parts = lines[5].strip().split()
    count_parts = lines[6].strip().split()

    has_elem = all(p[0].isalpha() for p in type_parts)
    if has_elem:
        elements = type_parts
        counts = [int(x) for x in count_parts]
    else:
        counts = [int(x) for x in type_parts]
        elements = [f"El{i+1}" for i in range(len(counts))]

    mode_line = lines[7].strip().lower()
    is_selective = mode_line.startswith("s")
    is_direct = mode_line.startswith("d") or (not is_selective and not mode_line.startswith("c") and mode_line)
    coord_line = 8 if (is_selective or is_direct or mode_line.startswith("c")) else 7
    coord_mode = "fractional" if is_direct else "cartesian"

    atoms = []
    idx = coord_line
    for i, (elem, count) in enumerate(zip(elements, counts)):
        for j in range(count):
            if idx >= len(lines):
                break
            parts = lines[idx].strip().split()
            idx += 1
            if len(parts) < 3:
                continue
            atom = _make_atom(elem, float(parts[0]), float(parts[1]), float(parts[2]))
            if is_selective and len(parts) >= 6:
                atom["selectiveDynamics"] = [
                    parts[3].upper().startswith("T"),
                    parts[4].upper().startswith("T"),
                    parts[5].upper().startswith("T"),
                ]
            atoms.append(atom)

    cell = _make_cell(a, b, c, coordinate_mode=coord_mode)
    return _make_structure(atoms, cell=cell, name=filename or comment,
                           source_format="poscar", source_software="vasp")


def _parse_gaussian_input(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse Gaussian input file (.gjf/.com) natively."""
    lines = content.split("\n")
    atoms = []
    in_atoms = False
    for line in lines:
        stripped = line.strip()
        if not in_atoms:
            # Look for charge/multiplicity line
            if stripped and all(p.lstrip("-").isdigit() for p in stripped.split()) and len(stripped.split()) == 2:
                in_atoms = True
                continue
        else:
            if stripped == "":
                break
            parts = stripped.split()
            if len(parts) >= 4:
                elem = parts[0]
                if elem[0].isalpha():
                    atoms.append(_make_atom(elem, float(parts[1]), float(parts[2]), float(parts[3])))
    return _make_structure(atoms, name=filename, source_format="gaussian", source_software="gaussian")


def _parse_orca_input(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse ORCA input file natively."""
    lines = content.split("\n")
    atoms = []
    in_atoms = False
    for line in lines:
        if "* xyz" in line.lower() or "* xyzfile" in line.lower():
            in_atoms = True
            continue
        if in_atoms and line.strip() == "*":
            break
        if in_atoms:
            parts = line.strip().split()
            if len(parts) >= 4 and parts[0][0].isalpha():
                atoms.append(_make_atom(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return _make_structure(atoms, name=filename, source_format="orca", source_software="orca")


# ---------------------------------------------------------------------------
# ASE-based parsers (optional dependency)
# ---------------------------------------------------------------------------

def _parse_with_ase(filepath: str, format_hint: str = None) -> Dict[str, Any]:
    """Parse structure using ASE (optional)."""
    try:
        from ase.io import read as ase_read
    except ImportError:
        raise ImportError("ASE is required for this format. Install: pip install ase")

    atoms_obj = ase_read(filepath, format=format_hint)
    atoms = []
    for i, atom in enumerate(atoms_obj):
        pos = atom.position
        atoms.append(_make_atom(str(atom.symbol), float(pos[0]), float(pos[1]), float(pos[2])))

    cell = None
    if any(atoms_obj.pbc):
        cell_vecs = atoms_obj.get_cell()
        cell = _make_cell(
            cell_vecs[0].tolist(), cell_vecs[1].tolist(), cell_vecs[2].tolist(),
            pbc=atoms_obj.pbc.tolist(),
        )

    return _make_structure(atoms, cell=cell, name=os.path.basename(filepath),
                           source_format=format_hint or "ase")


# ---------------------------------------------------------------------------
# pymatgen-based parsers (optional dependency)
# ---------------------------------------------------------------------------

def _parse_with_pymatgen(filepath: str, format_hint: str = None) -> Dict[str, Any]:
    """Parse structure using pymatgen (optional)."""
    try:
        from pymatgen.core import Structure as PymatgenStructure
        from pymatgen.io.cif import CifParser
    except ImportError:
        raise ImportError("pymatgen is required for this format. Install: pip install pymatgen")

    ext = format_hint or os.path.splitext(filepath)[1].lower()

    if ext in (".cif", "cif"):
        parser = CifParser(filepath)
        struct = parser.parse_structures()[0]
    else:
        struct = PymatgenStructure.from_file(filepath)

    atoms = []
    for i, site in enumerate(struct):
        atoms.append(_make_atom(str(site.specie), float(site.x), float(site.y), float(site.z)))

    lattice = struct.lattice
    cell = _make_cell(
        lattice.matrix[0].tolist(),
        lattice.matrix[1].tolist(),
        lattice.matrix[2].tolist(),
    )

    return _make_structure(atoms, cell=cell, name=os.path.basename(filepath),
                           source_format=ext, source_software="pymatgen")


# ---------------------------------------------------------------------------
# Supercell generation
# ---------------------------------------------------------------------------

def _generate_supercell(structure: Dict[str, Any], na: int, nb: int, nc: int) -> Dict[str, Any]:
    """Generate a supercell from a periodic structure."""
    if not structure.get("cell"):
        raise ValueError("Structure must have a cell for supercell generation")

    cell = structure["cell"]
    base_atoms = structure["atoms"]
    max_atoms = 10000
    total = len(base_atoms) * na * nb * nc
    if total > max_atoms:
        raise ValueError(f"Supercell would have {total} atoms (max {max_atoms})")

    super_atoms = []
    a = cell["a"]
    b = cell["b"]
    c = cell["c"]

    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                for atom in base_atoms:
                    dx = ia * a[0] + ib * b[0] + ic * c[0]
                    dy = ia * a[1] + ib * b[1] + ic * c[1]
                    dz = ia * a[2] + ib * b[2] + ic * c[2]
                    super_atoms.append(_make_atom(
                        atom["element"],
                        atom["x"] + dx,
                        atom["y"] + dy,
                        atom["z"] + dz,
                    ))

    super_cell = _make_cell(
        [a[0] * na, a[1] * na, a[2] * na],
        [b[0] * nb, b[1] * nb, b[2] * nb],
        [c[0] * nc, c[1] * nc, c[2] * nc],
        pbc=cell.get("pbc", [True, True, True]),
    )

    return _make_structure(
        super_atoms, cell=super_cell,
        name=f"{structure.get('name', 'supercell')} {na}x{nb}x{nc}",
        source_format=structure.get("metadata", {}).get("source", {}).get("format"),
    )


# ---------------------------------------------------------------------------
# Auto-detect format and parse
# ---------------------------------------------------------------------------

def _detect_format(filepath: str) -> str:
    """Detect file format from extension and name."""
    basename = os.path.basename(filepath).upper()
    ext = os.path.splitext(filepath)[1].lower()

    if basename in ("POSCAR", "CONTCAR"):
        return "poscar"
    if ext == ".cif":
        return "cif"
    if ext == ".xyz":
        return "xyz"
    if ext in (".gjf", ".com"):
        return "gaussian"
    if ext == ".inp":
        return "orca"  # default for .inp
    if basename.startswith("INCAR") or basename.startswith("KPOINTS"):
        return "vasp"
    return ext.strip(".") or "unknown"


def parse_file(filepath: str, format_hint: str = None, content: str = None) -> Dict[str, Any]:
    """Parse a structure file into OpenQCStructure JSON."""
    fmt = format_hint or _detect_format(filepath)
    name = os.path.basename(filepath)

    # Native parsers first (no dependency needed)
    if fmt == "xyz" and content:
        return _parse_xyz(content, name)
    if fmt == "poscar" and content:
        return _parse_poscar(content, name)
    if fmt == "gaussian" and content:
        return _parse_gaussian_input(content, name)
    if fmt == "orca" and content:
        return _parse_orca_input(content, name)

    # Try pymatgen for CIF and other crystal formats
    if fmt in ("cif",) or (not content):
        try:
            return _parse_with_pymatgen(filepath, fmt)
        except ImportError:
            pass

    # Try ASE as general fallback
    try:
        return _parse_with_ase(filepath, fmt)
    except ImportError:
        pass

    raise ValueError(
        f"Cannot parse format '{fmt}'. Install pymatgen or ASE for extended format support: "
        "pip install pymatgen ase"
    )


def main():
    """CLI entry point for structure bridge."""
    from openqc.bridge.protocol import read_request, write_response, write_error, get_command, get_args

    request = read_request()
    command = get_command(request)
    args = get_args(request)

    try:
        if command == "parse":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            fmt = args.get("format", "auto")
            content = args.get("content")
            if fmt == "auto":
                fmt = None
            result = parse_file(filepath, format_hint=fmt, content=content)
            write_response(result)

        elif command == "supercell":
            structure = args.get("structure")
            if not structure:
                raise ValueError("Missing 'structure' argument")
            na = args.get("na", 2)
            nb = args.get("nb", 2)
            nc = args.get("nc", 2)
            result = _generate_supercell(structure, na, nb, nc)
            write_response(result)

        elif command == "convert":
            filepath = args.get("path")
            if not filepath:
                raise ValueError("Missing 'path' argument")
            target_format = args.get("to", "xyz")
            # Parse first, then convert
            fmt = args.get("format", "auto")
            if fmt == "auto":
                fmt = None
            structure = parse_file(filepath, format_hint=fmt)
            # For now, output as the structure JSON (conversion to specific formats
            # can be added via ASE/pymatgen later)
            write_response(structure)

        else:
            write_error(f"Unknown command: {command}")

    except Exception as e:
        write_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
