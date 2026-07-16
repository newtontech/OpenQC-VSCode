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
import math
import re
from datetime import datetime, timezone
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
            "createdAt": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
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


def _parse_cp2k_input(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse CP2K input coordinates natively from &SUBSYS/&COORD blocks."""
    atoms = []
    cell_vectors: Dict[str, List[float]] = {}
    coordinate_mode = "cartesian"
    in_cell = False
    in_coord = False

    for raw_line in content.splitlines():
        line = raw_line.split("#", 1)[0].split("!", 1)[0].strip()
        upper = line.upper()
        if not line:
            continue

        if upper.startswith("&CELL"):
            in_cell = True
            continue
        if in_cell and upper.startswith("&END"):
            in_cell = False
            continue
        if upper.startswith("&COORD"):
            in_coord = True
            continue
        if in_coord and upper.startswith("&END"):
            in_coord = False
            continue

        if in_cell:
            parts = line.split()
            key = parts[0].upper()
            if key == "ABC" and len(parts) >= 4:
                a, b, c = (float(parts[1]), float(parts[2]), float(parts[3]))
                cell_vectors["a"] = [a, 0.0, 0.0]
                cell_vectors["b"] = [0.0, b, 0.0]
                cell_vectors["c"] = [0.0, 0.0, c]
            elif key in ("A", "B", "C") and len(parts) >= 4:
                cell_vectors[key.lower()] = [float(parts[1]), float(parts[2]), float(parts[3])]
            continue

        if in_coord:
            if upper in ("SCALED", "FRACTIONAL"):
                coordinate_mode = "fractional"
                continue
            if upper.startswith("UNIT"):
                continue
            parts = line.split()
            if len(parts) >= 4 and parts[0][0].isalpha():
                try:
                    atoms.append(_make_atom(parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                except ValueError:
                    continue

    cell = None
    if all(key in cell_vectors for key in ("a", "b", "c")):
        cell = _make_cell(
            cell_vectors["a"],
            cell_vectors["b"],
            cell_vectors["c"],
            coordinate_mode=coordinate_mode,
        )

    return _make_structure(atoms, cell=cell, name=filename, source_format="cp2k", source_software="cp2k")


def _parse_gamess_input(content: str, filename: str = "") -> Dict[str, Any]:
    """Parse GAMESS input coordinates natively from the $DATA group."""
    atoms = []
    in_data = False
    skipped_header_lines = 0

    for raw_line in content.splitlines():
        stripped = raw_line.strip()
        upper = stripped.upper()
        if not stripped:
            continue
        if "$DATA" in upper:
            in_data = True
            skipped_header_lines = 0
            continue
        if in_data and "$END" in upper:
            break
        if not in_data:
            continue
        if skipped_header_lines < 2:
            skipped_header_lines += 1
            continue

        parts = stripped.split()
        if len(parts) >= 5 and parts[0][0].isalpha() and _is_float(parts[1]):
            try:
                atoms.append(_make_atom(parts[0], float(parts[2]), float(parts[3]), float(parts[4])))
            except ValueError:
                continue

    return _make_structure(atoms, name=filename, source_format="gamess", source_software="gamess")


def _is_float(value: str) -> bool:
    try:
        float(value)
        return True
    except ValueError:
        return False


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
    coordinate_mode = cell.get("coordinateMode", "cartesian")

    for ia in range(na):
        for ib in range(nb):
            for ic in range(nc):
                for atom in base_atoms:
                    if coordinate_mode == "fractional":
                        super_atoms.append(_make_atom(
                            atom["element"],
                            (atom["x"] + ia) / na,
                            (atom["y"] + ib) / nb,
                            (atom["z"] + ic) / nc,
                        ))
                    else:
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
        coordinate_mode=coordinate_mode,
    )

    return _make_structure(
        super_atoms, cell=super_cell,
        name=f"{structure.get('name', 'supercell')} {na}x{nb}x{nc}",
        source_format=structure.get("metadata", {}).get("source", {}).get("format"),
    )


# ---------------------------------------------------------------------------
# Native structure rendering
# ---------------------------------------------------------------------------

def convert_structure(structure: Dict[str, Any], target_format: str) -> Dict[str, Any]:
    """Convert an OpenQCStructure dict into a native text format."""
    fmt = _normalize_target_format(target_format)
    if fmt == "openqc":
        return {
            "format": "openqc",
            "structure": structure,
            "atomCount": len(structure.get("atoms", [])),
        }

    renderers = {
        "xyz": _render_xyz,
        "pdb": _render_pdb,
        "poscar": _render_poscar,
        "cif": _render_cif,
    }
    if fmt not in renderers:
        raise ValueError(
            "Unsupported native conversion target. Supported targets: xyz, pdb, poscar, vasp, cif, openqc"
        )

    return {
        "format": fmt,
        "content": renderers[fmt](structure),
        "atomCount": len(structure.get("atoms", [])),
    }


def _normalize_target_format(target_format: str) -> str:
    fmt = (target_format or "xyz").strip().lower().lstrip(".")
    aliases = {
        "vasp": "poscar",
        "contcar": "poscar",
        "json": "openqc",
        "openqcstructure": "openqc",
    }
    return aliases.get(fmt, fmt)


def _render_xyz(structure: Dict[str, Any]) -> str:
    atoms = _cartesian_atoms(structure)
    lines = [str(len(atoms)), structure.get("name") or "OpenQC export"]
    lines.extend(
        f"{atom['element']} {_coord(atom['x'])} {_coord(atom['y'])} {_coord(atom['z'])}"
        for atom in atoms
    )
    return "\n".join(lines) + "\n"


def _render_pdb(structure: Dict[str, Any]) -> str:
    lines = []
    for idx, atom in enumerate(_cartesian_atoms(structure), start=1):
        element = str(atom["element"])[:2].rjust(2)
        lines.append(
            "HETATM"
            f"{idx:5d} {element.ljust(4)} MOL     1    "
            f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}"
            f"  1.00  0.00          {element}"
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def _render_poscar(structure: Dict[str, Any]) -> str:
    cell = structure.get("cell")
    if not cell:
        raise ValueError("POSCAR conversion requires a unit cell")

    grouped = _group_atoms_by_element(structure.get("atoms", []))
    coordinate_mode = cell.get("coordinateMode", "cartesian")
    lines = [
        structure.get("name") or "OpenQC export",
        "1.0",
        _vector_line(cell["a"]),
        _vector_line(cell["b"]),
        _vector_line(cell["c"]),
        " ".join(group["element"] for group in grouped),
        " ".join(str(len(group["atoms"])) for group in grouped),
        "Direct" if coordinate_mode == "fractional" else "Cartesian",
    ]

    for group in grouped:
        for atom in group["atoms"]:
            lines.append(f"{_coord(atom['x'])} {_coord(atom['y'])} {_coord(atom['z'])}")

    return "\n".join(lines) + "\n"


def _render_cif(structure: Dict[str, Any]) -> str:
    cell = structure.get("cell")
    if not cell:
        raise ValueError("CIF conversion requires a unit cell")

    a, b, c = cell["a"], cell["b"], cell["c"]
    atoms = _cartesian_atoms(structure)
    lines = [
        "data_openqc_export",
        f"_cell_length_a {_coord(_norm(a))}",
        f"_cell_length_b {_coord(_norm(b))}",
        f"_cell_length_c {_coord(_norm(c))}",
        f"_cell_angle_alpha {_coord(_angle_degrees(b, c))}",
        f"_cell_angle_beta {_coord(_angle_degrees(a, c))}",
        f"_cell_angle_gamma {_coord(_angle_degrees(a, b))}",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_Cartn_x",
        "_atom_site_Cartn_y",
        "_atom_site_Cartn_z",
    ]
    for idx, atom in enumerate(atoms, start=1):
        label = f"{atom['element']}{idx}"
        lines.append(
            f"{label} {atom['element']} {_coord(atom['x'])} {_coord(atom['y'])} {_coord(atom['z'])}"
        )
    return "\n".join(lines) + "\n"


def _cartesian_atoms(structure: Dict[str, Any]) -> List[Dict[str, Any]]:
    atoms = structure.get("atoms", [])
    cell = structure.get("cell")
    if not cell or cell.get("coordinateMode") != "fractional":
        return [dict(atom) for atom in atoms]

    return [
        {
            **atom,
            **dict(zip(("x", "y", "z"), _fractional_to_cartesian(
                [atom["x"], atom["y"], atom["z"]], cell
            ))),
        }
        for atom in atoms
    ]


def _fractional_to_cartesian(frac: List[float], cell: Dict[str, Any]) -> List[float]:
    u, v, w = frac
    a, b, c = cell["a"], cell["b"], cell["c"]
    return [
        u * a[0] + v * b[0] + w * c[0],
        u * a[1] + v * b[1] + w * c[1],
        u * a[2] + v * b[2] + w * c[2],
    ]


def _group_atoms_by_element(atoms: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    groups: List[Dict[str, Any]] = []
    index_by_element: Dict[str, int] = {}
    for atom in atoms:
        element = str(atom["element"])
        if element not in index_by_element:
            index_by_element[element] = len(groups)
            groups.append({"element": element, "atoms": []})
        groups[index_by_element[element]]["atoms"].append(atom)
    return groups


def _vector_line(vector: List[float]) -> str:
    return " ".join(_coord(value) for value in vector)


def _coord(value: float) -> str:
    return f"{float(value):.6f}" if math.isfinite(float(value)) else "0.000000"


def _norm(vector: List[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in vector))


def _angle_degrees(a: List[float], b: List[float]) -> float:
    denom = _norm(a) * _norm(b)
    if denom == 0:
        return 0.0
    cosine = sum(float(a[i]) * float(b[i]) for i in range(3)) / denom
    cosine = min(1.0, max(-1.0, cosine))
    return math.degrees(math.acos(cosine))


# ---------------------------------------------------------------------------
# Auto-detect format and parse
# ---------------------------------------------------------------------------

def _detect_format(filepath: str, content: Optional[str] = None) -> str:
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
        return _detect_inp_format(filepath, content)
    if basename.startswith("INCAR") or basename.startswith("KPOINTS"):
        return "vasp"
    return ext.strip(".") or "unknown"


def _detect_inp_format(filepath: str, content: Optional[str] = None) -> str:
    """Disambiguate common quantum-chemistry .inp dialects by content."""
    if content is None:
        try:
            with open(filepath, "r") as handle:
                content = handle.read(8192)
        except OSError:
            content = ""

    upper = content.upper()
    if re.search(r"^\s*&(?:GLOBAL|FORCE_EVAL|SUBSYS|COORD)\b", upper, re.MULTILINE):
        return "cp2k"
    if re.search(r"^\s*\$(?:CONTRL|DATA|BASIS|SYSTEM|SCF)\b", upper, re.MULTILINE):
        return "gamess"
    return "orca"


def parse_file(filepath: str, format_hint: str = None, content: str = None) -> Dict[str, Any]:
    """Parse a structure file into OpenQCStructure JSON."""
    normalized_hint = (format_hint or "").strip().lower()
    fmt = _detect_format(filepath, content) if normalized_hint in ("", "auto") else normalized_hint
    name = os.path.basename(filepath)

    if content is None and fmt in ("xyz", "poscar", "gaussian", "orca", "cp2k", "gamess"):
        try:
            with open(filepath, "r") as f:
                content = f.read()
        except OSError:
            content = None

    # Native parsers first (no dependency needed)
    if fmt == "xyz" and content:
        return _parse_xyz(content, name)
    if fmt == "poscar" and content:
        return _parse_poscar(content, name)
    if fmt == "gaussian" and content:
        return _parse_gaussian_input(content, name)
    if fmt == "orca" and content:
        return _parse_orca_input(content, name)
    if fmt == "cp2k" and content:
        return _parse_cp2k_input(content, name)
    if fmt == "gamess" and content:
        return _parse_gamess_input(content, name)

    parser_errors = []

    # Try pymatgen first for CIF, where it is usually the best parser.
    if fmt in ("cif",):
        try:
            return _parse_with_pymatgen(filepath, fmt)
        except Exception as exc:
            parser_errors.append(f"pymatgen: {exc}")

    # Try ASE as general fallback
    try:
        return _parse_with_ase(filepath, fmt)
    except Exception as exc:
        parser_errors.append(f"ASE: {exc}")

    # Try pymatgen after ASE for other file-backed formats.
    if fmt not in ("cif",):
        try:
            return _parse_with_pymatgen(filepath, fmt)
        except Exception as exc:
            parser_errors.append(f"pymatgen: {exc}")

    details = f" Attempts: {'; '.join(parser_errors)}" if parser_errors else ""

    raise ValueError(
        f"Cannot parse format '{fmt}'. Install pymatgen or ASE for extended format support: "
        f"pip install pymatgen ase.{details}"
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
            content = args.get("content")
            structure = parse_file(filepath, format_hint=fmt, content=content)
            write_response(convert_structure(structure, target_format))

        else:
            write_error(f"Unknown command: {command}")

    except Exception as e:
        write_error(str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
