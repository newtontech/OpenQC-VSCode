"""
ASE Complex Property Mapping

Backend for complex property mapping including Hubbard U parameters,
constraints, excited state methods, and pseudopotential strategies.
"""

from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path
import json
import sys
from dataclasses import dataclass, asdict
from ase import Atoms
import numpy as np


@dataclass
class HubbardParameters:
    """Hubbard U parameters for DFT+U calculations."""
    enabled: bool
    method: str  # e.g., "Dudarev", "Liechtenstein"
    u_values: Dict[str, float]  # element -> U value in eV
    j_values: Optional[Dict[str, float]] = None  # element -> J value in eV
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ConstraintInfo:
    """Information about atomic constraints."""
    enabled: bool
    constrained_atoms: List[int]  # atom indices
    constraint_type: str  # e.g., "fixed", "line", "plane"
    constraint_vectors: Optional[List[List[float]]] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ExcitedStateInfo:
    """Information about excited state methods."""
    enabled: bool
    method: str  # e.g., "TD-DFT", "GW", "BSE"
    n_states: Optional[int] = None
    n_occ: Optional[int] = None
    n_virt: Optional[int] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PseudopotentialInfo:
    """Pseudopotential/basis set information."""
    strategy: str  # e.g., "PAW", "NORMCONSERVING", "GTH"
    functional: str  # e.g., "PBE", "LDA"
    element_mapping: Dict[str, str]  # element -> PP name
    basis_set_mapping: Optional[Dict[str, str]] = None
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class ComplexPropertyMapper:
    """
    Maps complex properties between quantum chemistry codes.
    """
    
    def __init__(self):
        """Initialize the property mapper."""
        self.hubbard_mappings = {
            ("vasp", "vasp"): self._identity_hubbard,
            ("vasp", "qe"): self._vasp_to_qe_hubbard,
            ("vasp", "cp2k"): self._vasp_to_cp2k_hubbard,
            ("qe", "vasp"): self._qe_to_vasp_hubbard,
            ("qe", "qe"): self._identity_hubbard,
            ("qe", "cp2k"): self._qe_to_cp2k_hubbard,
            ("cp2k", "vasp"): self._cp2k_to_vasp_hubbard,
            ("cp2k", "qe"): self._cp2k_to_qe_hubbard,
            ("cp2k", "cp2k"): self._identity_hubbard,
        }
        
        self.pp_mappings = {
            ("vasp", "vasp"): self._identity_pp,
            ("vasp", "qe"): self._vasp_to_qe_pp,
            ("vasp", "cp2k"): self._vasp_to_cp2k_pp,
            ("qe", "vasp"): self._qe_to_vasp_pp,
            ("qe", "qe"): self._identity_pp,
            ("qe", "cp2k"): self._qe_to_cp2k_pp,
            ("cp2k", "vasp"): self._cp2k_to_vasp_pp,
            ("cp2k", "qe"): self._cp2k_to_qe_pp,
            ("cp2k", "cp2k"): self._identity_pp,
        }
    
    # Hubbard U parameter conversions
    def convert_hubbard(
        self,
        hubbard: HubbardParameters,
        source_code: str,
        target_code: str
    ) -> HubbardParameters:
        """
        Convert Hubbard U parameters between codes.
        
        Parameters
        ----------
        hubbard : HubbardParameters
            Source Hubbard parameters
        source_code : str
            Source code (vasp, qe, cp2k)
        target_code : str
            Target code (vasp, qe, cp2k)
        
        Returns
        -------
        HubbardParameters
            Converted Hubbard parameters
        """
        key = (source_code.lower(), target_code.lower())
        if key in self.hubbard_mappings:
            return self.hubbard_mappings[key](hubbard)
        return self._identity_hubbard(hubbard)
    
    def _identity_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Identity mapping for Hubbard parameters."""
        return HubbardParameters(
            enabled=hubbard.enabled,
            method=hubbard.method,
            u_values=hubbard.u_values.copy(),
            j_values=hubbard.j_values.copy() if hubbard.j_values else None
        )
    
    def _vasp_to_qe_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert VASP Hubbard to QE format."""
        # VASP uses LDAU, QE uses Hubbard U with different conventions
        method_map = {
            "Dudarev": "simplified",
            "Liechtenstein": "full"
        }
        return HubbardParameters(
            enabled=hubbard.enabled,
            method=method_map.get(hubbard.method, "simplified"),
            u_values=hubbard.u_values.copy(),
            j_values=hubbard.j_values.copy() if hubbard.j_values else None
        )
    
    def _vasp_to_cp2k_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert VASP Hubbard to CP2K format."""
        return HubbardParameters(
            enabled=hubbard.enabled,
            method="DUDAREV",  # CP2K uses this spelling
            u_values=hubbard.u_values.copy(),
            j_values=None  # CP2K simplified only
        )
    
    def _qe_to_vasp_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert QE Hubbard to VASP format."""
        method_map = {
            "simplified": "Dudarev",
            "full": "Liechtenstein"
        }
        return HubbardParameters(
            enabled=hubbard.enabled,
            method=method_map.get(hubbard.method, "Dudarev"),
            u_values=hubbard.u_values.copy(),
            j_values=hubbard.j_values.copy() if hubbard.j_values else None
        )
    
    def _qe_to_cp2k_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert QE Hubbard to CP2K format."""
        return HubbardParameters(
            enabled=hubbard.enabled,
            method="DUDAREV",
            u_values=hubbard.u_values.copy(),
            j_values=None
        )
    
    def _cp2k_to_vasp_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert CP2K Hubbard to VASP format."""
        return HubbardParameters(
            enabled=hubbard.enabled,
            method="Dudarev",
            u_values=hubbard.u_values.copy(),
            j_values=None
        )
    
    def _cp2k_to_qe_hubbard(self, hubbard: HubbardParameters) -> HubbardParameters:
        """Convert CP2K Hubbard to QE format."""
        return HubbardParameters(
            enabled=hubbard.enabled,
            method="simplified",
            u_values=hubbard.u_values.copy(),
            j_values=None
        )
    
    # Pseudopotential conversions
    def convert_pseudopotential(
        self,
        pp_info: PseudopotentialInfo,
        source_code: str,
        target_code: str
    ) -> PseudopotentialInfo:
        """
        Convert pseudopotential information between codes.
        
        Parameters
        ----------
        pp_info : PseudopotentialInfo
            Source pseudopotential info
        source_code : str
            Source code
        target_code : str
            Target code
        
        Returns
        -------
        PseudopotentialInfo
            Converted pseudopotential info
        """
        key = (source_code.lower(), target_code.lower())
        if key in self.pp_mappings:
            return self.pp_mappings[key](pp_info)
        return self._identity_pp(pp_info)
    
    def _identity_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Identity mapping for pseudopotentials."""
        return PseudopotentialInfo(
            strategy=pp_info.strategy,
            functional=pp_info.functional,
            element_mapping=pp_info.element_mapping.copy(),
            basis_set_mapping=pp_info.basis_set_mapping.copy() if pp_info.basis_set_mapping else None
        )
    
    def _vasp_to_qe_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert VASP POTCAR to QE UPF."""
        new_mapping = {}
        for elem, pp_name in pp_info.element_mapping.items():
            # Convert POTCAR naming to UPF naming
            new_mapping[elem] = f"{elem}.{pp_info.functional}.UPF"
        
        return PseudopotentialInfo(
            strategy="PAW" if "PAW" in pp_info.strategy else "NORMCONSERVING",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping=None
        )
    
    def _vasp_to_cp2k_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert VASP POTCAR to CP2K GTH."""
        new_mapping = {}
        for elem in pp_info.element_mapping:
            new_mapping[elem] = f"GTH-{pp_info.functional}"
        
        return PseudopotentialInfo(
            strategy="GTH",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping={elem: "DZVP-MOLOPT-GTH" for elem in pp_info.element_mapping}
        )
    
    def _qe_to_vasp_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert QE UPF to VASP POTCAR."""
        new_mapping = {}
        for elem in pp_info.element_mapping:
            new_mapping[elem] = elem  # VASP uses element names directly
        
        return PseudopotentialInfo(
            strategy="PAW",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping=None
        )
    
    def _qe_to_cp2k_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert QE UPF to CP2K GTH."""
        new_mapping = {}
        for elem in pp_info.element_mapping:
            new_mapping[elem] = f"GTH-{pp_info.functional}"
        
        return PseudopotentialInfo(
            strategy="GTH",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping={elem: "DZVP-MOLOPT-GTH" for elem in pp_info.element_mapping}
        )
    
    def _cp2k_to_vasp_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert CP2K GTH to VASP POTCAR."""
        new_mapping = {}
        for elem in pp_info.element_mapping:
            new_mapping[elem] = elem
        
        return PseudopotentialInfo(
            strategy="PAW",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping=None
        )
    
    def _cp2k_to_qe_pp(self, pp_info: PseudopotentialInfo) -> PseudopotentialInfo:
        """Convert CP2K GTH to QE UPF."""
        new_mapping = {}
        for elem in pp_info.element_mapping:
            new_mapping[elem] = f"{elem}.{pp_info.functional}.UPF"
        
        return PseudopotentialInfo(
            strategy="NORMCONSERVING",
            functional=pp_info.functional,
            element_mapping=new_mapping,
            basis_set_mapping=None
        )
    
    # Constraint handling
    def convert_constraints(
        self,
        constraints: ConstraintInfo,
        source_code: str,
        target_code: str
    ) -> ConstraintInfo:
        """
        Convert constraint information between codes.
        
        Parameters
        ----------
        constraints : ConstraintInfo
            Source constraint info
        source_code : str
            Source code
        target_code : str
            Target code
        
        Returns
        -------
        ConstraintInfo
            Converted constraint info
        """
        # Constraints are largely code-specific, but we can translate the concept
        if not constraints.enabled:
            return ConstraintInfo(
                enabled=False,
                constrained_atoms=[],
                constraint_type="none"
            )
        
        return ConstraintInfo(
            enabled=constraints.enabled,
            constrained_atoms=constraints.constrained_atoms.copy(),
            constraint_type=constraints.constraint_type,
            constraint_vectors=constraints.constraint_vectors.copy() if constraints.constraint_vectors else None
        )
    
    # Excited state method conversions
    def convert_excited_state(
        self,
        excited: ExcitedStateInfo,
        source_code: str,
        target_code: str
    ) -> ExcitedStateInfo:
        """
        Convert excited state method between codes.
        
        Parameters
        ----------
        excited : ExcitedStateInfo
            Source excited state info
        source_code : str
            Source code
        target_code : str
            Target code
        
        Returns
        -------
        ExcitedStateInfo
            Converted excited state info
        """
        if not excited.enabled:
            return ExcitedStateInfo(
                enabled=False,
                method="none"
            )
        
        # Method mapping between codes
        method_map = {
            ("TD-DFT", "vasp"): "TD-DFT",
            ("TD-DFT", "qe"): "turboTDDFT",
            ("TD-DFT", "cp2k"): "TDDFPT",
            ("GW", "vasp"): "GW",
            ("GW", "qe"): "gw",
            ("GW", "cp2k"): "GW",
        }
        
        key = (excited.method, target_code.lower())
        new_method = method_map.get(key, excited.method)
        
        return ExcitedStateInfo(
            enabled=excited.enabled,
            method=new_method,
            n_states=excited.n_states,
            n_occ=excited.n_occ,
            n_virt=excited.n_virt
        )
    
    # Generate code-specific input sections
    def generate_vasp_hubbard_section(self, hubbard: HubbardParameters) -> str:
        """Generate VASP LDAU section for INCAR."""
        if not hubbard.enabled:
            return ""
        
        lines = ["# Hubbard U parameters", "LDAU = .TRUE."]
        
        # LDAUTYPE: 2=Dudarev, 3=Liechtenstein
        ldau_type = 2 if hubbard.method.lower() == "dudarev" else 3
        lines.append(f"LDAUTYPE = {ldau_type}")
        
        # LDAUL: angular momentum for each element
        # Simplified: assume d-orbitals for transition metals
        ldaul_map = {"Fe": 2, "Co": 2, "Ni": 2, "Cu": 2, "Mn": 2, "Cr": 2, "V": 2, "Ti": 2}
        ldaul = [ldaul_map.get(e, -1) for e in hubbard.u_values.keys()]
        lines.append(f"LDAUL = {' '.join(map(str, ldaul))}")
        
        # LDAUU: U values
        ldauu = list(hubbard.u_values.values())
        lines.append(f"LDAUU = {' '.join(map(str, ldauu))}")
        
        # LDAUJ: J values
        if hubbard.j_values:
            ldauj = [hubbard.j_values.get(e, 0.0) for e in hubbard.u_values.keys()]
            lines.append(f"LDAUJ = {' '.join(map(str, ldauj))}")
        
        return "\n".join(lines)
    
    def generate_qe_hubbard_section(self, hubbard: HubbardParameters) -> str:
        """Generate QE Hubbard section for input."""
        if not hubbard.enabled:
            return ""
        
        lines = ["# Hubbard U parameters", "lda_plus_u = .true."]
        
        method = "simplified" if hubbard.method.lower() in ["dudarev", "simplified"] else "full"
        lines.append(f"lda_plus_u_kind = '{method}'")
        
        for elem, u_val in hubbard.u_values.items():
            lines.append(f"Hubbard_U({elem}) = {u_val}")
        
        if hubbard.j_values:
            for elem, j_val in hubbard.j_values.items():
                if elem in hubbard.u_values:
                    lines.append(f"Hubbard_J({elem}) = {j_val}")
        
        return "\n".join(lines)
    
    def generate_cp2k_hubbard_section(self, hubbard: HubbardParameters) -> str:
        """Generate CP2K Hubbard section for input."""
        if not hubbard.enabled:
            return ""
        
        lines = ["    &DFT+U", "      DFT_PLUS_U_TREATMENT DUDAREV"]
        
        for elem, u_val in hubbard.u_values.items():
            lines.append(f"      U_LOCAL {elem} {u_val}")
        
        lines.append("    &END DFT+U")
        
        return "\n".join(lines)


def main():
    """CLI interface for complex property mapping."""
    if len(sys.argv) < 2:
        print("Usage: python complex_properties.py <command> [args]")
        print("Commands:")
        print("  hubbard <source_code> <target_code> <hubbard_json>")
        print("  pp <source_code> <target_code> <pp_json>")
        print("  constraints <source_code> <target_code> <constraints_json>")
        print("  excited <source_code> <target_code> <excited_json>")
        print("  generate_vasp_hubbard <hubbard_json>")
        print("  generate_qe_hubbard <hubbard_json>")
        print("  generate_cp2k_hubbard <hubbard_json>")
        sys.exit(1)
    
    command = sys.argv[1]
    mapper = ComplexPropertyMapper()
    
    if command == "hubbard":
        if len(sys.argv) < 5:
            print("Error: source_code, target_code, and hubbard_json required")
            sys.exit(1)
        
        hubbard = HubbardParameters(**json.loads(sys.argv[4]))
        result = mapper.convert_hubbard(sys.argv[2], sys.argv[3], hubbard)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "pp":
        if len(sys.argv) < 5:
            print("Error: source_code, target_code, and pp_json required")
            sys.exit(1)
        
        pp_info = PseudopotentialInfo(**json.loads(sys.argv[4]))
        result = mapper.convert_pseudopotential(sys.argv[2], sys.argv[3], pp_info)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "constraints":
        if len(sys.argv) < 5:
            print("Error: source_code, target_code, and constraints_json required")
            sys.exit(1)
        
        constraints = ConstraintInfo(**json.loads(sys.argv[4]))
        result = mapper.convert_constraints(sys.argv[2], sys.argv[3], constraints)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "excited":
        if len(sys.argv) < 5:
            print("Error: source_code, target_code, and excited_json required")
            sys.exit(1)
        
        excited = ExcitedStateInfo(**json.loads(sys.argv[4]))
        result = mapper.convert_excited_state(sys.argv[2], sys.argv[3], excited)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "generate_vasp_hubbard":
        if len(sys.argv) < 3:
            print("Error: hubbard_json required")
            sys.exit(1)
        
        hubbard = HubbardParameters(**json.loads(sys.argv[2]))
        print(mapper.generate_vasp_hubbard_section(hubbard))
    
    elif command == "generate_qe_hubbard":
        if len(sys.argv) < 3:
            print("Error: hubbard_json required")
            sys.exit(1)
        
        hubbard = HubbardParameters(**json.loads(sys.argv[2]))
        print(mapper.generate_qe_hubbard_section(hubbard))
    
    elif command == "generate_cp2k_hubbard":
        if len(sys.argv) < 3:
            print("Error: hubbard_json required")
            sys.exit(1)
        
        hubbard = HubbardParameters(**json.loads(sys.argv[2]))
        print(mapper.generate_cp2k_hubbard_section(hubbard))
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
