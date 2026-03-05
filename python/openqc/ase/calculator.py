"""
ASE Calculator Backend

Provides ASE calculator wrapper for VASP, CP2K, QE.
Handles input generation, execution interface, and result parsing.
"""

from typing import Dict, Any, Optional, List
from pathlib import Path
import json
import sys
import os
import subprocess
import time
from dataclasses import dataclass, asdict
from ase import Atoms
from ase.io import write
import numpy as np


@dataclass
class CalculatorInput:
    """Calculator input parameters."""
    atoms: Dict[str, Any]
    calculator: str
    parameters: Dict[str, Any]
    pseudopotentials: Optional[Dict[str, str]] = None
    basis_set: Optional[str] = None
    charge: Optional[int] = None
    spin_multiplicity: Optional[int] = None
    work_dir: str = "."
    label: Optional[str] = None


@dataclass
class CalculatorConfig:
    """Calculator execution configuration."""
    executable: str
    command: Optional[str] = None
    environment: Optional[Dict[str, str]] = None
    max_memory: Optional[int] = None
    num_cores: Optional[int] = None
    queue: Optional[str] = None
    wall_time: Optional[int] = None


@dataclass
class CalculationResult:
    """Result of a calculation."""
    success: bool
    status: str
    work_dir: str
    input_files: List[str]
    output_files: List[str]
    energy: Optional[float] = None
    forces: Optional[List[List[float]]] = None
    stress: Optional[List[List[float]]] = None
    error: Optional[str] = None
    warnings: Optional[List[str]] = None
    metadata: Optional[Dict[str, Any]] = None
    execution_time: Optional[float] = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []
        if self.metadata is None:
            self.metadata = {}

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to dictionary."""
        return asdict(self)


class ASECalculatorWrapper:
    """
    ASE Calculator wrapper for quantum chemistry codes.
    
    Supports VASP, CP2K, Quantum ESPRESSO, Gaussian, ORCA.
    """
    
    def __init__(self):
        self.supported_calculators = {
            "vasp": self._generate_vasp_input,
            "cp2k": self._generate_cp2k_input,
            "qe": self._generate_qe_input,
            "gaussian": self._generate_gaussian_input,
            "orca": self._generate_orca_input,
        }
    
    def generate_input(self, input_data: CalculatorInput) -> CalculationResult:
        """Generate input files for a calculator."""
        try:
            calculator = input_data.calculator.lower()
            if calculator not in self.supported_calculators:
                return CalculationResult(
                    success=False,
                    status="failed",
                    work_dir=input_data.work_dir,
                    input_files=[],
                    output_files=[],
                    error=f"Unsupported calculator: {calculator}"
                )
            
            # Create work directory
            work_dir = Path(input_data.work_dir)
            work_dir.mkdir(parents=True, exist_ok=True)
            
            # Generate input files
            input_files = self.supported_calculators[calculator](input_data)
            
            return CalculationResult(
                success=True,
                status="completed",
                work_dir=str(work_dir),
                input_files=input_files,
                output_files=[],
                metadata={"calculator": calculator}
            )
            
        except Exception as e:
            return CalculationResult(
                success=False,
                status="failed",
                work_dir=input_data.work_dir,
                input_files=[],
                output_files=[],
                error=str(e)
            )
    
    def run_calculation(self, input_data: CalculatorInput, config: CalculatorConfig) -> CalculationResult:
        """Run a calculation."""
        start_time = time.time()
        
        try:
            # Generate input first
            gen_result = self.generate_input(input_data)
            if not gen_result.success:
                return gen_result
            
            calculator = input_data.calculator.lower()
            work_dir = Path(input_data.work_dir)
            
            # Build command
            if config.command:
                cmd = config.command.split()
            else:
                cmd = [config.executable]
            
            if config.num_cores and config.num_cores > 1:
                cmd = ["mpirun", "-np", str(config.num_cores)] + cmd
            
            # Set environment
            env = os.environ.copy()
            if config.environment:
                env.update(config.environment)
            if config.max_memory:
                env["MEMORY"] = str(config.max_memory)
            
            # Run calculation
            process = subprocess.run(
                cmd,
                cwd=str(work_dir),
                env=env,
                capture_output=True,
                text=True,
                timeout=config.wall_time if config.wall_time else None
            )
            
            execution_time = time.time() - start_time
            
            # Check for output files
            output_files = self._find_output_files(work_dir, calculator)
            
            # Parse results if successful
            energy = None
            forces = None
            stress = None
            
            if process.returncode == 0:
                energy, forces, stress = self._parse_results(work_dir, calculator)
            
            return CalculationResult(
                success=process.returncode == 0,
                status="completed" if process.returncode == 0 else "failed",
                work_dir=str(work_dir),
                input_files=gen_result.input_files,
                output_files=output_files,
                energy=energy,
                forces=forces,
                stress=stress,
                error=process.stderr if process.returncode != 0 else None,
                execution_time=execution_time,
                metadata={"return_code": process.returncode}
            )
            
        except subprocess.TimeoutExpired:
            return CalculationResult(
                success=False,
                status="failed",
                work_dir=input_data.work_dir,
                input_files=[],
                output_files=[],
                error=f"Calculation timed out after {config.wall_time} seconds",
                execution_time=time.time() - start_time
            )
        except Exception as e:
            return CalculationResult(
                success=False,
                status="failed",
                work_dir=input_data.work_dir,
                input_files=[],
                output_files=[],
                error=str(e),
                execution_time=time.time() - start_time
            )
    
    def read_results(self, work_dir: str, calculator: str) -> CalculationResult:
        """Read calculation results from working directory."""
        try:
            work_path = Path(work_dir)
            if not work_path.exists():
                return CalculationResult(
                    success=False,
                    status="failed",
                    work_dir=work_dir,
                    input_files=[],
                    output_files=[],
                    error=f"Working directory not found: {work_dir}"
                )
            
            # Find files
            input_files = list(work_path.glob("*.in")) + list(work_path.glob("*.inp"))
            input_files += list(work_path.glob("POSCAR")) + list(work_path.glob("INCAR"))
            
            output_files = self._find_output_files(work_path, calculator)
            
            # Parse results
            energy, forces, stress = self._parse_results(work_path, calculator)
            
            success = energy is not None or len(output_files) > 0
            
            return CalculationResult(
                success=success,
                status="completed" if success else "failed",
                work_dir=str(work_path),
                input_files=[str(f) for f in input_files],
                output_files=output_files,
                energy=energy,
                forces=forces,
                stress=stress,
                metadata={"calculator": calculator}
            )
            
        except Exception as e:
            return CalculationResult(
                success=False,
                status="failed",
                work_dir=work_dir,
                input_files=[],
                output_files=[],
                error=str(e)
            )
    
    def _generate_vasp_input(self, input_data: CalculatorInput) -> List[str]:
        """Generate VASP input files."""
        work_dir = Path(input_data.work_dir)
        atoms = self._dict_to_atoms(input_data.atoms)
        params = input_data.parameters
        
        # Write POSCAR
        write(str(work_dir / "POSCAR"), atoms, format="vasp")
        
        # Generate INCAR
        incar_content = self._generate_incar(params)
        with open(work_dir / "INCAR", "w") as f:
            f.write(incar_content)
        
        # Generate KPOINTS
        kpts = params.get("kpts", [4, 4, 4])
        kpoints_content = self._generate_kpoints(kpts)
        with open(work_dir / "KPOINTS", "w") as f:
            f.write(kpoints_content)
        
        return ["POSCAR", "INCAR", "KPOINTS"]
    
    def _generate_incar(self, params: Dict[str, Any]) -> str:
        """Generate VASP INCAR content."""
        lines = ["# VASP INCAR generated by OpenQC", ""]
        
        # Basic parameters
        if "xc" in params:
            xc_map = {"PBE": "PBE", "LDA": "CA", "PW91": "PW91"}
            val = params["xc"]
            lines.append(f"GGA = {xc_map.get(val, 'PBE')}")
        
        if "encut" in params:
            lines.append(f"ENCUT = {params['encut']}")
        
        if "ismear" in params:
            lines.append(f"ISMEAR = {params['ismear']}")
        
        if "sigma" in params:
            lines.append(f"SIGMA = {params['sigma']}")
        
        if "ibrion" in params:
            lines.append(f"IBRION = {params['ibrion']}")
        
        if "nsw" in params:
            lines.append(f"NSW = {params['nsw']}")
        
        if "ediff" in params:
            lines.append(f"EDIFF = {params['ediff']:.1e}")
        
        if "ediffg" in params:
            lines.append(f"EDIFFG = {params['ediffg']}")
        
        if "isif" in params:
            lines.append(f"ISIF = {params['isif']}")
        
        return "\n".join(lines)
    
    def _generate_kpoints(self, kpts: List[int]) -> str:
        """Generate VASP KPOINTS content."""
        return f"""Automatic mesh
0
Gamma
{kpts[0]} {kpts[1]} {kpts[2]}
0 0 0
"""
    
    def _generate_cp2k_input(self, input_data: CalculatorInput) -> List[str]:
        """Generate CP2K input file."""
        work_dir = Path(input_data.work_dir)
        atoms = self._dict_to_atoms(input_data.atoms)
        params = input_data.parameters
        label = input_data.label or "cp2k"
        
        run_type = params.get("run_type", "ENERGY")
        basis_file = params.get("basis_set_file", "BASIS_MOLOPT")
        pot_file = params.get("potential_file", "GTH_POTENTIALS")
        xc_func = params.get("xc", "PBE")
        max_scf = params.get("max_scf", 50)
        scf_eps = params.get("scf_eps", 1e-7)
        
        content = f"""&GLOBAL
  PROJECT {label}
  RUN_TYPE {run_type}
  PRINT_LEVEL MEDIUM
&END GLOBAL

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME {basis_file}
    POTENTIAL_FILE_NAME {pot_file}
    &XC
      &XC_FUNCTIONAL {xc_func}
      &END XC_FUNCTIONAL
    &END XC
    &SCF
      MAX_SCF {max_scf}
      EPS_SCF {scf_eps:.1e}
    &END SCF
  &END DFT
  &SUBSYS
    &CELL
"""
        
        # Add cell parameters
        if atoms.cell is not None and len(atoms.cell) > 0:
            content += f"      A {atoms.cell[0, 0]:.6f} {atoms.cell[0, 1]:.6f} {atoms.cell[0, 2]:.6f}\n"
            content += f"      B {atoms.cell[1, 0]:.6f} {atoms.cell[1, 1]:.6f} {atoms.cell[1, 2]:.6f}\n"
            content += f"      C {atoms.cell[2, 0]:.6f} {atoms.cell[2, 1]:.6f} {atoms.cell[2, 2]:.6f}\n"
        
        pbc_str = " ".join(["XYZ"[i] if p else "NONE" for i, p in enumerate(atoms.pbc)])
        content += f"      PERIODIC {pbc_str}\n"
        content += "    &END CELL\n"
        content += "    &COORD\n"
        
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            content += f"      {symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
        
        content += "    &END COORD\n"
        
        # Add kind sections
        basis_set = params.get("basis_set", "DZVP-MOLOPT-GTH")
        potential = params.get("potential", "GTH-PBE")
        unique_elements = set(atoms.get_chemical_symbols())
        for elem in unique_elements:
            content += f"""    &KIND {elem}
      BASIS_SET {basis_set}
      POTENTIAL {potential}
    &END KIND
"""
        
        content += """  &END SUBSYS
&END FORCE_EVAL
"""
        
        inp_file = f"{label}.inp"
        with open(work_dir / inp_file, "w") as f:
            f.write(content)
        
        return [inp_file]
    
    def _generate_qe_input(self, input_data: CalculatorInput) -> List[str]:
        """Generate Quantum ESPRESSO input file."""
        work_dir = Path(input_data.work_dir)
        atoms = self._dict_to_atoms(input_data.atoms)
        params = input_data.parameters
        label = input_data.label or "qe"
        
        calc_type = params.get("calculation", "scf")
        ecutwfc = params.get("ecutwfc", 40)
        ecutrho = params.get("ecutrho", 320)
        occupations = params.get("occupations", "smearing")
        smearing = params.get("smearing", "gaussian")
        degauss = params.get("degauss", 0.02)
        conv_thr = params.get("conv_thr", 1e-8)
        
        content = f"""&CONTROL
  calculation = '{calc_type}'
  prefix = '{label}'
  outdir = './tmp'
  pseudo_dir = './pseudo'
/

&SYSTEM
  ibrav = 0
  nat = {len(atoms)}
  ntyp = {len(set(atoms.get_chemical_symbols()))}
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  occupations = '{occupations}'
  smearing = '{smearing}'
  degauss = {degauss}
/

&ELECTRONS
  conv_thr = {conv_thr:.1e}
/

CELL_PARAMETERS angstrom
{atoms.cell[0, 0]:.10f} {atoms.cell[0, 1]:.10f} {atoms.cell[0, 2]:.10f}
{atoms.cell[1, 0]:.10f} {atoms.cell[1, 1]:.10f} {atoms.cell[1, 2]:.10f}
{atoms.cell[2, 0]:.10f} {atoms.cell[2, 1]:.10f} {atoms.cell[2, 2]:.10f}

ATOMIC_SPECIES
"""
        
        # Add atomic species
        for elem in set(atoms.get_chemical_symbols()):
            if input_data.pseudopotentials:
                pseudo = input_data.pseudopotentials.get(elem, f"{elem}.UPF")
            else:
                pseudo = f"{elem}.UPF"
            content += f"  {elem} 1.0 {pseudo}\n"
        
        content += "\nATOMIC_POSITIONS angstrom\n"
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            content += f"  {symbol} {pos[0]:.10f} {pos[1]:.10f} {pos[2]:.10f}\n"
        
        # K-points
        kpts = params.get("kpts", [4, 4, 4])
        content += f"\nK_POINTS automatic\n{kpts[0]} {kpts[1]} {kpts[2]} 0 0 0\n"
        
        inp_file = f"{label}.in"
        with open(work_dir / inp_file, "w") as f:
            f.write(content)
        
        return [inp_file]
    
    def _generate_gaussian_input(self, input_data: CalculatorInput) -> List[str]:
        """Generate Gaussian input file."""
        work_dir = Path(input_data.work_dir)
        atoms = self._dict_to_atoms(input_data.atoms)
        params = input_data.parameters
        label = input_data.label or "gaussian"
        
        charge = input_data.charge or 0
        mult = input_data.spin_multiplicity or 1
        
        method = params.get("method", "B3LYP")
        basis = params.get("basis", "6-31G(d)")
        
        route = f"#{method}/{basis}"
        if params.get("opt"):
            route += " opt"
        if params.get("scf"):
            route += f" scf={params['scf']}"
        
        content = f"{route}\n\n{label}\n\n{charge} {mult}\n"
        
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            content += f"{symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
        
        content += "\n"
        
        com_file = f"{label}.com"
        with open(work_dir / com_file, "w") as f:
            f.write(content)
        
        return [com_file]
    
    def _generate_orca_input(self, input_data: CalculatorInput) -> List[str]:
        """Generate ORCA input file."""
        work_dir = Path(input_data.work_dir)
        atoms = self._dict_to_atoms(input_data.atoms)
        params = input_data.parameters
        label = input_data.label or "orca"
        
        charge = input_data.charge or 0
        mult = input_data.spin_multiplicity or 1
        
        method = params.get("method", "B3LYP")
        basis = params.get("basis", "def2-SVP")
        scf_conv = params.get("scfconvergence", "TightSCF")
        
        content = f"! {method} {basis} {scf_conv}\n\n* xyz {charge} {mult}\n"
        
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
            content += f"  {symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
        
        content += "*\n"
        
        inp_file = f"{label}.inp"
        with open(work_dir / inp_file, "w") as f:
            f.write(content)
        
        return [inp_file]
    
    def _dict_to_atoms(self, atoms_dict: Dict[str, Any]) -> Atoms:
        """Convert dictionary to ASE Atoms."""
        positions = atoms_dict.get("positions", [])
        symbols = atoms_dict.get("chemical_symbols", [])
        cell = atoms_dict.get("cell")
        pbc = atoms_dict.get("pbc", [False, False, False])
        
        atoms = Atoms(symbols=symbols, positions=positions, pbc=pbc)
        if cell:
            atoms.set_cell(cell)
        
        return atoms
    
    def _find_output_files(self, work_dir: Path, calculator: str) -> List[str]:
        """Find output files in working directory."""
        output_files = []
        
        if calculator == "vasp":
            for f in ["OUTCAR", "vasprun.xml", "OSZICAR"]:
                if (work_dir / f).exists():
                    output_files.append(str(work_dir / f))
        elif calculator in ["cp2k", "qe", "orca"]:
            for f in work_dir.glob("*.out"):
                output_files.append(str(f))
        elif calculator == "gaussian":
            for f in work_dir.glob("*.log"):
                output_files.append(str(f))
        
        return output_files
    
    def _parse_results(self, work_dir: Path, calculator: str):
        """Parse calculation results from output files."""
        energy = None
        forces = None
        stress = None
        
        try:
            if calculator == "vasp":
                if (work_dir / "OUTCAR").exists():
                    energy, forces = self._parse_vasp_outcar(work_dir / "OUTCAR")
            elif calculator == "cp2k":
                for out_file in work_dir.glob("*.out"):
                    energy = self._parse_cp2k_output(out_file)
                    if energy:
                        break
            elif calculator == "qe":
                for out_file in work_dir.glob("*.out"):
                    energy = self._parse_qe_output(out_file)
                    if energy:
                        break
            elif calculator == "gaussian":
                for log_file in work_dir.glob("*.log"):
                    energy = self._parse_gaussian_output(log_file)
                    if energy:
                        break
            elif calculator == "orca":
                for out_file in work_dir.glob("*.out"):
                    energy = self._parse_orca_output(out_file)
                    if energy:
                        break
        except Exception:
            pass
        
        return energy, forces, stress
    
    def _parse_vasp_outcar(self, outcar_path: Path):
        """Parse energy and forces from VASP OUTCAR."""
        energy = None
        forces = None
        
        with open(outcar_path, "r") as f:
            lines = f.readlines()
        
        for i, line in enumerate(lines):
            if "free  energy   TOTEN" in line:
                try:
                    energy = float(line.split()[-2])
                except (ValueError, IndexError):
                    pass
            
            if "TOTAL-FORCE" in line and i + 1 < len(lines):
                forces = []
                for j in range(i + 2, len(lines)):
                    if "---" in lines[j] or lines[j].strip() == "":
                        break
                    parts = lines[j].split()
                    if len(parts) >= 6:
                        forces.append([float(parts[3]), float(parts[4]), float(parts[5])])
        
        return energy, forces
    
    def _parse_cp2k_output(self, out_file: Path):
        """Parse energy from CP2K output."""
        with open(out_file, "r") as f:
            for line in f:
                if "ENERGY|" in line:
                    try:
                        return float(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
        return None
    
    def _parse_qe_output(self, out_file: Path):
        """Parse energy from Quantum ESPRESSO output."""
        with open(out_file, "r") as f:
            for line in f:
                if "!    total energy" in line:
                    try:
                        return float(line.split()[-2]) * 13.605698  # Ry to eV
                    except (ValueError, IndexError):
                        pass
        return None
    
    def _parse_gaussian_output(self, log_file: Path):
        """Parse energy from Gaussian output."""
        with open(log_file, "r") as f:
            for line in f:
                if "SCF Done" in line:
                    try:
                        return float(line.split()[4])
                    except (ValueError, IndexError):
                        pass
        return None
    
    def _parse_orca_output(self, out_file: Path):
        """Parse energy from ORCA output."""
        with open(out_file, "r") as f:
            for line in f:
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        return float(line.split()[-1])
                    except (ValueError, IndexError):
                        pass
        return None


def main():
    """CLI interface for ASE calculator."""
    if len(sys.argv) < 2:
        print("Usage: python calculator.py <command> [args]")
        print("Commands:")
        print("  generate <input_json>")
        print("  run <input_json> <config_json>")
        print("  read <work_dir> <calculator>")
        sys.exit(1)
    
    command = sys.argv[1]
    calculator = ASECalculatorWrapper()
    
    if command == "generate":
        if len(sys.argv) < 3:
            print("Error: input JSON required")
            sys.exit(1)
        
        input_data = CalculatorInput(**json.loads(sys.argv[2]))
        result = calculator.generate_input(input_data)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "run":
        if len(sys.argv) < 4:
            print("Error: input and config JSON required")
            sys.exit(1)
        
        input_data = CalculatorInput(**json.loads(sys.argv[2]))
        config = CalculatorConfig(**json.loads(sys.argv[3]))
        result = calculator.run_calculation(input_data, config)
        print(json.dumps(result.to_dict(), indent=2))
    
    elif command == "read":
        if len(sys.argv) < 4:
            print("Error: work_dir and calculator required")
            sys.exit(1)
        
        result = calculator.read_results(sys.argv[2], sys.argv[3])
        print(json.dumps(result.to_dict(), indent=2))
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)


if __name__ == "__main__":
    main()
