"""
ASE Workflow Migration - MD/Optimization Parameter Conversion

Provides tools for migrating molecular dynamics and optimization workflows
between different quantum chemistry codes. Handles conversion of MD parameters,
optimization criteria, and code-specific settings.
"""

from typing import Dict, Any, Optional, List
from pathlib import Path
import json
from dataclasses import dataclass, asdict, field
from enum import Enum


class WorkflowType(Enum):
    """Types of computational workflows."""
    GEOMETRY_OPTIMIZATION = "geometry_optimization"
    MOLECULAR_DYNAMICS = "molecular_dynamics"
    FREQUENCY_CALCULATION = "frequency_calculation"
    TRANSITION_STATE = "transition_state"
    NEB = "nudged_elastic_band"


class ThermostatType(Enum):
    """Supported thermostat types."""
    NONE = "none"
    NOSE_HOOVER = "nose_hoover"
    LANGEVIN = "langevin"
    BERENDSEN = "berendsen"
    VELOCITY_SCALING = "velocity_scaling"


def _convert_enum_to_value(obj):
    """Convert Enum objects to their values for JSON serialization."""
    if isinstance(obj, Enum):
        return obj.value
    elif isinstance(obj, dict):
        return {k: _convert_enum_to_value(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [_convert_enum_to_value(item) for item in obj]
    else:
        return obj


@dataclass
class MDParameters:
    """Molecular dynamics parameters."""
    temperature: Optional[float] = None  # Kelvin
    pressure: Optional[float] = None  # GPa
    timestep: Optional[float] = None  # fs
    n_steps: Optional[int] = None
    thermostat: ThermostatType = ThermostatType.NONE
    thermostat_params: Optional[Dict[str, Any]] = None
    barostat: Optional[str] = None
    barostat_params: Optional[Dict[str, Any]] = None
    initial_velocities: Optional[str] = None
    velocity_temperature: Optional[float] = None


@dataclass
class OptimizationParameters:
    """Geometry optimization parameters."""
    convergence_energy: Optional[float] = None
    convergence_force: Optional[float] = None
    convergence_stress: Optional[float] = None
    max_steps: Optional[int] = None
    optimizer: Optional[str] = None
    trust_radius: Optional[float] = None
    line_search: Optional[str] = None
    constraints: Optional[List[str]] = None


@dataclass
class WorkflowParameters:
    """Complete workflow parameters."""
    workflow_type: WorkflowType
    md_params: Optional[MDParameters] = None
    opt_params: Optional[OptimizationParameters] = None
    code_specific: Optional[Dict[str, Any]] = None


class WorkflowMigrationResult:
    """Result of workflow migration."""
    
    def __init__(
        self,
        success: bool,
        data: Optional[Dict[str, Any]] = None,
        warnings: List[str] = None,
        metadata: Dict[str, Any] = None
    ):
        self.success = success
        self.data = data
        self.warnings = warnings or []
        self.metadata = metadata or {}
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'success': self.success,
            'data': self.data,
            'warnings': self.warnings,
            'metadata': self.metadata
        }


class WorkflowMigrator:
    """
    Migrate MD/Optimization workflows between quantum chemistry codes.
    """
    
    OPTIMIZER_MAPPING = {
        'vasp': {'cg': 'CG', 'rmm-diis': 'RMM-DIIS'},
        'cp2k': {'bfgs': 'BFGS', 'cg': 'CG', 'lbfgs': 'L-BFGS'},
        'qe': {'bfgs': 'bfgs', 'cg': 'cg', 'damp': 'damped-dynamics'},
        'gaussian': {'bfgs': 'CalcFC', 'rfo': 'RFO', 'gdiis': 'GDIIS'},
        'orca': {'bfgs': 'BFGS', 'nr': 'Newton-Raphson', 'cg': 'CG'},
        'lammps': {'cg': 'min_style cg', 'sd': 'min_style sd', 'fire': 'min_style fire'}
    }
    
    DEFAULT_CONVERGENCE = {
        'vasp': {'energy': 1e-4, 'force': 0.02, 'stress': 0.5},
        'cp2k': {'energy': 1e-3, 'force': 0.003, 'stress': 0.1},
        'qe': {'energy': 1e-4, 'force': 0.001, 'stress': 0.5},
        'gaussian': {'energy': 1e-6, 'force': 0.00045},
        'orca': {'energy': 5e-6, 'force': 0.002},
        'lammps': {'energy': 1e-4, 'force': 1e-6}
    }
    
    THERMOSTAT_MAPPING = {
        'vasp': {
            ThermostatType.NOSE_HOOVER: 'SMASS',
            ThermostatType.VELOCITY_SCALING: 'SCALEE'
        },
        'cp2k': {
            ThermostatType.NOSE_HOOVER: 'NOSE',
            ThermostatType.LANGEVIN: 'LANGEVIN'
        },
        'qe': {
            ThermostatType.VELOCITY_SCALING: 'velocity_scaling',
            ThermostatType.LANGEVIN: 'langevin',
            ThermostatType.NOSE_HOOVER: 'nose'
        },
        'lammps': {
            ThermostatType.NOSE_HOOVER: 'fix nvt',
            ThermostatType.LANGEVIN: 'fix langevin',
            ThermostatType.BERENDSEN: 'fix temp/berendsen'
        }
    }
    
    def __init__(self):
        self.supported_codes = ['vasp', 'cp2k', 'qe', 'gaussian', 'orca', 'lammps']
    
    def migrate_workflow(
        self,
        source_code: str,
        target_code: str,
        workflow_params: WorkflowParameters
    ) -> WorkflowMigrationResult:
        """Migrate workflow parameters from source code to target code."""
        source_code = source_code.lower()
        target_code = target_code.lower()
        
        if source_code not in self.supported_codes:
            return WorkflowMigrationResult(
                success=False,
                warnings=[f"Unsupported source code: {source_code}"]
            )
        
        if target_code not in self.supported_codes:
            return WorkflowMigrationResult(
                success=False,
                warnings=[f"Unsupported target code: {target_code}"]
            )
        
        warnings = []
        target_params = {}
        
        if workflow_params.workflow_type == WorkflowType.MOLECULAR_DYNAMICS:
            if workflow_params.md_params is None:
                workflow_params.md_params = MDParameters()
            target_params = self._migrate_md(
                source_code,
                target_code,
                workflow_params.md_params,
                warnings
            )
        elif workflow_params.workflow_type == WorkflowType.GEOMETRY_OPTIMIZATION:
            if workflow_params.opt_params is None:
                workflow_params.opt_params = OptimizationParameters()
            target_params = self._migrate_optimization(
                source_code,
                target_code,
                workflow_params.opt_params,
                warnings
            )
        
        code_notes = self._get_code_specific_notes(target_code, workflow_params)
        if code_notes:
            warnings.append(code_notes)
        
        return WorkflowMigrationResult(
            success=True,
            data=target_params,
            warnings=warnings,
            metadata={
                'source_code': source_code,
                'target_code': target_code,
                'workflow_type': workflow_params.workflow_type.value
            }
        )
    
    def _migrate_md(
        self,
        source_code: str,
        target_code: str,
        md_params: MDParameters,
        warnings: List[str]
    ) -> Dict[str, Any]:
        """Migrate MD parameters."""
        target_params = {}
        
        if md_params.temperature is not None:
            target_params['temperature'] = md_params.temperature
        
        if md_params.pressure is not None:
            target_params['pressure'] = md_params.pressure
            if target_code in ['gaussian', 'orca']:
                warnings.append(
                    f"{target_code.upper()} does not support NPT ensemble natively."
                )
        
        if md_params.timestep is not None:
            target_params['timestep'] = md_params.timestep
        
        if md_params.n_steps is not None:
            target_params['n_steps'] = md_params.n_steps
        
        if md_params.thermostat != ThermostatType.NONE:
            thermostat_map = self.THERMOSTAT_MAPPING.get(target_code, {})
            target_thermostat = thermostat_map.get(md_params.thermostat)
            
            if target_thermostat:
                target_params['thermostat'] = target_thermostat
            else:
                warnings.append(
                    f"Thermostat {md_params.thermostat.value} not directly supported in {target_code.upper()}."
                )
                available = list(thermostat_map.values())
                if available:
                    warnings.append(f"Available thermostats in {target_code.upper()}: {', '.join(available)}")
        
        if target_code == 'vasp':
            target_params.update(self._get_vasp_md_params(md_params))
        elif target_code == 'cp2k':
            target_params.update(self._get_cp2k_md_params(md_params))
        elif target_code == 'qe':
            target_params.update(self._get_qe_md_params(md_params))
        elif target_code == 'lammps':
            target_params.update(self._get_lammps_md_params(md_params))
        
        return target_params
    
    def _migrate_optimization(
        self,
        source_code: str,
        target_code: str,
        opt_params: OptimizationParameters,
        warnings: List[str]
    ) -> Dict[str, Any]:
        """Migrate optimization parameters."""
        target_params = {}
        
        if opt_params.convergence_energy is not None:
            target_params['convergence_energy'] = opt_params.convergence_energy
        else:
            defaults = self.DEFAULT_CONVERGENCE.get(target_code, {})
            if 'energy' in defaults:
                target_params['convergence_energy'] = defaults['energy']
        
        if opt_params.convergence_force is not None:
            target_params['convergence_force'] = opt_params.convergence_force
        else:
            defaults = self.DEFAULT_CONVERGENCE.get(target_code, {})
            if 'force' in defaults:
                target_params['convergence_force'] = defaults['force']
        
        if opt_params.convergence_stress is not None:
            target_params['convergence_stress'] = opt_params.convergence_stress
            if target_code in ['gaussian', 'orca']:
                warnings.append(
                    f"{target_code.upper()} does not support stress convergence for molecular systems."
                )
        
        if opt_params.max_steps is not None:
            target_params['max_steps'] = opt_params.max_steps
        
        if opt_params.optimizer is not None:
            optimizer_map = self.OPTIMIZER_MAPPING.get(target_code, {})
            target_optimizer = optimizer_map.get(opt_params.optimizer.lower())
            
            if target_optimizer:
                target_params['optimizer'] = target_optimizer
            else:
                warnings.append(
                    f"Optimizer {opt_params.optimizer} not directly supported in {target_code.upper()}."
                )
        
        if target_code == 'vasp':
            target_params.update(self._get_vasp_opt_params(opt_params))
        elif target_code == 'cp2k':
            target_params.update(self._get_cp2k_opt_params(opt_params))
        elif target_code == 'qe':
            target_params.update(self._get_qe_opt_params(opt_params))
        elif target_code == 'gaussian':
            target_params.update(self._get_gaussian_opt_params(opt_params))
        elif target_code == 'orca':
            target_params.update(self._get_orca_opt_params(opt_params))
        
        return target_params
    
    def _get_vasp_md_params(self, md_params: MDParameters) -> Dict[str, Any]:
        params = {}
        params['ibrion'] = 0
        params['potim'] = md_params.timestep if md_params.timestep else 1.0
        params['nsw'] = md_params.n_steps if md_params.n_steps else 100
        
        if md_params.thermostat == ThermostatType.NOSE_HOOVER:
            params['smaas'] = 0
        elif md_params.thermostat == ThermostatType.VELOCITY_SCALING:
            params['smaas'] = -1
        
        if md_params.temperature:
            params['tebeg'] = md_params.temperature
            params['teend'] = md_params.temperature
        
        return params
    
    def _get_cp2k_md_params(self, md_params: MDParameters) -> Dict[str, Any]:
        params = {}
        params['run_type'] = 'MD'
        params['timestep'] = md_params.timestep if md_params.timestep else 0.5
        params['steps'] = md_params.n_steps if md_params.n_steps else 1000
        
        if md_params.thermostat == ThermostatType.NOSE_HOOVER:
            params['thermostat'] = 'NOSE'
            params['thermostat_params'] = {'temperature': md_params.temperature}
        elif md_params.thermostat == ThermostatType.LANGEVIN:
            params['thermostat'] = 'LANGEVIN'
            params['thermostat_params'] = {'temperature': md_params.temperature}
        
        return params
    
    def _get_qe_md_params(self, md_params: MDParameters) -> Dict[str, Any]:
        params = {}
        params['calculation'] = 'md'
        params['dt'] = md_params.timestep if md_params.timestep else 20.0
        params['nstep'] = md_params.n_steps if md_params.n_steps else 100
        
        if md_params.temperature:
            params['tempw'] = md_params.temperature
        
        if md_params.thermostat == ThermostatType.VELOCITY_SCALING:
            params['ion_dynamics'] = 'verlet'
        elif md_params.thermostat == ThermostatType.LANGEVIN:
            params['ion_dynamics'] = 'langevin'
        
        return params
    
    def _get_lammps_md_params(self, md_params: MDParameters) -> Dict[str, Any]:
        params = {}
        params['run_style'] = 'verlet'
        params['timestep'] = md_params.timestep if md_params.timestep else 1.0
        params['run'] = md_params.n_steps if md_params.n_steps else 10000
        
        if md_params.thermostat == ThermostatType.NOSE_HOOVER:
            params['fix_nvt'] = {
                'temperature': md_params.temperature,
                'tchain': 1
            }
        elif md_params.thermostat == ThermostatType.LANGEVIN:
            params['fix_langevin'] = {
                'temperature': md_params.temperature,
                'damping': 100.0
            }
        
        return params
    
    def _get_vasp_opt_params(self, opt_params: OptimizationParameters) -> Dict[str, Any]:
        params = {}
        params['ibrion'] = 2
        if opt_params.optimizer:
            if opt_params.optimizer.lower() == 'rmm-diis':
                params['ibrion'] = 1
            elif opt_params.optimizer.lower() == 'damped':
                params['ibrion'] = 3
        
        if opt_params.convergence_energy:
            params['ediff'] = opt_params.convergence_energy
        if opt_params.convergence_force:
            params['ediffg'] = -opt_params.convergence_force
        
        params['nsw'] = opt_params.max_steps if opt_params.max_steps else 200
        params['potim'] = opt_params.trust_radius if opt_params.trust_radius else 0.5
        
        return params
    
    def _get_cp2k_opt_params(self, opt_params: OptimizationParameters) -> Dict[str, Any]:
        params = {}
        params['run_type'] = 'GEO_OPT'
        params['optimizer'] = opt_params.optimizer.upper() if opt_params.optimizer else 'BFGS'
        params['max_iter'] = opt_params.max_steps if opt_params.max_steps else 200
        
        if opt_params.convergence_force:
            params['rms_force'] = opt_params.convergence_force
        
        return params
    
    def _get_qe_opt_params(self, opt_params: OptimizationParameters) -> Dict[str, Any]:
        params = {}
        params['calculation'] = 'relax'
        params['ion_dynamics'] = 'bfgs'
        params['nstep'] = opt_params.max_steps if opt_params.max_steps else 100
        
        if opt_params.convergence_energy:
            params['etot_conv_thr'] = opt_params.convergence_energy
        if opt_params.convergence_force:
            params['forc_conv_thr'] = opt_params.convergence_force
        if opt_params.convergence_stress:
            params['calculation'] = 'vc-relax'
            params['press_conv_thr'] = opt_params.convergence_stress
        
        return params
    
    def _get_gaussian_opt_params(self, opt_params: OptimizationParameters) -> Dict[str, Any]:
        params = {}
        route = ['#P', 'B3LYP/6-31G*']
        opt_type = 'Opt'
        if opt_params.optimizer:
            if opt_params.optimizer.lower() in ['calcfc', 'rfo', 'gdiis']:
                opt_type = f"Opt={opt_params.optimizer.upper()}"
        
        route.append(opt_type)
        params['route'] = ' '.join(route)
        
        if opt_params.convergence_force:
            if opt_params.convergence_force < 0.0001:
                params['opt'] = 'VeryTight'
            elif opt_params.convergence_force < 0.001:
                params['opt'] = 'Tight'
        
        params['max_steps'] = opt_params.max_steps if opt_params.max_steps else 100
        
        return params
    
    def _get_orca_opt_params(self, opt_params: OptimizationParameters) -> Dict[str, Any]:
        params = {}
        params['method'] = 'B3LYP'
        params['basis'] = 'def2-SVP'
        params['job_type'] = 'Opt'
        
        if opt_params.convergence_energy:
            params['energy_threshold'] = opt_params.convergence_energy
        if opt_params.convergence_force:
            params['grad_threshold'] = opt_params.convergence_force
        
        params['maxiter'] = opt_params.max_steps if opt_params.max_steps else 100
        
        return params
    
    def _get_code_specific_notes(self, target_code: str, workflow_params: WorkflowParameters) -> Optional[str]:
        notes = []
        
        if target_code == 'vasp':
            if workflow_params.workflow_type == WorkflowType.MOLECULAR_DYNAMICS:
                notes.append("VASP requires POTCAR file with appropriate pseudopotentials for MD.")
        elif target_code == 'gaussian':
            if workflow_params.workflow_type == WorkflowType.MOLECULAR_DYNAMICS:
                notes.append("Gaussian MD is limited to Born-Oppenheimer MD.")
        elif target_code == 'lammps':
            notes.append("LAMMPS requires separate force field parameter file.")
        
        if notes:
            return f"{target_code.upper()} Notes: " + " ".join(notes)
        return None
    
    def extract_workflow_from_file(
        self,
        filepath: str,
        code_hint: Optional[str] = None
    ) -> WorkflowMigrationResult:
        """Extract workflow parameters from an input file."""
        code = code_hint or self._detect_code_from_file(filepath)
        
        if not code:
            return WorkflowMigrationResult(
                success=False,
                warnings=["Could not detect code from file"]
            )
        
        try:
            if code == 'vasp':
                params = self._extract_vasp_workflow(filepath)
            elif code == 'cp2k':
                params = self._extract_cp2k_workflow(filepath)
            elif code == 'qe':
                params = self._extract_qe_workflow(filepath)
            elif code == 'gaussian':
                params = self._extract_gaussian_workflow(filepath)
            elif code == 'orca':
                params = self._extract_orca_workflow(filepath)
            else:
                return WorkflowMigrationResult(
                    success=False,
                    warnings=[f"Workflow extraction not implemented for {code}"]
                )
            
            # Convert dataclasses to dict with Enum values converted
            data = asdict(params)
            data = _convert_enum_to_value(data)
            
            return WorkflowMigrationResult(
                success=True,
                data=data,
                metadata={'source_code': code, 'source_file': filepath}
            )
        except Exception as e:
            return WorkflowMigrationResult(
                success=False,
                warnings=[f"Failed to extract workflow: {str(e)}"]
            )
    
    def _detect_code_from_file(self, filepath: str) -> Optional[str]:
        path = Path(filepath)
        filename = path.name.upper()
        extension = path.suffix.lower()
        
        if filename in ['INCAR', 'POSCAR', 'CONTCAR']:
            return 'vasp'
        
        if extension == '.inp' and path.exists():
            with open(filepath, 'r') as f:
                content = f.read()
                if '&FORCE_EVAL' in content or '&CP2K_INPUT' in content:
                    return 'cp2k'
                if '!' in content and any(kw in content for kw in ['Opt', 'Freq', 'SP']):
                    return 'orca'
        
        if extension in ['.in', '.pw.in']:
            return 'qe'
        
        if extension in ['.com', '.gjf']:
            return 'gaussian'
        
        return None
    
    def _extract_vasp_workflow(self, filepath: str) -> WorkflowParameters:
        incar_params = {}
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if '=' in line and not line.startswith('#'):
                    key, value = line.split('=', 1)
                    incar_params[key.strip().upper()] = value.strip()
        
        ibrion = int(incar_params.get('IBRION', 0))
        nsw = int(incar_params.get('NSW', 0))
        
        if ibrion == 0 and nsw > 0:
            md_params = MDParameters(
                timestep=float(incar_params.get('POTIM', 1.0)),
                n_steps=nsw,
                temperature=float(incar_params.get('TEBEG', 300.0))
            )
            return WorkflowParameters(workflow_type=WorkflowType.MOLECULAR_DYNAMICS, md_params=md_params)
        elif ibrion in [1, 2, 3] and nsw > 0:
            opt_params = OptimizationParameters(
                max_steps=nsw,
                convergence_force=float(incar_params.get('EDIFFG', -0.02)),
                convergence_energy=float(incar_params.get('EDIFF', 1e-4))
            )
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION, opt_params=opt_params)
        else:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
    
    def _extract_cp2k_workflow(self, filepath: str) -> WorkflowParameters:
        with open(filepath, 'r') as f:
            content = f.read()
        
        if 'MD' in content or '&MD' in content:
            return WorkflowParameters(workflow_type=WorkflowType.MOLECULAR_DYNAMICS)
        elif 'GEO_OPT' in content or 'OPT' in content:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
        else:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
    
    def _extract_qe_workflow(self, filepath: str) -> WorkflowParameters:
        with open(filepath, 'r') as f:
            content = f.read()
        
        calculation = 'scf'
        for line in content.split('\n'):
            if 'calculation' in line.lower():
                calculation = line.split('=')[1].strip().strip("'\"")
                break
        
        if calculation == 'md':
            return WorkflowParameters(workflow_type=WorkflowType.MOLECULAR_DYNAMICS)
        elif calculation in ['relax', 'vc-relax']:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
        else:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
    
    def _extract_gaussian_workflow(self, filepath: str) -> WorkflowParameters:
        with open(filepath, 'r') as f:
            content = f.read()
        
        route_line = ''
        for line in content.split('\n'):
            if line.startswith('#'):
                route_line = line
                break
        
        if 'opt' in route_line.lower():
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
        elif 'freq' in route_line.lower():
            return WorkflowParameters(workflow_type=WorkflowType.FREQUENCY_CALCULATION)
        else:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
    
    def _extract_orca_workflow(self, filepath: str) -> WorkflowParameters:
        with open(filepath, 'r') as f:
            content = f.read()
        
        if 'opt' in content.lower():
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)
        elif 'freq' in content.lower():
            return WorkflowParameters(workflow_type=WorkflowType.FREQUENCY_CALCULATION)
        else:
            return WorkflowParameters(workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION)


def main():
    import sys
    
    if len(sys.argv) < 4:
        print("Usage: python workflow.py <source_code> <target_code> <workflow.json>")
        print("   or: python workflow.py extract <filepath> [code_hint]")
        sys.exit(1)
    
    migrator = WorkflowMigrator()
    
    if sys.argv[1] == 'extract':
        filepath = sys.argv[2]
        code_hint = sys.argv[3] if len(sys.argv) > 3 else None
        result = migrator.extract_workflow_from_file(filepath, code_hint)
    else:
        source_code = sys.argv[1]
        target_code = sys.argv[2]
        workflow_file = sys.argv[3]
        
        with open(workflow_file, 'r') as f:
            workflow_dict = json.load(f)
        
        workflow_type = WorkflowType(workflow_dict['workflow_type'])
        md_params = None
        opt_params = None
        
        if 'md_params' in workflow_dict and workflow_dict['md_params']:
            md_dict = workflow_dict['md_params']
            md_params = MDParameters(
                temperature=md_dict.get('temperature'),
                timestep=md_dict.get('timestep'),
                n_steps=md_dict.get('n_steps')
            )
        
        if 'opt_params' in workflow_dict and workflow_dict['opt_params']:
            opt_dict = workflow_dict['opt_params']
            opt_params = OptimizationParameters(
                max_steps=opt_dict.get('max_steps'),
                convergence_force=opt_dict.get('convergence_force'),
                convergence_energy=opt_dict.get('convergence_energy')
            )
        
        workflow_params = WorkflowParameters(
            workflow_type=workflow_type,
            md_params=md_params,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow(source_code, target_code, workflow_params)
    
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == '__main__':
    main()
