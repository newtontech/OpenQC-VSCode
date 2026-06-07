"""
Test Workflow Migration Module

Unit tests for the workflow migration module for ASE integration.
"""

import pytest
import tempfile
import os
from pathlib import Path

from core.openqc.workflow import (
    WorkflowMigrator,
    WorkflowParameters,
    MDParameters,
    OptimizationParameters,
    WorkflowType,
    ThermostatType
)


@pytest.fixture
def migrator():
    """Create workflow migrator instance."""
    return WorkflowMigrator()


@pytest.fixture
def temp_dir():
    """Create temporary directory for test files."""
    td = tempfile.mkdtemp()
    yield td
    import shutil
    shutil.rmtree(td, ignore_errors=True)


class TestMDParameters:
    """Test MD parameters dataclass."""
    
    def test_create_md_params(self):
        """Test creating MD parameters."""
        params = MDParameters(
            temperature=300.0,
            pressure=1.0,
            timestep=1.0,
            n_steps=1000,
            thermostat=ThermostatType.NOSE_HOOVER
        )
        assert params.temperature == 300.0
        assert params.pressure == 1.0
        assert params.timestep == 1.0
        assert params.n_steps == 1000
        assert params.thermostat == ThermostatType.NOSE_HOOVER


class TestOptimizationParameters:
    """Test optimization parameters dataclass."""
    
    def test_create_opt_params(self):
        """Test creating optimization parameters."""
        params = OptimizationParameters(
            convergence_energy=1e-4,
            convergence_force=0.02,
            max_steps=100,
            optimizer='bfgs'
        )
        assert params.convergence_energy == 1e-4
        assert params.convergence_force == 0.02
        assert params.max_steps == 100
        assert params.optimizer == 'bfgs'


class TestWorkflowMigrator:
    """Test workflow migrator class."""
    
    def test_initialization(self, migrator):
        """Test migrator initialization."""
        assert migrator is not None
        assert len(migrator.supported_codes) == 6
        assert 'vasp' in migrator.supported_codes
        assert 'cp2k' in migrator.supported_codes
        assert 'qe' in migrator.supported_codes
        assert 'gaussian' in migrator.supported_codes
        assert 'orca' in migrator.supported_codes
        assert 'lammps' in migrator.supported_codes
    
    def test_unsupported_source_code(self, migrator):
        """Test migration with unsupported source code."""
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION
        )
        result = migrator.migrate_workflow(
            'unsupported_code',
            'vasp',
            workflow_params
        )
        assert not result.success
        assert len(result.warnings) > 0
        assert 'Unsupported source code' in result.warnings[0]
    
    def test_unsupported_target_code(self, migrator):
        """Test migration with unsupported target code."""
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION
        )
        result = migrator.migrate_workflow(
            'vasp',
            'unsupported_code',
            workflow_params
        )
        assert not result.success
        assert len(result.warnings) > 0
        assert 'Unsupported target code' in result.warnings[0]
    
    def test_migrate_md_vasp_to_cp2k(self, migrator):
        """Test MD migration from VASP to CP2K."""
        md_params = MDParameters(
            temperature=300.0,
            timestep=1.0,
            n_steps=1000,
            thermostat=ThermostatType.NOSE_HOOVER
        )
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        result = migrator.migrate_workflow('vasp', 'cp2k', workflow_params)
        
        assert result.success
        assert 'temperature' in result.data
        assert result.data['temperature'] == 300.0
        assert 'timestep' in result.data
        assert result.data['timestep'] == 1.0
        assert 'n_steps' in result.data
        assert result.data['n_steps'] == 1000
        assert 'thermostat' in result.data
        assert result.data['thermostat'] == 'NOSE'
    
    def test_migrate_md_vasp_to_qe(self, migrator):
        """Test MD migration from VASP to QE."""
        md_params = MDParameters(
            temperature=300.0,
            timestep=1.0,
            n_steps=100,
            thermostat=ThermostatType.VELOCITY_SCALING
        )
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        result = migrator.migrate_workflow('vasp', 'qe', workflow_params)
        
        assert result.success
        assert 'temperature' in result.data
        assert 'timestep' in result.data
        assert 'n_steps' in result.data
        assert 'ion_dynamics' in result.data
        assert result.data['ion_dynamics'] == 'verlet'
    
    def test_migrate_opt_vasp_to_cp2k(self, migrator):
        """Test optimization migration from VASP to CP2K."""
        opt_params = OptimizationParameters(
            convergence_energy=1e-4,
            convergence_force=0.02,
            max_steps=100,
            optimizer='bfgs'
        )
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow('vasp', 'cp2k', workflow_params)
        
        assert result.success
        assert 'convergence_energy' in result.data
        assert 'convergence_force' in result.data
        assert 'max_steps' in result.data
        assert 'optimizer' in result.data
        assert result.data['optimizer'] == 'BFGS'
        assert 'run_type' in result.data
        assert result.data['run_type'] == 'GEO_OPT'
    
    def test_migrate_opt_vasp_to_gaussian(self, migrator):
        """Test optimization migration from VASP to Gaussian."""
        opt_params = OptimizationParameters(
            convergence_force=0.0003,
            max_steps=50
        )
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow('vasp', 'gaussian', workflow_params)
        
        assert result.success
        assert 'convergence_force' in result.data
        assert 'max_steps' in result.data
        assert 'route' in result.data
        assert result.data['opt'] == 'Tight'
    
    def test_migrate_default_convergence(self, migrator):
        """Test migration uses default convergence when not provided."""
        opt_params = OptimizationParameters(
            max_steps=100
        )
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow('vasp', 'cp2k', workflow_params)
        
        assert result.success
        assert 'convergence_energy' in result.data
        assert 'convergence_force' in result.data
        # Should use CP2K defaults
        assert result.data['convergence_energy'] == 1e-3
        assert result.data['convergence_force'] == 0.003
    
    def test_code_specific_warnings(self, migrator):
        """Test code-specific warnings are generated."""
        # Gaussian does not support NPT ensemble
        md_params = MDParameters(pressure=1.0)
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        result = migrator.migrate_workflow('vasp', 'gaussian', workflow_params)
        
        assert result.success
        assert len(result.warnings) > 0
        assert any('NPT' in w for w in result.warnings)
    
    def test_pressure_not_supported_warning(self, migrator):
        """Test pressure warning for Gaussian/ORCA."""
        md_params = MDParameters(pressure=1.0, temperature=300.0)
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        for target in ['gaussian', 'orca']:
            result = migrator.migrate_workflow('vasp', target, workflow_params)
            assert any('NPT' in w for w in result.warnings)
    
    def test_stress_convergence_warning(self, migrator):
        """Test stress convergence warning for Gaussian/ORCA."""
        opt_params = OptimizationParameters(convergence_stress=0.5)
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        for target in ['gaussian', 'orca']:
            result = migrator.migrate_workflow('vasp', target, workflow_params)
            assert any('stress' in w for w in result.warnings)


class TestWorkflowExtraction:
    """Test workflow extraction from files."""
    
    def test_detect_vascar_code(self, migrator, temp_dir):
        """Test detecting VASP INCAR file."""
        incar_path = os.path.join(temp_dir, 'INCAR')
        with open(incar_path, 'w') as f:
            f.write('IBRION = 2\nNSW = 100\n')
        
        code = migrator._detect_code_from_file(incar_path)
        assert code == 'vasp'
    
    def test_detect_gaussian_code(self, migrator, temp_dir):
        """Test detecting Gaussian input file."""
        com_path = os.path.join(temp_dir, 'test.com')
        with open(com_path, 'w') as f:
            f.write('#P B3LYP/6-31G* Opt\n\n0 1\nO\n')
        
        code = migrator._detect_code_from_file(com_path)
        assert code == 'gaussian'
    
    def test_detect_qe_code(self, migrator, temp_dir):
        """Test detecting Quantum ESPRESSO input file."""
        qe_path = os.path.join(temp_dir, 'pw.in')
        with open(qe_path, 'w') as f:
            f.write('&control\n calculation = \'relax\'\n/\n')
        
        code = migrator._detect_code_from_file(qe_path)
        assert code == 'qe'
    
    def test_extract_vasp_md(self, migrator, temp_dir):
        """Test extracting VASP MD workflow."""
        incar_path = os.path.join(temp_dir, 'INCAR')
        with open(incar_path, 'w') as f:
            f.write('''IBRION = 0
NSW = 1000
POTIM = 1.0
TEBEG = 300.0
TEEND = 300.0
''')
        
        result = migrator.extract_workflow_from_file(incar_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'molecular_dynamics'
        assert result.metadata['source_code'] == 'vasp'
        assert 'md_params' in result.data
        assert result.data['md_params']['temperature'] == 300.0
        assert result.data['md_params']['n_steps'] == 1000
    
    def test_extract_vasp_optimization(self, migrator, temp_dir):
        """Test extracting VASP optimization workflow."""
        incar_path = os.path.join(temp_dir, 'INCAR')
        with open(incar_path, 'w') as f:
            f.write('''IBRION = 2
NSW = 100
EDIFF = 1e-4
EDIFFG = -0.02
''')
        
        result = migrator.extract_workflow_from_file(incar_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'geometry_optimization'
        assert 'opt_params' in result.data
        assert result.data['opt_params']['convergence_energy'] == 1e-4
        assert result.data['opt_params']['convergence_force'] == -0.02
    
    def test_extract_gaussian_opt(self, migrator, temp_dir):
        """Test extracting Gaussian optimization workflow."""
        com_path = os.path.join(temp_dir, 'test.com')
        with open(com_path, 'w') as f:
            f.write('#P B3LYP/6-31G* Opt\n\n0 1\nO 0 0 0\n')
        
        result = migrator.extract_workflow_from_file(com_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'geometry_optimization'
        assert result.metadata['source_code'] == 'gaussian'
    
    def test_extract_gaussian_freq(self, migrator, temp_dir):
        """Test extracting Gaussian frequency calculation."""
        com_path = os.path.join(temp_dir, 'test.com')
        with open(com_path, 'w') as f:
            f.write('#P B3LYP/6-31G* Freq\n\n0 1\nO 0 0 0\n')
        
        result = migrator.extract_workflow_from_file(com_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'frequency_calculation'
    
    def test_extract_qe_relax(self, migrator, temp_dir):
        """Test extracting QE relaxation workflow."""
        qe_path = os.path.join(temp_dir, 'pw.in')
        with open(qe_path, 'w') as f:
            f.write('''&control
 calculation = 'relax'
 nstep = 100
/
&ions
 ion_dynamics = 'bfgs'
/
''')
        
        result = migrator.extract_workflow_from_file(qe_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'geometry_optimization'
    
    def test_extract_qe_md(self, migrator, temp_dir):
        """Test extracting QE MD workflow."""
        qe_path = os.path.join(temp_dir, 'pw.in')
        with open(qe_path, 'w') as f:
            f.write('''&control
 calculation = 'md'
 dt = 20.0
 nstep = 1000
 tempw = 300.0
/
''')
        
        result = migrator.extract_workflow_from_file(qe_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'molecular_dynamics'
    
    def test_extract_cp2k_md(self, migrator, temp_dir):
        """Test extracting CP2K MD workflow."""
        inp_path = os.path.join(temp_dir, 'test.inp')
        with open(inp_path, 'w') as f:
            f.write('''&FORCE_EVAL
 &MOTION
  &MD
   STEPS 1000
   TIMESTEP 0.5
  &END MD
 &END MOTION
&END FORCE_EVAL
''')
        
        result = migrator.extract_workflow_from_file(inp_path)
        
        assert result.success
        assert result.data['workflow_type'] == 'molecular_dynamics'
        assert result.metadata['source_code'] == 'cp2k'


class TestOptimizerMapping:
    """Test optimizer parameter mappings."""
    
    def test_vasp_optimizer_mapping(self, migrator):
        """Test VASP optimizer mapping."""
        opt_params = OptimizationParameters(optimizer='cg')
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow('cp2k', 'vasp', workflow_params)
        
        assert result.success
        assert result.data['optimizer'] == 'CG'
        assert result.data['ibrion'] == 2
    
    def test_gaussian_optimizer_mapping(self, migrator):
        """Test Gaussian optimizer mapping."""
        opt_params = OptimizationParameters(optimizer='rfo')
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.GEOMETRY_OPTIMIZATION,
            opt_params=opt_params
        )
        
        result = migrator.migrate_workflow('vasp', 'gaussian', workflow_params)
        
        assert result.success
        assert 'Opt=RFO' in result.data['route']


class TestThermostatMapping:
    """Test thermostat parameter mappings."""
    
    def test_nose_hoover_mapping(self, migrator):
        """Test Nose-Hoover thermostat mapping."""
        md_params = MDParameters(thermostat=ThermostatType.NOSE_HOOVER)
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        # VASP
        result_vasp = migrator.migrate_workflow('cp2k', 'vasp', workflow_params)
        assert result_vasp.success
        assert result_vasp.data['thermostat'] == 'SMASS'
        assert result_vasp.data['smaas'] == 0
        
        # CP2K
        result_cp2k = migrator.migrate_workflow('vasp', 'cp2k', workflow_params)
        assert result_cp2k.success
        assert result_cp2k.data['thermostat'] == 'NOSE'
    
    def test_langevin_mapping(self, migrator):
        """Test Langevin thermostat mapping."""
        md_params = MDParameters(thermostat=ThermostatType.LANGEVIN)
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS,
            md_params=md_params
        )
        
        # QE
        result_qe = migrator.migrate_workflow('vasp', 'qe', workflow_params)
        assert result_qe.success
        assert result_qe.data['ion_dynamics'] == 'langevin'
        
        # LAMMPS
        result_lammps = migrator.migrate_workflow('vasp', 'lammps', workflow_params)
        assert result_lammps.success
        assert result_lammps.data['fix_langevin']


class TestCodeSpecificNotes:
    """Test code-specific notes generation."""
    
    def test_vasp_md_notes(self, migrator):
        """Test VASP MD notes."""
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS
        )
        result = migrator.migrate_workflow('cp2k', 'vasp', workflow_params)
        
        assert any('POTCAR' in w for w in result.warnings)
    
    def test_gaussian_md_notes(self, migrator):
        """Test Gaussian MD notes."""
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS
        )
        result = migrator.migrate_workflow('vasp', 'gaussian', workflow_params)
        
        assert any('Born-Oppenheimer' in w for w in result.warnings)
    
    def test_lammps_md_notes(self, migrator):
        """Test LAMMPS MD notes."""
        workflow_params = WorkflowParameters(
            workflow_type=WorkflowType.MOLECULAR_DYNAMICS
        )
        result = migrator.migrate_workflow('vasp', 'lammps', workflow_params)
        
        assert any('force field' in w for w in result.warnings)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
