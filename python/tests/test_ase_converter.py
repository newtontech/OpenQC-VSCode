"""
Test ASE Converter Backend

Unit tests for the Python ASE converter module.
"""

import pytest
import tempfile
import os
from pathlib import Path
from ase import Atoms
from ase.build import molecule, bulk
import numpy as np

# Import after creating the module
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

try:
    from openqc.ase.converter import ASEConverter, ConversionResult
    from openqc.ase.utils import validate_atoms, get_supported_formats
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


@pytest.mark.skipif(not ASE_AVAILABLE, reason="ASE module not available")
class TestASEConverter:
    """Test ASE Converter class."""

    def setup_method(self):
        """Setup test fixtures."""
        self.converter = ASEConverter()
        self.temp_dir = tempfile.mkdtemp()

    def teardown_method(self):
        """Cleanup test fixtures."""
        import shutil
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_instantiation(self):
        """Test converter instantiation."""
        assert self.converter is not None
        assert hasattr(self.converter, 'supported_formats')

    def test_supported_formats(self):
        """Test that all expected formats are supported."""
        formats = self.converter.supported_formats
        
        assert 'vasp' in formats
        assert 'cp2k' in formats
        assert 'qe' in formats
        assert 'gaussian' in formats

    def test_read_xyz_to_atoms(self):
        """Test reading XYZ file to ASE Atoms."""
        xyz_content = """3
Water molecule
O  0.000000  0.000000  0.000000
H  0.758602  0.000000  0.504284
H -0.758602  0.000000  0.504284
"""
        xyz_path = os.path.join(self.temp_dir, 'test.xyz')
        with open(xyz_path, 'w') as f:
            f.write(xyz_content)

        result = self.converter.read_to_atoms(xyz_path)

        assert result.success
        assert isinstance(result.data, Atoms)
        assert len(result.data) == 3

    def test_read_nonexistent_file(self):
        """Test reading a non-existent file."""
        result = self.converter.read_to_atoms('/nonexistent/file.xyz')
        
        assert not result.success
        assert 'File not found' in result.error


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
