import unittest

from openqc.bridge.check_backend import _derive_capabilities


class CheckBackendCapabilityTest(unittest.TestCase):
    def test_structure_parsing_reports_native_degraded_without_optional_parsers(self):
        packages = {
            "ase": {"available": False},
            "pymatgen": {"available": False},
            "cclib": {"available": False},
            "dpdata": {"available": False},
            "spglib": {"available": False},
        }
        external_tools = {
            "multiwfn": {"available": False},
            "c2x": {"available": False},
            "openbabel": {"available": False},
        }

        capabilities = _derive_capabilities(packages, external_tools)
        structure = capabilities["structureParsing"]

        self.assertEqual(structure["status"], "degraded")
        self.assertIn("ase", structure["requires"])
        self.assertIn("pymatgen", structure["requires"])
        self.assertIn("Native", structure["detail"])

    def test_external_analysis_reports_openbabel_requirement(self):
        packages = {
            "ase": {"available": False},
            "pymatgen": {"available": False},
            "cclib": {"available": False},
            "dpdata": {"available": False},
            "spglib": {"available": False},
        }
        external_tools = {
            "multiwfn": {"available": False},
            "c2x": {"available": False},
            "openbabel": {"available": False},
        }

        capabilities = _derive_capabilities(packages, external_tools)
        external = capabilities["externalAnalysis"]

        self.assertEqual(external["status"], "missing")
        self.assertIn("Open Babel", external["requires"])
        self.assertIn("Open Babel", external["detail"])


if __name__ == "__main__":
    unittest.main()
