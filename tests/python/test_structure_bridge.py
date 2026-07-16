import unittest
import tempfile
from pathlib import Path

from openqc.bridge.structure_bridge import _generate_supercell, convert_structure, parse_file


class StructureBridgeSupercellTest(unittest.TestCase):
    def test_fractional_supercell_keeps_fractional_coordinates(self):
        structure = {
            "name": "Si",
            "atoms": [
                {"element": "Si", "x": 0.0, "y": 0.0, "z": 0.0},
                {"element": "Si", "x": 0.25, "y": 0.25, "z": 0.25},
            ],
            "cell": {
                "a": [2.715, 2.715, 0.0],
                "b": [0.0, 2.715, 2.715],
                "c": [2.715, 0.0, 2.715],
                "pbc": [True, True, True],
                "coordinateMode": "fractional",
            },
            "metadata": {"source": {"format": "poscar"}},
        }

        supercell = _generate_supercell(structure, 2, 1, 1)

        self.assertEqual(supercell["cell"]["coordinateMode"], "fractional")
        self.assertEqual(supercell["cell"]["a"], [5.43, 5.43, 0.0])
        self.assertEqual(len(supercell["atoms"]), 4)
        self.assertEqual(supercell["atoms"][0]["x"], 0.0)
        self.assertEqual(supercell["atoms"][1]["x"], 0.125)
        self.assertEqual(supercell["atoms"][2]["x"], 0.5)
        self.assertEqual(supercell["atoms"][3]["x"], 0.625)

    def test_file_only_xyz_uses_native_parser_without_optional_dependencies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "water.xyz"
            xyz_path.write_text(
                "3\nwater\nO 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n",
                encoding="utf-8",
            )

            structure = parse_file(str(xyz_path))

        self.assertEqual(structure["kind"], "molecule")
        self.assertEqual(len(structure["atoms"]), 3)
        self.assertEqual(structure["atoms"][0]["element"], "O")

    def test_file_only_poscar_uses_native_parser_without_optional_dependencies(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            poscar_path = Path(tmpdir) / "POSCAR"
            poscar_path.write_text(
                "\n".join(
                    [
                        "Si",
                        "1.0",
                        "5.43 0 0",
                        "0 5.43 0",
                        "0 0 5.43",
                        "Si",
                        "2",
                        "Direct",
                        "0 0 0",
                        "0.25 0.25 0.25",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            structure = parse_file(str(poscar_path))

        self.assertEqual(structure["kind"], "periodic")
        self.assertEqual(structure["cell"]["coordinateMode"], "fractional")
        self.assertEqual(len(structure["atoms"]), 2)

    def test_cp2k_inp_autodetects_and_uses_native_coord_parser(self):
        fixture = Path(__file__).parents[1] / "fixtures" / "cp2k" / "H2O.inp"

        structure = parse_file(str(fixture))

        self.assertEqual(structure["metadata"]["source"]["format"], "cp2k")
        self.assertEqual(structure["metadata"]["source"]["software"], "cp2k")
        self.assertEqual(structure["kind"], "periodic")
        self.assertEqual(len(structure["atoms"]), 3)
        self.assertEqual(structure["atoms"][0]["element"], "O")
        self.assertEqual(structure["cell"]["a"], [10.0, 0.0, 0.0])

    def test_gamess_inp_autodetects_and_uses_native_data_parser(self):
        fixture = Path(__file__).parents[1] / "fixtures" / "gamess" / "H2O.inp"

        structure = parse_file(str(fixture))

        self.assertEqual(structure["metadata"]["source"]["format"], "gamess")
        self.assertEqual(structure["metadata"]["source"]["software"], "gamess")
        self.assertEqual(structure["kind"], "molecule")
        self.assertEqual(len(structure["atoms"]), 3)
        self.assertEqual(structure["atoms"][0]["element"], "O")

    def test_convert_xyz_to_pdb_returns_target_content_not_openqc_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xyz_path = Path(tmpdir) / "water.xyz"
            xyz_path.write_text(
                "3\nwater\nO 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n",
                encoding="utf-8",
            )

            structure = parse_file(str(xyz_path))
            converted = convert_structure(structure, "pdb")

        self.assertEqual(converted["format"], "pdb")
        self.assertEqual(converted["atomCount"], 3)
        self.assertIn("HETATM", converted["content"])
        self.assertIn("END", converted["content"])
        self.assertNotIn("schemaVersion", converted)

    def test_convert_fractional_poscar_to_xyz_uses_cartesian_coordinates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            poscar_path = Path(tmpdir) / "POSCAR"
            poscar_path.write_text(
                "\n".join(
                    [
                        "Si",
                        "1.0",
                        "5.43 0 0",
                        "0 5.43 0",
                        "0 0 5.43",
                        "Si",
                        "2",
                        "Direct",
                        "0 0 0",
                        "0.25 0.25 0.25",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            structure = parse_file(str(poscar_path))
            converted = convert_structure(structure, "xyz")

        self.assertEqual(converted["format"], "xyz")
        self.assertIn("Si 1.357500 1.357500 1.357500", converted["content"])


if __name__ == "__main__":
    unittest.main()
