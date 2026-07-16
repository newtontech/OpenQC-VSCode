import sys
import tempfile
import types
import unittest
from pathlib import Path

from openqc.bridge.output_bridge import _parse_output_native, main


class OutputBridgeTest(unittest.TestCase):
    def test_native_output_parser_extracts_scf_energy_series(self):
        results = _parse_output_native(
            "\n".join(
                [
                    " SCF Done:  E(RHF) =  -75.500000 A.U. after 1 cycles",
                    " SCF Done:  E(RHF) =  -76.100000 A.U. after 2 cycles",
                ]
            ),
            "gaussian",
        )

        self.assertTrue(results["success"])
        self.assertEqual(results["software"], "gaussian")
        self.assertEqual(results["finalEnergy"], {"value": -76.1, "unit": "hartree"})
        self.assertEqual(results["scfEnergies"], [-75.5, -76.1])

    def test_native_output_parser_extracts_orca_final_energy_without_cclib(self):
        results = _parse_output_native(
            "\n".join(
                [
                    "* O   R   C   A *",
                    "FINAL SINGLE POINT ENERGY     -75.983271231",
                    "ORCA TERMINATED NORMALLY",
                ]
            ),
            "auto",
        )

        self.assertTrue(results["success"])
        self.assertEqual(results["software"], "orca")
        self.assertEqual(results["finalEnergy"], {"value": -75.983271231, "unit": "hartree"})
        self.assertEqual(results["scfEnergies"], [-75.983271231])

    def test_native_output_parser_extracts_cp2k_energy_without_cclib(self):
        results = _parse_output_native(
            "\n".join(
                [
                    "CP2K| version string",
                    "ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):      -76.234567890",
                    "ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):      -76.345678901",
                ]
            ),
            "auto",
        )

        self.assertTrue(results["success"])
        self.assertEqual(results["software"], "cp2k")
        self.assertEqual(results["finalEnergy"], {"value": -76.345678901, "unit": "hartree"})
        self.assertEqual(results["scfEnergies"], [-76.23456789, -76.345678901])

    def test_native_output_parser_extracts_qe_energy_without_cclib(self):
        results = _parse_output_native(
            "\n".join(
                [
                    "Program PWSCF v.7.2 starts",
                    "!    total energy              =     -15.25000000 Ry",
                ]
            ),
            "auto",
        )

        self.assertTrue(results["success"])
        self.assertEqual(results["software"], "qe")
        self.assertEqual(results["finalEnergy"]["unit"], "eV")
        self.assertAlmostEqual(results["finalEnergy"]["value"], -207.4868201256585)
        self.assertAlmostEqual(results["scfEnergies"][0], -207.4868201256585)

    def test_native_output_parser_extracts_gamess_energy_without_cclib(self):
        results = _parse_output_native(
            "\n".join(
                [
                    "GAMESS VERSION = 00",
                    "          TOTAL ENERGY =      -75.654321000",
                ]
            ),
            "auto",
        )

        self.assertTrue(results["success"])
        self.assertEqual(results["software"], "gamess")
        self.assertEqual(results["finalEnergy"], {"value": -75.654321, "unit": "hartree"})
        self.assertEqual(results["scfEnergies"], [-75.654321])

    def test_native_output_parser_reports_unsupported_for_junk_text(self):
        results = _parse_output_native(
            "this is not a quantum chemistry output file\njust ordinary notes\n",
            "auto",
        )

        self.assertFalse(results["success"])
        self.assertEqual(results["software"], "auto")
        self.assertNotIn("finalEnergy", results)
        self.assertNotIn("scfEnergies", results)
        self.assertIn("No calculation output data", results["warnings"][0])

    def test_parse_command_falls_back_when_cclib_cannot_parse_file(self):
        fake_cclib = types.SimpleNamespace(
            io=types.SimpleNamespace(
                ccread=lambda _path: (_ for _ in ()).throw(ValueError("not a cclib file"))
            )
        )
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "gaussian.log"
                output_path.write_text(
                    " SCF Done:  E(RHF) =  -76.000000 A.U. after 3 cycles\n",
                    encoding="utf-8",
                )

                results = _run_output_bridge_parse(output_path)
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertTrue(results["success"])
        self.assertFalse(results["data"]["cclibAvailable"])
        self.assertEqual(results["data"]["finalEnergy"], {"value": -76.0, "unit": "hartree"})
        self.assertIn("cclib parser unavailable or failed", results["data"]["warnings"][0])

    def test_parse_command_falls_back_when_cclib_returns_none(self):
        fake_cclib = types.SimpleNamespace(io=types.SimpleNamespace(ccread=lambda _path: None))
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "gaussian.log"
                output_path.write_text(
                    " SCF Done:  E(RHF) =  -77.000000 A.U. after 3 cycles\n",
                    encoding="utf-8",
                )

                results = _run_output_bridge_parse(output_path)
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertFalse(results["data"]["cclibAvailable"])
        self.assertEqual(results["data"]["finalEnergy"], {"value": -77.0, "unit": "hartree"})

    def test_parse_command_does_not_claim_success_for_unrecognized_output(self):
        fake_cclib = types.SimpleNamespace(io=types.SimpleNamespace(ccread=lambda _path: None))
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "notes.out"
                output_path.write_text(
                    "ordinary lab notes\nnot a calculation output\n",
                    encoding="utf-8",
                )

                results = _run_output_bridge_parse(output_path)
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertTrue(results["success"])
        self.assertFalse(results["data"]["success"])
        self.assertFalse(results["data"]["cclibAvailable"])
        self.assertIn("No calculation output data", results["data"]["warnings"][0])
        self.assertIn("cclib parser unavailable or failed", results["data"]["warnings"][1])

    def test_summarize_command_uses_native_fallback_when_cclib_fails(self):
        fake_cclib = types.SimpleNamespace(io=types.SimpleNamespace(ccread=lambda _path: None))
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "orca.out"
                output_path.write_text(
                    "* O   R   C   A *\nFINAL SINGLE POINT ENERGY     -75.983271231\n",
                    encoding="utf-8",
                )

                results = _run_output_bridge_command("summarize", output_path, software="auto")
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertTrue(results["success"])
        self.assertEqual(results["data"]["software"], "orca")
        self.assertTrue(results["data"]["success"])
        self.assertEqual(results["data"]["finalEnergy"], {"value": -75.983271231, "unit": "hartree"})
        self.assertEqual(results["data"]["scfStepCount"], 1)
        self.assertFalse(results["data"]["cclibAvailable"])

    def test_trajectory_command_extracts_cclib_atomcoord_frames(self):
        fake_data = types.SimpleNamespace(
            metadata={"package": "gaussian"},
            atomnos=[1, 8],
            atomcoords=[
                [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                [[0.0, 0.0, 0.1], [0.0, 0.0, 1.1]],
            ],
            scfenergies=[-75.0, -76.0],
        )
        fake_cclib = types.SimpleNamespace(io=types.SimpleNamespace(ccread=lambda _path: fake_data))
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "gaussian.log"
                output_path.write_text("fake gaussian optimization\n", encoding="utf-8")

                results = _run_output_bridge_command("trajectory", output_path)
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertTrue(results["success"])
        data = results["data"]
        self.assertTrue(data["supported"])
        self.assertEqual(data["sourceFile"], str(output_path))
        self.assertEqual(data["software"], "gaussian")
        self.assertEqual(data["frameCount"], 2)
        self.assertEqual(data["energies"], [-75.0, -76.0])
        self.assertEqual(data["frames"][0]["schemaVersion"], "openqc.structure.v1")
        self.assertEqual(data["frames"][0]["atoms"][1]["element"], "O")
        self.assertEqual(data["frames"][1]["atoms"][1]["z"], 1.1)

    def test_trajectory_command_returns_structured_unsupported_without_frames(self):
        fake_data = types.SimpleNamespace(metadata={"package": "orca"}, scfenergies=[-10.0])
        fake_cclib = types.SimpleNamespace(io=types.SimpleNamespace(ccread=lambda _path: fake_data))
        previous = sys.modules.get("cclib")
        sys.modules["cclib"] = fake_cclib
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_path = Path(tmpdir) / "orca.out"
                output_path.write_text("fake single point\n", encoding="utf-8")

                results = _run_output_bridge_command("trajectory", output_path)
        finally:
            if previous is None:
                sys.modules.pop("cclib", None)
            else:
                sys.modules["cclib"] = previous

        self.assertTrue(results["success"])
        self.assertFalse(results["data"]["supported"])
        self.assertEqual(results["data"]["frames"], [])
        self.assertEqual(results["data"]["frameCount"], 0)
        self.assertIn("No atom coordinate trajectory", results["data"]["warnings"][0])


def _run_output_bridge_parse(output_path: Path):
    return _run_output_bridge_command("parse", output_path)


def _run_output_bridge_command(command: str, output_path: Path, software: str = "gaussian"):
    import io
    import json

    previous_stdin = sys.stdin
    previous_stdout = sys.stdout
    request = {
        "command": command,
        "args": {"path": str(output_path), "software": software},
    }
    captured_stdout = io.StringIO()
    try:
        sys.stdin = io.StringIO(json.dumps(request))
        sys.stdout = captured_stdout
        main()
    finally:
        sys.stdin = previous_stdin
        sys.stdout = previous_stdout

    return json.loads(captured_stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
