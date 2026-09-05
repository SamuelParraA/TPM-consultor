import csv
import unittest
from pathlib import Path

from app import attributes, strip_prefix

ROOT = Path(__file__).resolve().parents[1]


class ParserTests(unittest.TestCase):
    def test_gff_attributes(self):
        parsed = attributes("ID=gene:Solyc01g000010.4;Alias=ACS7;Note=ethylene enzyme")
        self.assertEqual(parsed["Alias"], "ACS7")
        self.assertEqual(strip_prefix(parsed["ID"]), "Solyc01g000010.4")


class PackagingTests(unittest.TestCase):
    def test_dockerfile_copy_inputs_exist(self):
        dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
        for source in ("app.py", "config.example.json", "sample_metadata.csv", "web"):
            self.assertIn(source, dockerfile)
            self.assertTrue((ROOT / source).exists(), f"Dockerfile requiere {source}")

    def test_default_sample_metadata(self):
        with (ROOT / "sample_metadata.csv").open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 18)
        self.assertEqual(set(rows[0]), {"sample", "timepoint", "replicate", "order"})
        self.assertEqual(len({row["sample"] for row in rows}), 18)

if __name__ == "__main__":
    unittest.main()
