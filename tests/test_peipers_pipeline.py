import csv
import subprocess
import sys
import unittest
from collections import Counter
from pathlib import Path

from app import load_metadata

ROOT = Path(__file__).resolve().parents[1]
PIPELINE = ROOT / "peipers_rnaseq"
RUNS = PIPELINE / "manifests" / "PRJNA1347747_runs.tsv"
SAMPLES = PIPELINE / "manifests" / "PRJNA1347747_samples.tsv"


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PeipersPipelineTests(unittest.TestCase):
    def test_project_validator(self):
        command = [
            sys.executable, str(PIPELINE / "scripts" / "validate_project.py"),
            "--runs", str(RUNS), "--samples", str(SAMPLES),
            "--expected-project", "PRJNA1347747",
        ]
        completed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        self.assertIn("18 corridas", completed.stdout)

    def test_design_is_balanced(self):
        rows = read_tsv(RUNS)
        counts = Counter((row["genotype"], row["treatment"]) for row in rows)
        expected = {
            (genotype, treatment): 3
            for genotype in ("WT", "epi", "ACCD")
            for treatment in ("Control", "Salt")
        }
        self.assertEqual(counts, expected)
        self.assertEqual({row["replicate"] for row in rows}, {"R1", "R2", "R3"})

    def test_consultor_has_six_condition_labels(self):
        labels = Counter(row["timepoint"] for row in read_tsv(SAMPLES))
        self.assertEqual(len(labels), 6)
        self.assertTrue(all(count == 3 for count in labels.values()))

    def test_metadata_is_accepted_by_current_consultor(self):
        sample_ids = [row["SampleID"] for row in read_tsv(SAMPLES)]
        loaded = load_metadata(SAMPLES, sample_ids)
        labels = Counter(row["timepoint"] for row in loaded)
        self.assertEqual(len(loaded), 18)
        self.assertEqual(len(labels), 6)
        self.assertTrue(all(count == 3 for count in labels.values()))

    def test_reference_is_pinned(self):
        lines = (PIPELINE / "manifests" / "reference.sha256").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 2)
        self.assertIn("2956f6f4c76c09eaef65d8732d9fa02d3682a0bff35a3b5027473d224f7beb32", lines[0])
        self.assertIn("61bdf7c0a607c723fa75c3b3139690489aedf23c0dcc230e15050b6edfc9f820", lines[1])

    def test_shell_scripts_are_strict(self):
        scripts = [PIPELINE / "run_pipeline.sh", *sorted((PIPELINE / "scripts").glob("*.sh"))]
        self.assertGreaterEqual(len(scripts), 10)
        for path in scripts:
            text = path.read_text(encoding="utf-8")
            self.assertTrue(text.startswith("#!/usr/bin/env bash"), path)
            self.assertIn("set -Eeuo pipefail", text, path)

    def test_source_tree_contains_no_raw_reads(self):
        forbidden = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".sra", ".sqlite")
        offending = [path for path in PIPELINE.rglob("*") if path.is_file() and path.name.lower().endswith(forbidden)]
        self.assertEqual(offending, [])


if __name__ == "__main__":
    unittest.main()
