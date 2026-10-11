"""
Automated tests for mutation_detective_final.py

Run them from the folder that holds both files:

    python3 -m unittest -v test_mutation_detective
    (or)  python3 test_mutation_detective.py

The tests use only the Python standard library. The picture tests are skipped
automatically if matplotlib is not installed.
"""

import contextlib
import io
import os
import random
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import mutation_detective_final as md

PROGRAM = md.__file__

try:
    import matplotlib  # noqa: F401
    HAVE_MATPLOTLIB = True
except ImportError:
    HAVE_MATPLOTLIB = False


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def category_of(normal, mutated, start=None):
    return md.analyze_sequences(normal, mutated, start=start)["category"]


def long_gene_pair():
    """A 1,200-base gene with a missense change at codon 10 and a 3-base deletion
    near codon 300, far apart so the alignment has to span the gap between them."""
    rng = random.Random(7)
    pool = [c for c, aa in md.CODON_TABLE.items() if aa != "*" and c != "ATG"]
    codons = ["ATG"] + [rng.choice(pool) for _ in range(398)] + ["TAA"]
    mutated = list(codons)
    mutated[10] = "GGG" if codons[10] != "GGG" else "GCT"
    del mutated[300]
    return "".join(codons), "".join(mutated)


def rebuild_mutated(ops, normal, mutated):
    """Use only the opcodes to rebuild the mutated sequence (checks they are valid)."""
    return "".join(normal[i1:i2] if tag == "equal" else mutated[j1:j2]
                   for tag, i1, i2, j1, j2 in ops)


def run_report(normal, mutated, **kwargs):
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        md.analyze_mutation(normal, mutated, **kwargs)
    return buffer.getvalue()


# --------------------------------------------------------------------------
# Cleaning and translation
# --------------------------------------------------------------------------
class CleanAndTranslateTests(unittest.TestCase):
    def test_clean_uppercases_and_strips_whitespace(self):
        self.assertEqual(md.clean_sequence(" atg gct\n"), "ATGGCT")

    def test_clean_converts_rna_to_dna(self):
        self.assertEqual(md.clean_sequence("AUGGCU"), "ATGGCT")

    def test_clean_rejects_invalid_letters(self):
        with self.assertRaises(ValueError):
            md.clean_sequence("ATGXCT")

    def test_clean_rejects_empty_input(self):
        for text in ("", "   \n"):
            with self.assertRaises(ValueError):
                md.clean_sequence(text)

    def test_genetic_code_is_complete(self):
        self.assertEqual(len(md.CODON_TABLE), 64)
        self.assertEqual(sorted(c for c, aa in md.CODON_TABLE.items() if aa == "*"),
                         ["TAA", "TAG", "TGA"])
        self.assertEqual(set(md.CODON_TABLE.values()) - {"*"}, set(md.AMINO_ACIDS) - {"*"})

    def test_translate_simple(self):
        self.assertEqual(md.translate_dna("ATGGAATTTTAA"), "MEF*")

    def test_translate_stops_at_first_stop_codon(self):
        self.assertEqual(md.translate_dna("ATGTAAGGGCCC"), "M*")

    def test_translate_ignores_leftover_bases(self):
        self.assertEqual(md.translate_dna("ATGGCTGA"), "MA")

    def test_translate_without_a_full_codon(self):
        self.assertEqual(md.translate_dna("AT"), "")


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------
class AlignmentTests(unittest.TestCase):
    def changes(self, ops):
        return [op for op in ops if op[0] != "equal"]

    def test_identical_sequences_have_no_changes(self):
        self.assertEqual(self.changes(md.align("ATGGCT", "ATGGCT")), [])

    def test_single_substitution(self):
        ops = self.changes(md.align("ATGGCTGAA", "ATGGCCGAA"))
        self.assertEqual(ops, [("replace", 5, 6, 5, 6)])

    def test_single_insertion(self):
        ops = self.changes(md.align("ATGGCTGAA", "ATGGCTAGAA"))
        self.assertEqual(len(ops), 1)
        self.assertEqual(ops[0][0], "insert")

    def test_single_deletion(self):
        ops = self.changes(md.align("ATGGCTGAAAAG", "ATGGCTGAAAG"))
        self.assertEqual(len(ops), 1)
        self.assertEqual(ops[0][0], "delete")

    def test_substitution_next_to_a_deletion_stays_two_changes(self):
        ops = self.changes(md.align("ATGGCTGAAAAGCTGTAA", "ATGGCCGAAAGCTGTAA"))
        self.assertEqual(sorted(op[0] for op in ops), ["delete", "replace"])

    def test_distant_changes_in_a_long_sequence(self):
        normal, mutated = long_gene_pair()
        ops = self.changes(md.align(normal, mutated))
        self.assertEqual(sorted(op[0] for op in ops), ["delete", "replace"])

    def test_opcodes_rebuild_the_mutated_sequence(self):
        rng = random.Random(1)
        for _ in range(150):
            normal = "".join(rng.choice("ACGT") for _ in range(rng.randint(30, 300)))
            letters = list(normal)
            for _ in range(rng.randint(1, 4)):
                position = rng.randrange(len(letters))
                kind = rng.choice("sid")
                if kind == "s":
                    letters[position] = rng.choice(
                        [b for b in "ACGT" if b != letters[position]])
                elif kind == "i":
                    letters.insert(position, rng.choice("ACGT"))
                elif len(letters) > 10:
                    del letters[position]
            mutated = "".join(letters)
            ops = md.align(normal, mutated)
            self.assertEqual(rebuild_mutated(ops, normal, mutated), mutated)
            self.assertEqual(sum(o[2] - o[1] for o in ops), len(normal))
            self.assertEqual(sum(o[4] - o[3] for o in ops), len(mutated))

    def test_fallback_to_difflib_for_large_problems(self):
        normal, mutated = "ATGGCTGAAAAGCTGTAA", "ATGGCCGAAAGCTGTAA"
        with mock.patch.object(md, "MAX_CELLS", 10):
            ops = md.align(normal, mutated)
        self.assertEqual(rebuild_mutated(ops, normal, mutated), mutated)


# --------------------------------------------------------------------------
# Reading frame
# --------------------------------------------------------------------------
class ReadingFrameTests(unittest.TestCase):
    def test_auto_detect_finds_the_start_codon(self):
        position, note = md.choose_start("CCATGGCTTAA")
        self.assertEqual(position, 2)
        self.assertIn("base 3", note)

    def test_auto_detect_prefers_the_longest_orf(self):
        position, _ = md.choose_start("ATGTAAATGGCTGAAAAGCTGTAA")
        self.assertEqual(position, 6)

    def test_no_atg_falls_back_to_base_one(self):
        position, note = md.choose_start("CCCGGGTTT")
        self.assertEqual(position, 0)
        self.assertIn("no ATG", note)

    def test_user_start_is_respected(self):
        position, note = md.choose_start("ATGTAAATGGCTGAAAAGCTGTAA", start=1)
        self.assertEqual(position, 0)
        self.assertIn("user", note)

    def test_user_start_that_is_not_atg_gives_a_warning(self):
        _, note = md.choose_start("CCCATGGCT", start=1)
        self.assertIn("warning", note)

    def test_user_start_out_of_range(self):
        for bad in (0, 50):
            with self.assertRaises(ValueError):
                md.choose_start("ATGGCT", start=bad)

    def test_map_position_shifts_after_an_upstream_insertion(self):
        normal, mutated = "GGATGGCTGAATAA", "GGCATGGCTGAATAA"
        ops = md.align(normal, mutated)
        self.assertEqual(md.map_position(ops, 2), 3)

    def test_locate_labels_regions(self):
        self.assertIn("before the start", md.locate(1, 3, 12))
        self.assertIn("after the stop", md.locate(12, 3, 12))
        self.assertEqual(md.locate(3, 3, 12), "codon 1")
        self.assertEqual(md.locate(9, 3, 12), "codon 3")


# --------------------------------------------------------------------------
# Classification: one case for every verdict
# --------------------------------------------------------------------------
CASES = [
    ("identical", "ATGGCTGAAAAGCTGTAA", "ATGGCTGAAAAGCTGTAA", "NONE", None),
    ("silent", "ATGGCTGAAAAGCTGTAA", "ATGGCCGAAAAGCTGTAA", "SILENT", None),
    ("missense (sickle-cell style)", "ATGGTGCACCTGACTCCTGAGGAGAAGTAA",
     "ATGGTGCACCTGACTCCTGTGGAGAAGTAA", "MISSENSE", None),
    ("nonsense", "ATGGCTGAAAAGCTGTAA", "ATGGCTTAAAAGCTGTAA", "NONSENSE", None),
    ("stop-loss with read-through", "ATGTAACCCTGA", "ATGCAACCCTGA", "STOP-LOSS", None),
    ("start loss", "ATGGCTGAATAA", "ATAGCTGAATAA", "START LOSS", None),
    ("frameshift (1-base deletion)", "ATGGCTGAAAAGCTGTAA", "ATGGCTGAAAGCTGTAA",
     "FRAMESHIFT", None),
    ("frameshift (1-base insertion)", "ATGGCTGAAAAGCTGTAA", "ATGGCTGGAAAAGCTGTAA",
     "FRAMESHIFT", None),
    ("in-frame deletion", "ATGGCTGAAAAGCTGTAA", "ATGGCTAAGCTGTAA", "IN-FRAME INDEL", None),
    ("in-frame insertion", "ATGGCTGAAAAGCTGTAA", "ATGGCTGAAGGGAAGCTGTAA",
     "IN-FRAME INDEL", None),
    ("cancelling indels", "ATGGCTGAAAAGCTGCCCTAA", "ATGGCTGAAAGCTGCCCGTAA",
     "TEMPORARY FRAMESHIFT", None),
    ("change after the stop codon", "ATGGCTTAAGGG", "ATGGCTTAAGGA", "SILENT", None),
    ("indel after the stop codon", "ATGGCTTAAGGGCCC", "ATGGCTTAAGGCCC",
     "NO PROTEIN CHANGE", None),
    ("insertion before the start codon", "GGATGGCTGAATAA", "GGCATGGCTGAATAA",
     "NO PROTEIN CHANGE", None),
    ("missense with a 5' leader", "CCATGGTGCACCTGACTCCTGAGGAGAAGTAA",
     "CCATGGTGCACCTGACTCCTGTGGAGAAGTAA", "MISSENSE", None),
    ("start loss with a 5' leader", "CCATGGCTGAATAA", "CCATAGCTGAATAA", "START LOSS", None),
    ("auto frame picks the longest ORF", "ATGTAAATGGCTGAAAAGCTGTAA",
     "ATGTATATGGCTGAAAAGCTGTAA", "NO PROTEIN CHANGE", None),
    ("--start overrides the auto frame", "ATGTAAATGGCTGAAAAGCTGTAA",
     "ATGTATATGGCTGAAAAGCTGTAA", "STOP-LOSS", 1),
    ("silent + missense together", "ATGGCTGAAAAGCTGTAA", "ATGGCCGAAAATCTGTAA",
     "MISSENSE", None),
    ("RNA input", "AUGGAAUUUUAA", "AUGGUAUUUUAA", "MISSENSE", None),
    ("long gene, distant changes", *long_gene_pair(), "IN-FRAME INDEL", None),
]


class ClassificationTests(unittest.TestCase):
    def test_every_verdict(self):
        for name, normal, mutated, expected, start in CASES:
            with self.subTest(name):
                self.assertEqual(category_of(normal, mutated, start), expected)

    def test_results_dictionary_is_complete(self):
        result = md.analyze_sequences("ATGGAATTTTAA", "ATGGTATTTTAA")
        for key in ("normal", "mutated", "ops", "start", "mutated_start", "frame_note",
                    "normal_protein", "mutated_protein", "coding_end", "category",
                    "effect"):
            self.assertIn(key, result)
        self.assertEqual(result["normal_protein"], "MEF*")
        self.assertEqual(result["mutated_protein"], "MVF*")

    def test_nonsense_shortens_the_protein(self):
        result = md.analyze_sequences("ATGGCTGAAAAGCTGTAA", "ATGGCTTAAAAGCTGTAA")
        self.assertEqual(result["mutated_protein"], "MA*")
        self.assertIn("premature stop", result["effect"])

    def test_frameshift_reports_where_the_protein_first_differs(self):
        result = md.analyze_sequences("ATGGCTGAAAAGCTGTAA", "ATGGCTGAAAGCTGTAA")
        self.assertIn("residue 4", result["effect"])


# --------------------------------------------------------------------------
# Text report
# --------------------------------------------------------------------------
class ReportTests(unittest.TestCase):
    def test_missense_report(self):
        text = run_report("ATGGAATTTTAA", "ATGGTATTTTAA")
        for expected in ("MISSENSE", "Substitution at base 5: A -> T", "[codon 2]",
                         "MEF*", "MVF*", "Residue 2: GAA (Glu) -> GTA (Val)",
                         "Protein changed? YES"):
            self.assertIn(expected, text)

    def test_marker_points_at_the_changed_amino_acid(self):
        lines = run_report("ATGGCTGAAAAGCTGTAA", "ATGGCTTAAAAGCTGTAA").splitlines()
        normal_line = next(l for l in lines if l.startswith("Normal protein:"))
        marker_line = lines[lines.index(normal_line) + 1]
        mutated_line = next(l for l in lines if l.startswith("Mutated protein:"))
        column = marker_line.index("^")           # first marked position
        self.assertEqual(normal_line[column], "E")  # Glu: the first residue that changed
        self.assertEqual(mutated_line[column], "*")  # now a STOP
        # both protein rows must start in the same column
        self.assertEqual(normal_line.index("MAEKL"), mutated_line.index("MA*"))

    def test_silent_report_says_protein_unchanged(self):
        text = run_report("ATGGCTGAAAAGCTGTAA", "ATGGCCGAAAAGCTGTAA")
        self.assertIn("Protein changed? NO", text)
        self.assertIn("SILENT", text)

    def test_identical_report(self):
        text = run_report("ATGGCTTAA", "ATGGCTTAA")
        self.assertIn("No DNA mutation detected.", text)

    def test_start_loss_mentions_the_next_atg(self):
        text = run_report("ATGGCTATGAAATAA", "ATAGCTATGAAATAA")
        self.assertIn("next ATG", text)

    def test_long_gene_report_is_windowed(self):
        normal, mutated = long_gene_pair()
        text = run_report(normal, mutated)
        self.assertIn("showing residues", text)
        self.assertIn("IN-FRAME", text)

    def test_invalid_sequence_raises(self):
        with self.assertRaises(ValueError):
            md.analyze_mutation("ATGXXX", "ATGGCT")


# --------------------------------------------------------------------------
# Picture
# --------------------------------------------------------------------------
@unittest.skipUnless(HAVE_MATPLOTLIB, "matplotlib is not installed")
class ImageTests(unittest.TestCase):
    def make_image(self, normal, mutated, **kwargs):
        folder = tempfile.mkdtemp()
        path = os.path.join(folder, "result.png")
        with contextlib.redirect_stdout(io.StringIO()):
            md.analyze_mutation(normal, mutated, image=path, **kwargs)
        self.addCleanup(lambda: (os.path.exists(path) and os.remove(path), os.rmdir(folder)))
        return path

    def assert_valid_png(self, path):
        self.assertTrue(os.path.exists(path))
        self.assertGreater(os.path.getsize(path), 5000)
        with open(path, "rb") as handle:
            self.assertEqual(handle.read(8), b"\x89PNG\r\n\x1a\n")

    def test_missense_image(self):
        self.assert_valid_png(self.make_image("ATGGAATTTTAA", "ATGGTATTTTAA"))

    def test_frameshift_image_with_a_gap(self):
        self.assert_valid_png(self.make_image("ATGGCTGAAAAGCTGTAA", "ATGGCTGAAAGCTGTAA"))

    def test_nonsense_image_with_a_shorter_protein(self):
        self.assert_valid_png(self.make_image("ATGGCTGAAAAGCTGTAA", "ATGGCTTAAAAGCTGTAA"))

    def test_identical_sequences_image(self):
        self.assert_valid_png(self.make_image("ATGGCTTAA", "ATGGCTTAA"))

    def test_long_gene_image(self):
        self.assert_valid_png(self.make_image(*long_gene_pair()))

    def test_no_protein_image(self):
        self.assert_valid_png(self.make_image("CCCGG", "CCCGA"))


class MissingMatplotlibTests(unittest.TestCase):
    def test_program_keeps_working_without_matplotlib(self):
        results = md.analyze_sequences("ATGGAATTTTAA", "ATGGTATTTTAA")
        buffer = io.StringIO()
        with mock.patch.dict(sys.modules, {"matplotlib": None}), \
                contextlib.redirect_stdout(buffer):
            saved = md.save_comparison_image(results, "never_written.png")
        self.assertFalse(saved)
        self.assertIn("matplotlib is not installed", buffer.getvalue())
        self.assertFalse(os.path.exists("never_written.png"))


# --------------------------------------------------------------------------
# Command line
# --------------------------------------------------------------------------
class CommandLineTests(unittest.TestCase):
    def run_program(self, *args, input_text=None):
        return subprocess.run([sys.executable, PROGRAM, *args], capture_output=True,
                              text=True, input=input_text, timeout=60)

    def test_arguments(self):
        result = self.run_program("ATGGAATTTTAA", "ATGGTATTTTAA")
        self.assertEqual(result.returncode, 0)
        self.assertIn("MISSENSE", result.stdout)

    def test_start_option(self):
        result = self.run_program("ATGTAAATGGCTGAAAAGCTGTAA",
                                  "ATGTATATGGCTGAAAAGCTGTAA", "--start", "1")
        self.assertIn("STOP-LOSS", result.stdout)

    def test_invalid_sequence_is_reported_politely(self):
        result = self.run_program("ATGXXX", "ATGGCT")
        self.assertIn("Input error", result.stdout)
        self.assertNotIn("Traceback", result.stderr)

    def test_interactive_mode(self):
        result = self.run_program(input_text="ATGGAATTTTAA\nATGGTATTTTAA\n\n")
        self.assertIn("MISSENSE", result.stdout)

    @unittest.skipUnless(HAVE_MATPLOTLIB, "matplotlib is not installed")
    def test_image_option(self):
        with tempfile.TemporaryDirectory() as folder:
            path = os.path.join(folder, "out.png")
            result = self.run_program("ATGGAATTTTAA", "ATGGTATTTTAA", "--image", path)
            self.assertIn("Image saved", result.stdout)
            self.assertTrue(os.path.exists(path))


if __name__ == "__main__":
    unittest.main(verbosity=2)
