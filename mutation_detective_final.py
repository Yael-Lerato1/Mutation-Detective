"""
MUTATION DETECTIVE
==================
Compares a normal DNA sequence with a mutated copy, finds the mutation, translates
both into proteins (chains of amino acids) and tells you whether the mutation
changes the protein.

HOW TO RUN IT
-------------
    python mutation_detective_final.py                     # asks you to type the sequences
    python mutation_detective_final.py NORMAL MUTATED      # sequences given directly
    python mutation_detective_final.py NORMAL MUTATED --start 3
    python mutation_detective_final.py NORMAL MUTATED --image result.png

    Example:
    python mutation_detective_final.py ATGGCTGAAAAGCTGTAA ATGGCTTAAAAGCTGTAA

Sequences may use the letters A, T, C and G (RNA 'U' is accepted and converted to T).
--start N is optional: it is the position of the start codon in the normal sequence.
If you leave it out, the program finds the start codon itself.
--image FILE.png is optional: it also saves a picture of the comparison. This is the
only feature that needs an extra package (python3 -m pip install matplotlib).

WHAT HAPPENS, STEP BY STEP
--------------------------
When you run the program, main() calls analyze_mutation(), which runs
analyze_sequences() (steps 1-5), then print_report() (step 6) and, if you asked for
a picture, save_comparison_image() (step 7):

 1. CLEAN      clean_sequence() upper-cases the text and rejects anything that is not A, T, C or G.
 2. ALIGN      align() lines the two sequences up, like two rows of beads, and marks every difference as a substitution (one base swapped), an insertion
               (extra bases) or a deletion (missing bases). It first skips the parts that already match, then uses a scoring method (_gotoh) that prefers
               one small insertion/deletion over a messy pile of substitutions.
 3. FIND FRAME choose_start() looks for the start codon (ATG), because a protein is read in groups of three bases (codons) beginning there. It picks the
               start of the longest open reading frame, or uses your --start value. 
               map_position() finds the same spot in the mutated sequence.
 4. TRANSLATE  translate_dna() reads each codon with the genetic code table (CODON_TABLE) and builds the protein until it reaches a STOP codon.This is done for both sequences.
 5. CLASSIFY   classify_effect() compares the two proteins and names the effect.
 6. REPORT     print_report() prints what changed and where, both proteins with a ^marker under each difference, the changed codons, and the verdict.
 7. PICTURE    (optional, --image) save_comparison_image() draws the aligned DNA with the changed bases highlighted, the two proteins amino acid by amino acid, and the verdict, and saves it as an image file.

THE VERDICTS YOU MAY SEE
------------------------
 SILENT                 DNA changed, but the amino acids did not (several codons code for the same amino acid).
 MISSENSE               One amino acid is swapped for a different one.
 NONSENSE               A codon turns into a STOP, so the protein is cut short.
 STOP-LOSS              The STOP codon is destroyed, so the protein runs on too long.
 START LOSS             The start codon (ATG) is destroyed.
 FRAMESHIFT             Bases were added or removed, not in multiples of 3, so every codon after that point is read wrongly.
 TEMPORARY FRAMESHIFT   An insertion and a deletion cancel out in length, but the codons between them are misread.
 IN-FRAME INDEL         Bases were added or removed in multiples of 3, so amino acids are gained or lost but the rest stays correct.
 NO PROTEIN CHANGE      DNA changed, but outside the part that makes the protein (before the start codon or after the stop codon).
 NONE                   The sequences are identical.

HOW THE FILE IS ORGANISED (top to bottom)
-----------------------------------------
 - Data tables:        CODON_TABLE (codon -> amino acid), AMINO_ACIDS (letter -> name) and the scoring settings used for aligning.
 - Basic helpers:      clean_sequence, translate_dna
 - Alignment:          align and its helpers (_gotoh, _columns_to_ops)
 - Reading frame:      choose_start, map_position, cds_ops, locate
 - Describe/classify:  describe_mutations, classify_effect and small helpers
 - Analysis:           analyze_sequences does steps 1-5 and returns a dictionary
 - Report:             print_report prints the results (step 6)
 - Picture:            ops_to_columns, save_comparison_image (step 7, optional)
 - analyze_mutation:   runs the analysis, the report and the optional picture
 - main():             reads the command line (or asks you to type) and runs it all

LIMITATIONS
-----------
Reads the forward DNA strand only, uses the standard genetic code, and treats the
sequence as one continuous coding region (no introns or splicing).
"""
import argparse
import textwrap
from difflib import SequenceMatcher

# Standard genetic code: DNA codons -> amino acids
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMINO_ACIDS = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp",
    "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "STOP",
}

# Alignment scoring (affine gaps: opening a gap costs more than extending it)
MATCH, MISMATCH = 2, -3
GAP_OPEN, GAP_EXTEND = -6, -1
BAND = 25               # alignment band around the diagonal
MAX_CELLS = 1_500_000   # beyond this, fall back to difflib
NEG = -10 ** 9


# --------------------------------------------------------------------------
# Basic helpers
# --------------------------------------------------------------------------
def clean_sequence(sequence):
    """Normalize and validate a DNA sequence (RNA 'U' is converted to 'T')."""
    sequence = "".join(sequence.upper().split()).replace("U", "T")

    if not sequence:
        raise ValueError("DNA sequence cannot be empty.")

    if any(base not in "ATCG" for base in sequence):
        raise ValueError("Use only DNA bases A, T, C and G.")

    return sequence


def translate_dna(sequence):
    """Translate DNA from the first base in groups of three, up to the first stop."""
    protein = []

    for i in range(0, len(sequence) - 2, 3):
        amino_acid = CODON_TABLE[sequence[i:i + 3]]
        protein.append(amino_acid)
        if amino_acid == "*":
            break

    return "".join(protein)


# --------------------------------------------------------------------------
# Alignment (affine-gap, banded Gotoh)
# --------------------------------------------------------------------------
def _columns_to_ops(cols, i, j):
    """Group alignment columns into difflib-style opcodes."""
    ops = []
    for ca, cb in cols:
        if ca == "-":
            tag = "insert"
        elif cb == "-":
            tag = "delete"
        elif ca == cb:
            tag = "equal"
        else:
            tag = "replace"
        di = 0 if ca == "-" else 1
        dj = 0 if cb == "-" else 1
        if ops and ops[-1][0] == tag:
            t, i1, i2, j1, j2 = ops[-1]
            ops[-1] = (t, i1, i2 + di, j1, j2 + dj)
        else:
            ops.append((tag, i, i + di, j, j + dj))
        i += di
        j += dj
    return ops


def _gotoh(x, y):
    """Banded global alignment with affine gap penalties.

    Returns a list of columns (char_from_x or '-', char_from_y or '-'),
    or None if the problem is too large for the band approach.
    """
    n, m = len(x), len(y)
    if (n + 1) * (abs(m - n) + 2 * BAND + 1) > MAX_CELLS:
        return None

    lo = min(0, m - n) - BAND
    hi = max(0, m - n) + BAND

    M, X, Y = {(0, 0): 0}, {}, {}      # score matrices (match / delete / insert)
    pM, pX, pY = {}, {}, {}            # back-pointers (0=M, 1=X, 2=Y)

    def best(cands):
        # max() keeps the first of equal values, so ties prefer M, then X, then Y
        return max(cands, key=lambda c: c[0])

    for i in range(n + 1):
        for j in range(max(0, i + lo), min(m, i + hi) + 1):
            if i == 0 and j == 0:
                continue
            if i > 0 and j > 0:
                s = MATCH if x[i - 1] == y[j - 1] else MISMATCH
                v, st = best([(M.get((i - 1, j - 1), NEG), 0),
                              (X.get((i - 1, j - 1), NEG), 1),
                              (Y.get((i - 1, j - 1), NEG), 2)])
                if v > NEG // 2:
                    M[(i, j)] = v + s
                    pM[(i, j)] = st
            if i > 0:   # x[i-1] against a gap: a deletion in the mutated sequence
                v, st = best([(M.get((i - 1, j), NEG) + GAP_OPEN, 0),
                              (X.get((i - 1, j), NEG) + GAP_EXTEND, 1),
                              (Y.get((i - 1, j), NEG) + GAP_OPEN, 2)])
                if v > NEG // 2:
                    X[(i, j)] = v
                    pX[(i, j)] = st
            if j > 0:   # gap against y[j-1]: an insertion in the mutated sequence
                v, st = best([(M.get((i, j - 1), NEG) + GAP_OPEN, 0),
                              (Y.get((i, j - 1), NEG) + GAP_EXTEND, 2),
                              (X.get((i, j - 1), NEG) + GAP_OPEN, 1)])
                if v > NEG // 2:
                    Y[(i, j)] = v
                    pY[(i, j)] = st

    _, state = best([(M.get((n, m), NEG), 0), (X.get((n, m), NEG), 1),
                     (Y.get((n, m), NEG), 2)])
    cols, i, j = [], n, m
    while i > 0 or j > 0:
        if state == 0:
            cols.append((x[i - 1], y[j - 1]))
            state = pM[(i, j)]
            i, j = i - 1, j - 1
        elif state == 1:
            cols.append((x[i - 1], "-"))
            state = pX[(i, j)]
            i -= 1
        else:
            cols.append(("-", y[j - 1]))
            state = pY[(i, j)]
            j -= 1
    cols.reverse()
    return cols


def align(a, b):
    """Return difflib-style opcodes describing how to turn a into b.

    The shared prefix and suffix are trimmed first (so an indel in a repeat is
    reported at its 3'-most position), and only the differing core is aligned.
    """
    limit = min(len(a), len(b))
    p = 0
    while p < limit and a[p] == b[p]:
        p += 1
    s = 0
    while s < limit - p and a[-1 - s] == b[-1 - s]:
        s += 1

    x, y = a[p:len(a) - s], b[p:len(b) - s]
    ops = []
    if p:
        ops.append(("equal", 0, p, 0, p))

    if x or y:
        if not x or not y:
            cols = [(c, "-") for c in x] + [("-", c) for c in y]
            ops += _columns_to_ops(cols, p, p)
        else:
            cols = _gotoh(x, y)
            if cols is None:  # too large: fall back to difflib
                ops += [(t, i1 + p, i2 + p, j1 + p, j2 + p)
                        for t, i1, i2, j1, j2
                        in SequenceMatcher(None, x, y, autojunk=False).get_opcodes()]
            else:
                ops += _columns_to_ops(cols, p, p)

    if s:
        ops.append(("equal", len(a) - s, len(a), len(b) - s, len(b)))
    return ops


# --------------------------------------------------------------------------
# Reading frame
# --------------------------------------------------------------------------
def choose_start(normal, start=None):
    """Pick the start codon position (0-based) and describe how it was chosen."""
    if start is not None:
        if not 1 <= start <= len(normal) - 2:
            raise ValueError(f"--start must be between 1 and {len(normal) - 2}.")
        note = f"start set by user at base {start}"
        if normal[start - 1:start + 2] != "ATG":
            note += f" (warning: that codon is {normal[start - 1:start + 2]}, not ATG)"
        return start - 1, note

    best_pos, best_len = None, -1
    pos = normal.find("ATG")
    while pos != -1:
        length = len(translate_dna(normal[pos:]))
        if length > best_len:        # strict '>' keeps the earliest on ties
            best_pos, best_len = pos, length
        pos = normal.find("ATG", pos + 1)

    if best_pos is None:
        return 0, "no ATG found, so reading from base 1"
    return best_pos, f"longest open reading frame starts at base {best_pos + 1}"


def map_position(ops, pos):
    """Position in the mutated sequence that corresponds to normal[pos]."""
    for tag, i1, i2, j1, j2 in ops:
        if tag != "insert" and i1 <= pos < i2:
            if tag == "equal" or (tag == "replace" and i2 - i1 == j2 - j1):
                return j1 + (pos - i1)
            return min(j1 + (pos - i1), j2)   # that base itself was deleted/rewritten
    return ops[-1][4] if ops else 0


def cds_ops(ops, s, j):
    """Re-express opcodes in coordinates of normal[s:] and mutated[j:]."""
    out = []
    for tag, i1, i2, j1, j2 in ops:
        if i2 <= s and j2 <= j:
            continue                           # wholly before the start codon
        ni1, ni2 = max(i1, s) - s, max(i2, s) - s
        nj1, nj2 = max(j1, j) - j, max(j2, j) - j
        if tag == "equal" and ni1 == ni2:
            continue
        out.append((tag, ni1, ni2, nj1, nj2))
    return out


def locate(pos, start, coding_end):
    """Describe where a base position lies relative to the coding region."""
    if pos < start:
        return "before the start codon (5' non-coding)"
    if pos >= coding_end:
        return "after the stop codon (3' non-coding)"
    return f"codon {(pos - start) // 3 + 1}"


# --------------------------------------------------------------------------
# Describing and classifying
# --------------------------------------------------------------------------
def describe_mutations(ops, normal, mutated, start, coding_end):
    """Turn opcodes into readable descriptions with their location."""
    descriptions = []

    for tag, i1, i2, j1, j2 in ops:
        if tag == "equal":
            continue

        old, new = normal[i1:i2], mutated[j1:j2]
        where_pos = max(i1 - 1, 0) if tag == "insert" else i1
        place = locate(where_pos, start, coding_end)

        if tag == "replace" and len(old) == len(new):
            where = f"base {i1 + 1}" if len(old) == 1 else f"bases {i1 + 1}-{i2}"
            text = f"Substitution at {where}: {old} -> {new}"
        elif tag == "replace":
            text = (f"Replacement at bases {i1 + 1}-{i2}: {old} -> {new} "
                    f"({len(old)} bases replaced by {len(new)})")
        elif tag == "insert":
            text = f"Insertion after base {i1}: {new} ({len(new)} base(s))"
        else:
            where = f"base {i1 + 1}" if len(old) == 1 else f"bases {i1 + 1}-{i2}"
            text = f"Deletion at {where}: {old} ({len(old)} base(s))"
        descriptions.append(f"{text}  [{place}]")

    return descriptions or ["No DNA mutation detected."]


def frame_offsets(opcodes):
    """Reading-frame offset (mutated minus normal length so far) after each change."""
    return [j2 - i2 for tag, i1, i2, j1, j2 in opcodes if tag != "equal"]


def has_length_change(opcodes):
    return any((j2 - j1) != (i2 - i1)
               for tag, i1, i2, j1, j2 in opcodes if tag != "equal")


def first_difference(a, b):
    """Index of the first position where two strings differ, or None."""
    for i in range(min(len(a), len(b))):
        if a[i] != b[i]:
            return i
    return None if len(a) == len(b) else min(len(a), len(b))


def marker_line(a, b):
    """'^' under each position where two protein strings differ."""
    return "".join(
        "^" if i >= len(a) or i >= len(b) or a[i] != b[i] else " "
        for i in range(max(len(a), len(b)))
    )


def codon_changes(normal, mutated, normal_protein):
    """Changed codons within the coding region (same-length comparisons)."""
    changes = []
    for i in range(len(normal_protein)):
        a, b = normal[i * 3:i * 3 + 3], mutated[i * 3:i * 3 + 3]
        if len(b) == 3 and a != b:
            changes.append((i + 1, a, b))
    return changes


def classify_effect(normal, mutated, normal_protein, mutated_protein, opcodes):
    """Return (category, explanation). Inputs start at the start codon."""
    if normal == mutated:
        return "NONE", "No mutation: the DNA sequences are identical."

    offsets = frame_offsets(opcodes)
    indel = has_length_change(opcodes)
    net_shift = (len(mutated) - len(normal)) % 3 != 0
    transient_shift = (not net_shift) and any(o % 3 != 0 for o in offsets)

    # 1. The protein is unchanged
    if normal_protein == mutated_protein:
        if not indel and not codon_changes(normal, mutated, normal_protein):
            return "SILENT", (
                "SILENT: the DNA changed, but only outside the coding region "
                "(after the stop codon), so the protein is identical."
            )
        if not indel:
            return "SILENT", (
                "SILENT: the DNA changed, but every altered codon still codes for "
                "the same amino acid, so the protein is identical."
            )
        return "NO PROTEIN CHANGE", (
            "NO PROTEIN CHANGE: the DNA changed, but the translated protein is "
            "identical (the change lies after the stop codon)."
        )

    # 2. Start codon destroyed
    if normal[:3] == "ATG" and mutated[:3] != "ATG":
        return "START LOSS", (
            "START LOSS: the start codon (ATG, Met) is altered, so translation "
            "may not begin at the normal place and the protein may not be made."
        )

    # 3. Frameshift (permanent or temporary)
    if net_shift:
        d = first_difference(normal_protein, mutated_protein)
        return "FRAMESHIFT", (
            "FRAMESHIFT: the net change in length is not a multiple of 3, so every "
            f"codon after the mutation is read in a shifted frame (protein first "
            f"differs at residue {d + 1})."
        )
    if transient_shift:
        return "TEMPORARY FRAMESHIFT", (
            "TEMPORARY FRAMESHIFT: separate insertions/deletions cancel out in total "
            "length, but the codons between them are read in the wrong frame, "
            "altering that stretch of the protein."
        )

    # 4. Premature stop codon (nonsense)
    if "*" in mutated_protein:
        n_stop = normal_protein.find("*")
        m_stop = mutated_protein.find("*")
        # An in-frame indel moves the normal stop by the number of codons gained
        # or lost; only a stop earlier than that shifted position is premature.
        shift = (len(mutated) - len(normal)) // 3
        if n_stop == -1 or m_stop < n_stop + shift:
            return "NONSENSE", (
                f"NONSENSE: a premature stop codon appears at residue {m_stop + 1}, "
                f"shortening the protein from {len(normal_protein.rstrip('*'))} to "
                f"{len(mutated_protein.rstrip('*'))} amino acids."
            )

    # 5. Stop codon lost (read-through)
    n_stop = normal_protein.find("*")
    if (n_stop != -1 and len(mutated_protein) > n_stop
            and mutated_protein[n_stop] != "*"
            and mutated_protein[:n_stop] == normal_protein[:n_stop]):
        return "STOP-LOSS", (
            "STOP-LOSS: the original stop codon is changed, so translation reads "
            f"on past it and extends the protein to {len(mutated_protein.rstrip('*'))} "
            f"amino acids (normal: {len(normal_protein.rstrip('*'))})."
        )

    # 6. In-frame insertion/deletion
    if indel:
        diff = (len(mutated) - len(normal)) // 3
        word = "added" if diff > 0 else "removed"
        return "IN-FRAME INDEL", (
            f"IN-FRAME INSERTION/DELETION: the reading frame is kept, and about "
            f"{abs(diff)} amino acid(s) are {word}."
        )

    # 7. Plain amino acid substitution(s)
    return "MISSENSE", (
        "MISSENSE: one or more amino acids are replaced by different ones. "
        "The protein keeps its length but may fold or work differently."
    )


# --------------------------------------------------------------------------
# Main analysis and report
# --------------------------------------------------------------------------
def shorten(seq, width=80):
    return seq if len(seq) <= width else f"{seq[:width]}... ({len(seq)} bases)"


def print_protein_comparison(normal_protein, mutated_protein, window=100):
    """Print both proteins with a '^' marker line; window long proteins."""
    d = first_difference(normal_protein, mutated_protein)
    start = 0
    if max(len(normal_protein), len(mutated_protein)) > window and d is not None:
        start = max(0, d - 20)
    end = start + window
    a, b = normal_protein[start:end], mutated_protein[start:end]
    if start or max(len(normal_protein), len(mutated_protein)) > end:
        print(f"(showing residues {start + 1}-{start + max(len(a), len(b))})")
    print("Normal protein: ", a or "(no full codon)")
    print(" " * 17 + marker_line(a, b))
    print("Mutated protein:", b or "(no full codon)")


def analyze_sequences(normal_input, mutated_input, start=None):
    """Run the whole analysis and return the results as a dictionary (nothing is printed)."""
    normal = clean_sequence(normal_input)
    mutated = clean_sequence(mutated_input)

    ops = align(normal, mutated)
    s, frame_note = choose_start(normal, start)
    j = map_position(ops, s)

    normal_cds, mutated_cds = normal[s:], mutated[j:]
    normal_protein = translate_dna(normal_cds)
    mutated_protein = translate_dna(mutated_cds)
    coding_end = s + len(normal_protein) * 3 if normal_protein.endswith("*") else len(normal)

    if normal == mutated:
        category, effect = "NONE", "No mutation: the DNA sequences are identical."
    elif normal_cds == mutated_cds:
        category = "NO PROTEIN CHANGE"
        effect = ("NO PROTEIN CHANGE: the DNA changed, but only before the start codon "
                  "(5' non-coding region), so the protein is identical.")
    else:
        category, effect = classify_effect(
            normal_cds, mutated_cds, normal_protein, mutated_protein,
            cds_ops(ops, s, j))

    return {
        "normal": normal, "mutated": mutated, "ops": ops,
        "start": s, "mutated_start": j, "frame_note": frame_note,
        "normal_cds": normal_cds, "mutated_cds": mutated_cds,
        "normal_protein": normal_protein, "mutated_protein": mutated_protein,
        "coding_end": coding_end, "category": category, "effect": effect,
    }


def print_report(r):
    """Print the text report for the results of analyze_sequences()."""
    normal, mutated = r["normal"], r["mutated"]
    normal_cds, mutated_cds = r["normal_cds"], r["mutated_cds"]
    normal_protein, mutated_protein = r["normal_protein"], r["mutated_protein"]
    category = r["category"]

    print("\n=== MUTATION DETECTIVE ===")
    print("Normal DNA:  ", shorten(normal))
    print("Mutated DNA: ", shorten(mutated))
    print(f"Reading frame: {r['frame_note']}")

    print("\n--- Mutation identification ---")
    for mutation in describe_mutations(r["ops"], normal, mutated, r["start"], r["coding_end"]):
        print("-", mutation)

    print("\n--- Protein translation (* = STOP) ---")
    print_protein_comparison(normal_protein, mutated_protein)

    if max(len(normal_protein), len(mutated_protein)) <= 40:
        print("\nNormal amino acids:")
        print("-".join(AMINO_ACIDS[aa] for aa in normal_protein) or "(none)")
        print("\nMutated amino acids:")
        print("-".join(AMINO_ACIDS[aa] for aa in mutated_protein) or "(none)")

    if len(normal_cds) == len(mutated_cds):
        changes = codon_changes(normal_cds, mutated_cds, normal_protein)
        if changes:
            print("\n--- Codon-level changes ---")
            for num, a, b in changes:
                aa_a, aa_b = CODON_TABLE[a], CODON_TABLE[b]
                tag = "same amino acid" if aa_a == aa_b else "amino acid changes"
                print(f"Residue {num}: {a} ({AMINO_ACIDS[aa_a]}) -> "
                      f"{b} ({AMINO_ACIDS[aa_b]})  [{tag}]")
    elif category in ("IN-FRAME INDEL", "FRAMESHIFT", "TEMPORARY FRAMESHIFT"):
        d = first_difference(normal_protein, mutated_protein)
        if d is not None:
            print(f"\nProtein first differs at residue {d + 1}.")

    if category == "START LOSS":
        k = mutated_cds.find("ATG")
        if k == -1:
            print("\nNo other ATG exists in the mutated sequence, so no protein "
                  "would be made.")
        else:
            alt = translate_dna(mutated_cds[k:])
            print(f"\nIf translation restarted at the next ATG (base "
                  f"{r['mutated_start'] + k + 1} of the mutated sequence), the protein "
                  f"would be {len(alt.rstrip('*'))} amino acids long "
                  f"(normal: {len(normal_protein.rstrip('*'))}).")

    print("\n--- Effect on protein ---")
    print(f"Protein changed? {'NO' if normal_protein == mutated_protein else 'YES'}")
    print(r["effect"])

    if len(normal_cds) % 3 != 0 or len(mutated_cds) % 3 != 0:
        print("\nNote: at least one coding region has leftover bases that do not "
              "form a complete codon.")

    print("\nNote: this analysis uses the standard genetic code and a forward-strand "
          "reading frame.")


# --------------------------------------------------------------------------
# Picture of the comparison (optional, needs matplotlib)
# --------------------------------------------------------------------------
def ops_to_columns(ops, normal, mutated):
    """Aligned columns as (normal base, mutated base, normal index, mutated index).

    A missing base is shown as '-' with index None.
    """
    cols = []
    for tag, i1, i2, j1, j2 in ops:
        a, b = normal[i1:i2], mutated[j1:j2]
        for k in range(max(len(a), len(b))):
            cols.append((
                a[k] if k < len(a) else "-",
                b[k] if k < len(b) else "-",
                i1 + k if k < len(a) else None,
                j1 + k if k < len(b) else None,
            ))
    return cols


def save_comparison_image(r, path):
    """Draw the aligned DNA, the two proteins and the verdict, and save as an image.

    Returns True on success, False if matplotlib is not installed.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
    except ImportError:
        print("\nImage not created: matplotlib is not installed. "
              "Install it with:  python3 -m pip install matplotlib")
        return False

    normal, mutated, ops = r["normal"], r["mutated"], r["ops"]
    s, j = r["start"], r["mutated_start"]
    np_, mp_ = r["normal_protein"], r["mutated_protein"]
    category = r["category"]

    # Colours
    tints = ("#D6E6F7", "#EAF2FB")       # alternating codon shading
    noncoding, highlight, edge = "#F2F2F2", "#F6B8B0", "#C0392B"
    neutral, stop_fill = "#E4EDF6", "#C9CED4"
    verdict_colors = {"SILENT": "#2E7D32", "NO PROTEIN CHANGE": "#2E7D32",
                      "NONE": "#555555", "MISSENSE": "#D9822B"}
    verdict_color = verdict_colors.get(category, edge)

    # ---- which part of the DNA to show (at most 48 columns) ----
    cols = ops_to_columns(ops, normal, mutated)
    first = next((k for k, c in enumerate(cols) if c[0] != c[1]), None)
    max_cols = 48
    lo = 0
    if len(cols) > max_cols:
        lo = min(max(0, (first or 0) - 20), len(cols) - max_cols)
    hi = min(len(cols), lo + max_cols)
    window = cols[lo:hi]

    # ---- which proteins residues to show (at most 16) ----
    d = first_difference(np_, mp_)
    total_res = max(len(np_), len(mp_))
    max_res = 16
    plo = 0
    if total_res > max_res:
        plo = min(max(0, (d or 0) - 5), total_res - max_res)
    phi = min(total_res, plo + max_res)
    n_res = phi - plo

    # ---- figure geometry: every axis uses the same inches-per-unit scale ----
    label = 5.0
    units_dna = len(window) + label + 0.5
    units_prot = max(n_res, 3) * 3 + label
    total_units = max(units_dna, units_prot)
    fig_w = max(8.0, total_units * 0.24)
    unit_in = fig_w * 0.96 / total_units
    fs_dna = max(7, min(14, unit_in * 72 * 0.75))
    fs_res = max(8, min(18, 2.8 * unit_in * 72 * 0.45))

    # ---- verdict text ----
    wrap_width = max(40, int(fig_w * 0.96 / 0.078))
    head, sep, tail = r["effect"].partition(": ")
    effect_text = tail if sep and head == head.upper() else r["effect"]
    effect_text = effect_text[:1].upper() + effect_text[1:]
    effect_lines = textwrap.wrap(effect_text, width=wrap_width)
    mutation_lines = describe_mutations(ops, normal, mutated, s, r["coding_end"])
    extra = len(mutation_lines) - 4
    mutation_lines = [textwrap.shorten(m, width=wrap_width) for m in mutation_lines[:4]]
    if extra > 0:
        mutation_lines.append(f"... and {extra} more")

    title_h, dna_h, prot_h = 0.85, 1.9, 2.3
    verdict_h = 1.05 + 0.25 * len(effect_lines) + 0.4 + 0.25 * len(mutation_lines)
    fig_h = title_h + dna_h + prot_h + verdict_h

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

    def add_axes(top_in, height_in, units):
        axes = fig.add_axes([0.02, 1 - (top_in + height_in) / fig_h,
                             units / total_units * 0.96, height_in / fig_h])
        axes.set_xlim(0, units)
        axes.axis("off")
        return axes

    fig.text(0.02, 1 - 0.35 / fig_h, "Mutation Detective", fontsize=17,
             fontweight="bold", va="center")
    fig.text(0.02, 1 - 0.65 / fig_h, f"Reading frame: {r['frame_note']}",
             fontsize=9, color="#666666", va="center")

    # ---- DNA panel ----
    ax = add_axes(title_h, dna_h, units_dna)
    ax.set_ylim(0, 4)
    caption = "DNA, aligned"
    if len(cols) > max_cols:
        caption += f" (columns {lo + 1}-{hi} of {len(cols)})"
    ax.text(0, 3.85, caption, fontsize=10, fontweight="bold", va="center")
    ax.text(label - 0.5, 2.7, "Normal", ha="right", va="center", fontsize=10)
    ax.text(label - 0.5, 1.3, "Mutated", ha="right", va="center", fontsize=10)

    for x, (ca, cb, ia, jb) in enumerate(window):
        changed = ca != cb
        for base, idx, y, seq_start, coding_len in (
                (ca, ia, 2.2, s, len(np_) * 3), (cb, jb, 0.8, j, len(mp_) * 3)):
            k = None if idx is None else idx - seq_start
            if changed:
                fill = highlight
            elif k is None:
                fill = "white"
            elif 0 <= k < coding_len:
                fill = tints[(k // 3) % 2]
            else:
                fill = noncoding
            ax.add_patch(Rectangle((label + x, y), 1, 1, facecolor=fill,
                                   edgecolor=edge if changed else "white",
                                   linewidth=1.6 if changed else 1.2,
                                   zorder=3 if changed else 2))
            ax.text(label + x + 0.5, y + 0.5, base, ha="center", va="center",
                    family="monospace", fontsize=fs_dna,
                    fontweight="bold" if changed else "normal",
                    color="#222222" if base != "-" else "#999999")
            if (y == 2.2 and k is not None and 0 <= k < coding_len and k % 3 == 0):
                ax.text(label + x + 0.5, 3.45, str(k // 3 + 1), ha="center",
                        va="center", fontsize=7, color="#777777")
    ax.text(0, 0.3, "red = changed base     '-' = missing base     "
            "blue shading = codons (small numbers count them)",
            fontsize=8, color="#666666", va="center")

    # ---- Protein panel ----
    ax = add_axes(title_h + dna_h, prot_h, units_prot)
    ax.set_ylim(0, 4.6)
    ax.text(0, 4.35, "Protein, amino acid by amino acid", fontsize=10,
            fontweight="bold", va="center")
    ax.text(label - 0.5, 2.8, "Normal", ha="right", va="center", fontsize=10)
    ax.text(label - 0.5, 1.1, "Mutated", ha="right", va="center", fontsize=10)

    if n_res == 0:
        ax.text(label, 2.0, "(no protein: no complete codon)", fontsize=10, color="#666666")
    for i in range(n_res):
        idx = plo + i
        a = np_[idx] if idx < len(np_) else None
        b = mp_[idx] if idx < len(mp_) else None
        x = label + i * 3
        ax.text(x + 1.4, 3.75, str(idx + 1), ha="center", va="center",
                fontsize=8, color="#777777")
        for aa, y in ((a, 2.2), (b, 0.5)):
            if aa is None:
                ax.add_patch(Rectangle((x, y), 2.8, 1.2, fill=False,
                                       edgecolor="#999999", linestyle="--"))
                ax.text(x + 1.4, y + 0.6, "none", ha="center", va="center",
                        fontsize=8, color="#999999")
                continue
            if a != b:
                fill = highlight
            else:
                fill = stop_fill if aa == "*" else neutral
            ax.add_patch(Rectangle((x, y), 2.8, 1.2, facecolor=fill,
                                   edgecolor=edge if a != b else "white", linewidth=1.6))
            ax.text(x + 1.4, y + 0.78, aa, ha="center", va="center",
                    fontsize=fs_res, fontweight="bold", family="monospace")
            ax.text(x + 1.4, y + 0.25, AMINO_ACIDS[aa], ha="center", va="center",
                    fontsize=max(6, fs_res * 0.5), color="#444444")

    # ---- Verdict panel ----
    ax = add_axes(title_h + dna_h + prot_h, verdict_h, total_units)
    ax.set_ylim(0, verdict_h)
    ax.set_xlim(0, 1)
    cursor = verdict_h - 0.2
    ax.text(0, cursor, category, fontsize=16, fontweight="bold", color=verdict_color, va="center")
    cursor -= 0.38
    changed_text = "NO" if np_ == mp_ else "YES"
    ax.text(0, cursor, f"Protein changed? {changed_text}", fontsize=11, fontweight="bold", va="center")
    cursor -= 0.35
    for line in effect_lines:
        ax.text(0, cursor, line, fontsize=10, va="center")
        cursor -= 0.25
    cursor -= 0.15
    ax.text(0, cursor, "Mutations found:", fontsize=10, fontweight="bold", va="center")
    cursor -= 0.25
    for line in mutation_lines:
        ax.text(0.01, cursor, "- " + line, fontsize=10, va="center", family="monospace")
        cursor -= 0.25

    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"\nImage saved to {path}")
    return True


def analyze_mutation(normal_input, mutated_input, start=None, image=None):
    """Analyze, print the report and (optionally) save a comparison picture."""
    results = analyze_sequences(normal_input, mutated_input, start=start)
    print_report(results)
    if image:
        save_comparison_image(results, image)
    return results


def main():
    parser = argparse.ArgumentParser(description="Compare normal and mutated DNA.")
    parser.add_argument("normal", nargs="?", help="normal DNA sequence")
    parser.add_argument("mutated", nargs="?", help="mutated DNA sequence")
    parser.add_argument("--start", type=int,
                        help="1-based position of the start codon in the normal sequence")
    parser.add_argument("--image", metavar="FILE.png",
                        help="also save a picture of the comparison (needs matplotlib)")
    args = parser.parse_args()

    try:
        if args.normal and args.mutated:
            analyze_mutation(args.normal, args.mutated, start=args.start, image=args.image)
        else:
            print("Enter DNA sequences using A, T, C and G.")
            normal_dna = input("Normal DNA sequence: ")
            mutated_dna = input("Mutated DNA sequence: ")
            start_text = input("Start codon position (blank = auto-detect): ").strip()
            analyze_mutation(normal_dna, mutated_dna,
                             start=int(start_text) if start_text else None,
                             image=args.image)
    except ValueError as error:
        print("Input error:", error)


if __name__ == "__main__":
    main()
