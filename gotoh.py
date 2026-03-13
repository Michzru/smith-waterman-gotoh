import argparse
import os
import numpy as np
from pathlib import Path

# EDNAFULL matrix for DNA match/mismatch
EDNAFULL = {
    "A": {"A": 5, "C": -4, "G": -4, "T": -4, "R": 1, "Y": -4, "N": -1},
    "C": {"A": -4, "C": 5, "G": -4, "T": -4, "R": -4, "Y": 1, "N": -1},
    "G": {"A": -4, "C": -4, "G": 5, "T": -4, "R": 1, "Y": -4, "N": -1},
    "T": {"A": -4, "C": -4, "G": -4, "T": 5, "R": -4, "Y": 1, "N": -1},
    "R": {"A": 1, "C": -4, "G": 1, "T": -4, "R": 1, "Y": -4, "N": -1},
    "N": {"A": -1, "C": -1, "G": -1, "T": -1, "R": -1, "Y": -1, "N": -1}
}

# OR

# EDNAMAT matrix for DNA match/mismatch
EDNAMAT = {
"A":{"A":5,"C":-4,"G":-4,"T":-4},
"C":{"A":-4,"C":5,"G":-4,"T":-4},
"G":{"A":-4,"C":-4,"G":5,"T":-4},
"T":{"A":-4,"C":-4,"G":-4,"T":5}
}

# or

# BLOSUM62 SUBSTITUTION MATRIX for protein sequences
BLOSUM62 = {
"A":{"A":4,"R":-1,"N":-2,"D":-2,"C":0,"Q":-1,"E":-1,"G":0,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":0,"W":-3,"Y":-2,"V":0},
"R":{"A":-1,"R":5,"N":0,"D":-2,"C":-3,"Q":1,"E":0,"G":-2,"H":0,"I":-3,"L":-2,"K":2,"M":-1,"F":-3,"P":-2,"S":-1,"T":-1,"W":-3,"Y":-2,"V":-3},
"N":{"A":-2,"R":0,"N":6,"D":1,"C":-3,"Q":0,"E":0,"G":0,"H":1,"I":-3,"L":-3,"K":0,"M":-2,"F":-3,"P":-2,"S":1,"T":0,"W":-4,"Y":-2,"V":-3},
"D":{"A":-2,"R":-2,"N":1,"D":6,"C":-3,"Q":0,"E":2,"G":-1,"H":-1,"I":-3,"L":-4,"K":-1,"M":-3,"F":-3,"P":-1,"S":0,"T":-1,"W":-4,"Y":-3,"V":-3},
"C":{"A":0,"R":-3,"N":-3,"D":-3,"C":9,"Q":-3,"E":-4,"G":-3,"H":-3,"I":-1,"L":-1,"K":-3,"M":-1,"F":-2,"P":-3,"S":-1,"T":-1,"W":-2,"Y":-2,"V":-1},
"Q":{"A":-1,"R":1,"N":0,"D":0,"C":-3,"Q":5,"E":2,"G":-2,"H":0,"I":-3,"L":-2,"K":1,"M":0,"F":-3,"P":-1,"S":0,"T":-1,"W":-2,"Y":-1,"V":-2},
"E":{"A":-1,"R":0,"N":0,"D":2,"C":-4,"Q":2,"E":5,"G":-2,"H":0,"I":-3,"L":-3,"K":1,"M":-2,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
"G":{"A":0,"R":-2,"N":0,"D":-1,"C":-3,"Q":-2,"E":-2,"G":6,"H":-2,"I":-4,"L":-4,"K":-2,"M":-3,"F":-3,"P":-2,"S":0,"T":-2,"W":-2,"Y":-3,"V":-3},
"H":{"A":-2,"R":0,"N":1,"D":-1,"C":-3,"Q":0,"E":0,"G":-2,"H":8,"I":-3,"L":-3,"K":-1,"M":-2,"F":-1,"P":-2,"S":-1,"T":-2,"W":-2,"Y":2,"V":-3},
"I":{"A":-1,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-3,"E":-3,"G":-4,"H":-3,"I":4,"L":2,"K":-3,"M":1,"F":0,"P":-3,"S":-2,"T":-1,"W":-3,"Y":-1,"V":3},
"L":{"A":-1,"R":-2,"N":-3,"D":-4,"C":-1,"Q":-2,"E":-3,"G":-4,"H":-3,"I":2,"L":4,"K":-2,"M":2,"F":0,"P":-3,"S":-2,"T":-1,"W":-2,"Y":-1,"V":1},
"K":{"A":-1,"R":2,"N":0,"D":-1,"C":-3,"Q":1,"E":1,"G":-2,"H":-1,"I":-3,"L":-2,"K":5,"M":-1,"F":-3,"P":-1,"S":0,"T":-1,"W":-3,"Y":-2,"V":-2},
"M":{"A":-1,"R":-1,"N":-2,"D":-3,"C":-1,"Q":0,"E":-2,"G":-3,"H":-2,"I":1,"L":2,"K":-1,"M":5,"F":0,"P":-2,"S":-1,"T":-1,"W":-1,"Y":-1,"V":1},
"F":{"A":-2,"R":-3,"N":-3,"D":-3,"C":-2,"Q":-3,"E":-3,"G":-3,"H":-1,"I":0,"L":0,"K":-3,"M":0,"F":6,"P":-4,"S":-2,"T":-2,"W":1,"Y":3,"V":-1},
"P":{"A":-1,"R":-2,"N":-2,"D":-1,"C":-3,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-3,"L":-3,"K":-1,"M":-2,"F":-4,"P":7,"S":-1,"T":-1,"W":-4,"Y":-3,"V":-2},
"S":{"A":1,"R":-1,"N":1,"D":0,"C":-1,"Q":0,"E":0,"G":0,"H":-1,"I":-2,"L":-2,"K":0,"M":-1,"F":-2,"P":-1,"S":4,"T":1,"W":-3,"Y":-2,"V":-2},
"T":{"A":0,"R":-1,"N":0,"D":-1,"C":-1,"Q":-1,"E":-1,"G":-2,"H":-2,"I":-1,"L":-1,"K":-1,"M":-1,"F":-2,"P":-1,"S":1,"T":5,"W":-2,"Y":-2,"V":0},
"W":{"A":-3,"R":-3,"N":-4,"D":-4,"C":-2,"Q":-2,"E":-3,"G":-2,"H":-2,"I":-3,"L":-2,"K":-3,"M":-1,"F":1,"P":-4,"S":-3,"T":-2,"W":11,"Y":2,"V":-3},
"Y":{"A":-2,"R":-2,"N":-2,"D":-3,"C":-2,"Q":-1,"E":-2,"G":-3,"H":2,"I":-1,"L":-1,"K":-2,"M":-1,"F":3,"P":-3,"S":-2,"T":-2,"W":2,"Y":7,"V":-1},
"V":{"A":0,"R":-3,"N":-3,"D":-3,"C":-1,"Q":-2,"E":-2,"G":-3,"H":-3,"I":3,"L":1,"K":-2,"M":1,"F":-1,"P":-2,"S":-2,"T":0,"W":-3,"Y":-1,"V":4}
}

# PAM250 substitution matrix for protein sequences
PAM250 = {
"A":{"A":2,"R":-2,"N":0,"D":0,"C":-2,"Q":0,"E":0,"G":1,"H":-1,"I":-1,"L":-2,"K":-1,"M":-1,"F":-3,"P":1,"S":1,"T":1,"W":-6,"Y":-3,"V":0},
"R":{"A":-2,"R":6,"N":0,"D":-1,"C":-4,"Q":1,"E":-1,"G":-3,"H":2,"I":-2,"L":-3,"K":3,"M":0,"F":-4,"P":0,"S":0,"T":-1,"W":2,"Y":-4,"V":-2},
"N":{"A":0,"R":0,"N":2,"D":2,"C":-4,"Q":1,"E":1,"G":0,"H":2,"I":-2,"L":-3,"K":1,"M":-2,"F":-3,"P":0,"S":1,"T":0,"W":-4,"Y":-2,"V":-2},
"D":{"A":0,"R":-1,"N":2,"D":4,"C":-5,"Q":2,"E":3,"G":1,"H":1,"I":-2,"L":-4,"K":0,"M":-3,"F":-6,"P":-1,"S":0,"T":0,"W":-7,"Y":-4,"V":-2},
"C":{"A":-2,"R":-4,"N":-4,"D":-5,"C":12,"Q":-5,"E":-5,"G":-3,"H":-3,"I":-2,"L":-6,"K":-5,"M":-5,"F":-4,"P":-3,"S":0,"T":-2,"W":-8,"Y":0,"V":-2},
"Q":{"A":0,"R":1,"N":1,"D":2,"C":-5,"Q":4,"E":2,"G":-1,"H":3,"I":-2,"L":-2,"K":1,"M":-1,"F":-5,"P":0,"S":-1,"T":-1,"W":-5,"Y":-4,"V":-2},
"E":{"A":0,"R":-1,"N":1,"D":3,"C":-5,"Q":2,"E":4,"G":0,"H":1,"I":-2,"L":-3,"K":0,"M":-2,"F":-5,"P":-1,"S":0,"T":0,"W":-7,"Y":-4,"V":-2},
"G":{"A":1,"R":-3,"N":0,"D":1,"C":-3,"Q":-1,"E":0,"G":5,"H":-2,"I":-3,"L":-4,"K":-2,"M":-3,"F":-5,"P":0,"S":1,"T":0,"W":-7,"Y":-5,"V":-1},
"H":{"A":-1,"R":2,"N":2,"D":1,"C":-3,"Q":3,"E":1,"G":-2,"H":6,"I":-2,"L":-2,"K":0,"M":-2,"F":-2,"P":0,"S":-1,"T":-1,"W":-3,"Y":0,"V":-2},
"I":{"A":-1,"R":-2,"N":-2,"D":-2,"C":-2,"Q":-2,"E":-2,"G":-3,"H":-2,"I":5,"L":2,"K":-2,"M":2,"F":1,"P":-2,"S":-1,"T":0,"W":-5,"Y":-1,"V":4},
"L":{"A":-2,"R":-3,"N":-3,"D":-4,"C":-6,"Q":-2,"E":-3,"G":-4,"H":-2,"I":2,"L":6,"K":-3,"M":4,"F":2,"P":-3,"S":-3,"T":-2,"W":-2,"Y":-1,"V":2},
"K":{"A":-1,"R":3,"N":1,"D":0,"C":-5,"Q":1,"E":0,"G":-2,"H":0,"I":-2,"L":-3,"K":5,"M":0,"F":-5,"P":-1,"S":0,"T":0,"W":-3,"Y":-4,"V":-2},
"M":{"A":-1,"R":0,"N":-2,"D":-3,"C":-5,"Q":-1,"E":-2,"G":-3,"H":-2,"I":2,"L":4,"K":0,"M":6,"F":0,"P":-2,"S":-2,"T":-1,"W":-4,"Y":-2,"V":2},
"F":{"A":-3,"R":-4,"N":-3,"D":-6,"C":-4,"Q":-5,"E":-5,"G":-5,"H":-2,"I":1,"L":2,"K":-5,"M":0,"F":9,"P":-5,"S":-3,"T":-3,"W":0,"Y":7,"V":-1},
"P":{"A":1,"R":0,"N":0,"D":-1,"C":-3,"Q":0,"E":-1,"G":0,"H":0,"I":-2,"L":-3,"K":-1,"M":-2,"F":-5,"P":6,"S":1,"T":0,"W":-6,"Y":-5,"V":-1},
"S":{"A":1,"R":0,"N":1,"D":0,"C":0,"Q":-1,"E":0,"G":1,"H":-1,"I":-1,"L":-3,"K":0,"M":-2,"F":-3,"P":1,"S":2,"T":1,"W":-2,"Y":-3,"V":-1},
"T":{"A":1,"R":-1,"N":0,"D":0,"C":-2,"Q":-1,"E":0,"G":0,"H":-1,"I":0,"L":-2,"K":0,"M":-1,"F":-3,"P":0,"S":1,"T":3,"W":-5,"Y":-3,"V":0},
"W":{"A":-6,"R":2,"N":-4,"D":-7,"C":-8,"Q":-5,"E":-7,"G":-7,"H":-3,"I":-5,"L":-2,"K":-3,"M":-4,"F":0,"P":-6,"S":-2,"T":-5,"W":17,"Y":0,"V":-6},
"Y":{"A":-3,"R":-4,"N":-2,"D":-4,"C":0,"Q":-4,"E":-4,"G":-5,"H":0,"I":-1,"L":-1,"K":-4,"M":-2,"F":7,"P":-5,"S":-3,"T":-3,"W":0,"Y":10,"V":-2},
"V":{"A":0,"R":-2,"N":-2,"D":-2,"C":-2,"Q":-2,"E":-2,"G":-1,"H":-2,"I":4,"L":2,"K":-2,"M":2,"F":-1,"P":-1,"S":-1,"T":0,"W":-6,"Y":-2,"V":4}
}

parser = argparse.ArgumentParser(
    description="Smith-Waterman local alignment for DNA or protein sequences"
)

parser.add_argument(
    "-t", "--type",
    choices=["dna", "prt"],
    required=True,
    help="Type of sequence (dna or prt)"
)

parser.add_argument(
    "-m", "--matrix",
    choices=["EDNAMAT", "EDNAFULL", "BLOSUM62", "PMA250"],
    required=True,
    help="Substitution matrix (EDNAMAT, EDNAFULL, BLOSUM62, PMA250)"
)

parser.add_argument(
    "-go", "--gap_open",
    type=float,
    required=True,
    help="Value of gap open penalization"
)

parser.add_argument(
    "-ge", "--gap_extend",
    type=float,
    required=True,
    help="Value of gap extend penalization"
)

parser.add_argument(
    "seq1",
    help="Path to first FASTA sequence"
)

parser.add_argument(
    "seq2",
    help="Path to second FASTA sequence"
)

args = parser.parse_args()

sequence_type = args.type
matrix_type = args.matrix
sequence_1_path = args.seq1
sequence_2_path = args.seq2
name1 = Path(sequence_1_path).stem
name2 = Path(sequence_2_path).stem
output_filename = f"outputs/{name1}_{name2}_alignment.txt"
file_path = Path(output_filename)
file_path.parent.mkdir(parents=True, exist_ok=True)
u = -abs(args.gap_extend)
v = -abs(args.gap_open)

matrix_cost = None

# Select type of matrix
if sequence_type == 'dna':
    if matrix_type not in ['EDNAMAT', 'EDNAFULL']:
        exit("For DNA you must use EDNAMAT or EDNAFULL")
    if matrix_type == 'EDNAMAT':
        matrix_cost = EDNAMAT
    else:
        matrix_cost = EDNAFULL

elif sequence_type == 'prt':
    if matrix_type not in ['BLOSUM62', 'PAM250']:
        exit("For PROTEIN you must use BLOSUM62 or PAM250")
    if matrix_type == 'BLOSUM62':
        matrix_cost = BLOSUM62
    else:
        matrix_cost = PAM250

# Check paths
if not os.path.exists(sequence_1_path):
    exit(f"File not found: {sequence_1_path}")
elif not os.path.exists(sequence_2_path):
    exit(f"File not found: {sequence_2_path}")

# FASTA LOADER
def load_fasta(path):
    sequence = ""

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                continue

            sequence += line.upper()

    return sequence

# Load seq from fasta files (remove header and white spaces)
seq1 = load_fasta(sequence_1_path)
seq2 = load_fasta(sequence_2_path)

# Validate seq using predefined alphabet
DNA_ALPHABET = set("ACGTRYN")
PROTEIN_ALPHABET = set("ARNDCQEGHILKMFPSTWYV")

def is_dna(seq):
    return set(seq).issubset(DNA_ALPHABET)

def is_protein(seq):
    return set(seq).issubset(PROTEIN_ALPHABET)

# making sure the sequence correspond to at least one alphabet
if sequence_type == "dna":
    if not is_dna(seq1) or not is_dna(seq2):
        exit("One of the sequences is not valid DNA (allowed: A,C,G,T, R, Y, N)")
elif sequence_type == "prt":
    if not is_protein(seq1) or not is_protein(seq2):
        exit("One of the sequences is not valid protein sequence")

print("Parameters:")
print(f"\tSequence type: {sequence_type}")
print(f"\tMatrix type: {matrix_type}")
print(f"\tSequence 1 path: {sequence_1_path}")
print(f"\tSequence 2 path: {sequence_2_path}")

print("\nRunning Smith Waterman using Gotoh optimization...")

# space to align sequences with empty spaces in matrix
n1 = " " + seq1
n2 = " " + seq2

# Initialization of matrices for affined penalties

# D matrix (main) - storing highest possible score for i, j
d_matrix = np.zeros((len(n2), len(n1)), dtype=float)

# P matrix (vertical gaps) - storing score for padding ending with gap in seq 1 (indel)
# initialized with -inf to eliminate illegal paths on the start as empty + -
p_matrix = np.full((len(n2), len(n1)), -np.inf, dtype=float)

# Q matrix (horizontal) - storing score for padding ending with gap in seq 2 (indel)
# initialized with -inf to eliminate illegal paths on the start as empty + -
q_matrix = np.full((len(n2), len(n1)), -np.inf, dtype=float)


# From i, j = 1, 1 we look into each matrix until the last cell
# Using dynamic programming to calculate matrices
for i in range(1, len(n2)):
    for j in range(1, len(n1)):
        # Zaciname v P = pocitame to min(Di-1, j + u + v al P i-j, j + u

        # calculate P, new vertical gap, or extending gap, depends on where we came from
        p_matrix[i,j] = max((p_matrix[i - 1, j] + u), (d_matrix[i-1, j] + v))

        # calculate Q, new horizontal gap, or extending gap, depends on where we came from
        q_matrix[i,j] = max((q_matrix[i, j - 1] + u), (d_matrix[i, j - 1] + v))

        # Calculate D score, selecting either new gap, ending gap, or indel
        cost = matrix_cost[n2[i]][n1[j]]
        d_matrix[i][j] = max(0, p_matrix[i, j], q_matrix[i, j], d_matrix[i - 1, j - 1] + cost)

# Find the maximum score in the matrix to extract indices
start = np.unravel_index(np.argmax(d_matrix), d_matrix.shape)

i, j = start
# extract final score
score = d_matrix[i, j]

# TRACEBACK
a1 = []
a2 = []

# statistics
gaps = 0
state = 'D' # the start state, we are in D matrix

while d_matrix[i,j] != 0:

    # State D (match/mismatch or translation into gap)
    if state == 'D':

        # match mismatch score
        cost = matrix_cost[n2[i]][n1[j]]

        # Diagonal move (match or mismatch)
        if d_matrix[i,j] == d_matrix[i-1,j-1] + cost:
            a1.append(n1[j])
            a2.append(n2[i])
            i -= 1
            j -= 1

        # Vertical move into P
        elif d_matrix[i,j] == p_matrix[i,j]:
            state = 'P'

        # Horizontal move into Q
        else:
            state = 'Q'

    # State P (vertical gap in seq 1)
    elif state == 'P':
        # adding gap into seq 1
        a1.append('-')
        a2.append(n2[i])

        gaps += 1

        # Does gap end?
        if p_matrix[i,j] == d_matrix[i-1,j] + v:
            state = 'D'

        # move only in seq 2
        i -= 1

    # State Q (horizontal gap in seq 2)
    else:
        a1.append(n1[j])
        a2.append('-')

        gaps += 1

        # Does gap end?
        if q_matrix[i,j] == d_matrix[i,j-1] + v:
            state = 'D'

        # Move only in seq 1
        j -= 1

# The seqes were collected from back to start
a1.reverse()
a2.reverse()

sequence_1 = ""
sequence_2 = ""
sequence_1 += "".join(a1)
sequence_2 += "".join(a2)

# Some statistics
is_dna = True
print(f"\n{30*"#"} {"DNA" if sequence_type == "dna" else "PROTEIN"} STATISTICS {30*"#"}")
print("#" + 78*"-" + "#")

# Matrix Used: Match Mismatch
print(f"# Used matrix: {matrix_type}")

# Gap Penalties
print(f"# Gap penalty open: {v}")
print(f"# Gap penalty extend: {u}")

# Length
print("#" + 78*"-" + "#")
print(f"# Length: {len(sequence_1)}")

# Identity
print("#" + 78*"-" + "#")
count = 0
count = sum(1 for i in range(len(a1)) if a1[i] == a2[i])
print(f"# Identity: {count}/{len(a1)} ({count/len(a1)*100:.2f}%)")

# Similarity
similarity = 0
if args.type == 'prt':
    for x, y in zip(a1, a2):
        if x != '-' and y != '-':
            if x == y:
                similarity += 1
            elif matrix_cost[x][y] > 0:
                similarity += 1
    print(f"# Similarity: {similarity}/{len(a1)} ({similarity / len(a1) * 100:.2f}%)")

else:
    similarity = count
    print(f"# Similarity: {similarity}/{len(a1)} ({similarity / len(a1) * 100:.2f}%)")

# Gaps
print(f"# Number of gaps: {gaps} ({gaps/len(a1)*100:.2f}%)")

# Score
print("#" + 78*"-" + "#")
print(f"# Final Score: {score}")

# Alignment
limit = 50
match_line = "".join("|" if a1[i] == a2[i] else " " for i in range(len(a1)))
print("#" + 78*"-" + "#")
print(f"# Sequence 1 (first 50): {sequence_1[:limit]}...")
print(f"#                        {match_line[:limit]}...")
print(f"# Sequence 2 (first 50): {sequence_2[:limit]}...")


# To format matrix when storing in txt file
def format_matrix(matrix, seq1, seq2, title):
    width = 8
    output = f"\n--- {title} MATRIX ---\n"

    header = " " * width
    for char in seq1:
        header += f"{char:>{width}}"
    output += header + "\n"

    for i, row in enumerate(matrix):
        row_str = f"{seq2[i]:>{width}}"
        for val in row:
            if val == -np.inf:
                row_str += f"{'-inf':>{width}}"
            else:
                row_str += f"{val:>{width}.1f}"
        output += row_str + "\n"
    return output


with open(output_filename, "w", encoding="utf-8") as out_file:
    # Statistics
    out_file.write(f"Alignment of {name1} and {name2}\n")
    out_file.write(f"Matrix used: {matrix_type}\n")
    out_file.write(f"Gap penalty open: {v}\n")
    out_file.write(f"Gap penalty extend: {u}\n")
    out_file.write(f"Length: {len(sequence_1)}\n")
    out_file.write(f"Identity: {count}/{len(a1)} ({count / len(a1) * 100:.2f}%)\n")
    out_file.write(f"Similarity: {similarity}/{len(a1)} ({similarity / len(a1) * 100:.2f}%)\n")
    out_file.write(f"Score: {score}\n")
    out_file.write("-" * 60 + "\n")

    # Alignement
    out_file.write(f"Seq1: {sequence_1}\n")
    out_file.write(f"      {match_line}\n")
    out_file.write(f"Seq2: {sequence_2}\n")
    out_file.write("-" * 60 + "\n")

    # Matrices
    out_file.write(format_matrix(d_matrix, n1, n2, "D (Main Score)"))
    out_file.write(format_matrix(p_matrix, n1, n2, "P (Vertical Gaps)"))
    out_file.write(format_matrix(q_matrix, n1, n2, "Q (Horizontal Gaps)"))
print("#" + 78*"-" + "#")

print(f"\nDone! Everything can be found in: {output_filename}")