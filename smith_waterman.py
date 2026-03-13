import numpy as np

#fasta_file = open("data/X67616.1.fasta")

#string = fasta_file.read()

#fasta_file.close()

#rows = string.split("\n")

#nucleotide_seq = "".join(row.strip() for row in rows[1:])
#print(nucleotide_seq)

nucleotide_seq_input_1 = " " + "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG"#input("Enter nucleotide sequence: ")
nucleotide_seq_input_2 = " " + "CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG" #input("Enter nucleotide sequence: ")

#print(nucleotide_seq_input)

matrix = np.zeros((len(nucleotide_seq_input_2), len(nucleotide_seq_input_1)), dtype=int)

MATCH = 5
MISMATCH = -4
GAP = -1

score = 0
score_position = ()

for i in range(1, len(nucleotide_seq_input_2)):
    for j in range(1, len(nucleotide_seq_input_1)):

        if nucleotide_seq_input_2[i] == nucleotide_seq_input_1[j]:
            value_1 = max(matrix[i - 1][j - 1] + MATCH, 0)
        else:
            value_1 = max(matrix[i - 1][j - 1] + MISMATCH, 0)


        value_2 = max(matrix[i][j - 1] + GAP, matrix[i - 1][j] + GAP, 0)

        matrix[i][j] = max(value_1, value_2)

        if matrix[i][j] > score:
            score = matrix[i][j]
            score_position = (i, j)

print(matrix)
print(score)
print(score_position)

align_1 = ""
align_2 = ""

i, j = score_position

while matrix[i][j] != 0:
    current_score = matrix[i][j]

    diff = MATCH if nucleotide_seq_input_2[i] == nucleotide_seq_input_1[j] else MISMATCH

    if current_score == matrix[i - 1][j - 1] + diff:
        align_1 += nucleotide_seq_input_1[j]
        align_2 += nucleotide_seq_input_2[i]
        i, j = i - 1, j - 1
    elif current_score - GAP == matrix[i-1][j]:
        align_1 += nucleotide_seq_input_1[j]
        align_2 += '-'
        i, j = i - 1, j
    else:
        align_1 += '.'
        align_2 += nucleotide_seq_input_2[i]
        i, j = i, j - 1

print(align_1[::-1])
print(align_2[::-1])

