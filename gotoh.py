import numpy as np


# ednamat = {
#     'A': {'A': 10, 'C': -1, 'G': -3, 'T': -4},
#     'C': {'A': -1, 'C': 7, 'G': -5, 'T': -3},
#     'G': {'A': -3, 'C': -5, 'G': 9, 'T': 0},
#     'T': {'A': -4, 'C': -3, 'G': 0, 'T': 8}
# }
# cost = ednamat[n2[i]][n1[j]]
#fasta_file = open("data/X67616.1.fasta")

#string = fasta_file.read()

#fasta_file.close()

#rows = string.split("\n")

#nucleotide_seq = "".join(row.strip() for row in rows[1:])
#print(nucleotide_seq)

n1 = " " + "ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG"#input("Enter nucleotide sequence: ")
n2 = " " + "CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG" #input("Enter nucleotide sequence: ")

n = len(n1)
m = len(n2)

#print(nucleotide_seq_input)

u = -1
v = -5
match = 5
mismatch = -4

#prefer = 'D' # 'D' - diagonal, 'P' - to prefer gap in b seq, 'Q' - to prefer gap in a seq

# Inicializacia D matrix. prvy riadok a prvy stlpec
d_matrix = np.zeros((len(n2), len(n1)), dtype=float)
#d_matrix[0, 1] = v
# for j in range(1, len(n1)):
#     d_matrix[0, j] = v + u + (j-1) * u
#
# for i in range(1, len(n2)):
#     d_matrix[i, 0] = v  + u + (i-1) * u

# Inicializacia P matrix. Prvy riadok su nekonecno lebo je nam zbytocne zacinat seq empty a sparovat to s medzerou. teda dlhsia sekvencia bude mat pri inicializacii nekonecno, prvy stlpec nuly vratane empty empty
p_matrix = np.full((len(n2), len(n1)), -np.inf, dtype=float)
# p_matrix[:, 0] = 0
# p_matrix[0, 1:] = -np.inf

# Inicializacia Q matrix. To iste ale pre stlpec. prvy riadok nuly
q_matrix = np.full((len(n2), len(n1)), -np.inf, dtype=float)
# q_matrix[0, :] = 0
# q_matrix[1:, 0] = -np.inf

# Ideme od i, j = 1, 1 a pozerame sa do kazdej matice na to miesto
for i in range(1, len(n2)):
    for j in range(1, len(n1)):
        # Zaciname v P = pocitame to min(Di-1, j + u + v al P i-j, j + u
        p_matrix[i,j] = max((p_matrix[i - 1, j] + u), (d_matrix[i-1, j] + v + u))
        p_matrix[i, j] = max(0, p_matrix[i, j])

        # Teraz to iste pre Q
        q_matrix[i,j] = max((q_matrix[i, j - 1] + u), (d_matrix[i, j - 1] + v + u))
        q_matrix[i, j] = max(0, q_matrix[i, j])

        # Nakoniec pozrieme D a najdeme ten min
        cost = match if n2[i] == n1[j] else mismatch
        d_matrix[i][j] = max(0, p_matrix[i, j], q_matrix[i, j], d_matrix[i - 1, j - 1] + cost)

# Until the last cell


# --- FORMÁTOVANÝ VÝPIS ---
width = 5  # Šírka stĺpca pre pekné zarovnanie

# P matrix display
# 1. Horný riadok s písmenami sekvencie 1
header = " " * width  # Prázdne miesto pre bočnú lištu
for char in n1:
    header += f"{char:>{width}}"
print(header)

# 2. Riadky matice s bočným písmenom zo sekvencie 2
for i, row in enumerate(p_matrix):
    # Začneme písmenom zo sekvencie 2
    row_str = f"{n2[i]:>{width}}"
    # Pridáme čísla z matice
    for val in row:
        row_str += f"{val:>{width}}"
    print(row_str)

# Q matrix display
# 1. Horný riadok s písmenami sekvencie 1
header = " " * width  # Prázdne miesto pre bočnú lištu
for char in n1:
    header += f"{char:>{width}}"
print(header)

# 2. Riadky matice s bočným písmenom zo sekvencie 2
for i, row in enumerate(q_matrix):
    # Začneme písmenom zo sekvencie 2
    row_str = f"{n2[i]:>{width}}"
    # Pridáme čísla z matice
    for val in row:
        row_str += f"{val:>{width}}"
    print(row_str)

# D matric display
# 1. Horný riadok s písmenami sekvencie 1
header = " " * width  # Prázdne miesto pre bočnú lištu
for char in n1:
    header += f"{char:>{width}}"
print(header)

# 2. Riadky matice s bočným písmenom zo sekvencie 2
for i, row in enumerate(d_matrix):
    # Začneme písmenom zo sekvencie 2
    row_str = f"{n2[i]:>{width}}"
    # Pridáme čísla z matice
    for val in row:
        row_str += f"{val:>{width}}"
    print(row_str)


start = np.unravel_index(np.argmax(d_matrix), d_matrix.shape)

i, j = start
print(f"Score: {d_matrix[i, j]}")
print(start)
# traceback
a1 = []
a2 = []

state = 'D'

while d_matrix[i,j] != 0:

    if state == 'D':

        cost = match if n2[i] == n1[j] else mismatch

        if d_matrix[i,j] == d_matrix[i-1,j-1] + cost:
            a1.append(n1[j])
            a2.append(n2[i])
            i -= 1
            j -= 1

        elif d_matrix[i,j] == p_matrix[i,j]:
            state = 'P'

        else:
            state = 'Q'

    elif state == 'P':
        a1.append('-')
        a2.append(n2[i])

        if p_matrix[i,j] == d_matrix[i-1,j] + v + u:
            state = 'D'

        i -= 1

    else:  # Q
        a1.append(n1[j])
        a2.append('-')

        if q_matrix[i,j] == d_matrix[i,j-1] + v + u:
            state = 'D'

        j -= 1
a1.reverse()
a2.reverse()
print(a1)
print(a2)


#
#
# matrix = np.zeros((len(nucleotide_seq_input_2), len(nucleotide_seq_input_1)), dtype=int)
#
# MATCH = 5
# MISMATCH = -4
# GAP = -1
#
# score = 0
# score_position = ()
#
# for i in range(1, len(nucleotide_seq_input_2)):
#     for j in range(1, len(nucleotide_seq_input_1)):
#
#         if nucleotide_seq_input_2[i] == nucleotide_seq_input_1[j]:
#             value_1 = max(matrix[i - 1][j - 1] + MATCH, 0)
#         else:
#             value_1 = max(matrix[i - 1][j - 1] + MISMATCH, 0)
#
#
#         value_2 = max(matrix[i][j - 1] + GAP, matrix[i - 1][j] + GAP, 0)
#
#         matrix[i][j] = max(value_1, value_2)
#
#         if matrix[i][j] > score:
#             score = matrix[i][j]
#             score_position = (i, j)
#
# print(matrix)
# print(score)
# print(score_position)
#
# align_1 = ""
# align_2 = ""
#
# i, j = score_position
#
# while matrix[i][j] != 0:
#     current_score = matrix[i][j]
#
#     diff = MATCH if nucleotide_seq_input_2[i] == nucleotide_seq_input_1[j] else MISMATCH
#
#     if current_score == matrix[i - 1][j - 1] + diff:
#         align_1 += nucleotide_seq_input_1[j]
#         align_2 += nucleotide_seq_input_2[i]
#         i, j = i - 1, j - 1
#     elif current_score - GAP == matrix[i-1][j]:
#         align_1 += nucleotide_seq_input_1[j]
#         align_2 += '-'
#         i, j = i - 1, j
#     else:
#         align_1 += '.'
#         align_2 += nucleotide_seq_input_2[i]
#         i, j = i, j - 1
#
# print(align_1[::-1])
# print(align_2[::-1])
#
