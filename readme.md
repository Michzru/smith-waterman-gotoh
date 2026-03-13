# Smith-Waterman Sequence Alignment (Gotoh Optimization)

Implementation of the **Smith–Waterman local alignment algorithm** with **affine gap penalties using Gotoh optimization**.
The tool supports both **DNA and protein sequence alignment** and allows choosing between multiple substitution matrices.

## Features

* Local sequence alignment using Smith-Waterman
* Affine gap penalties (gap open + gap extend)
* Gotoh optimization with three matrices:

  * **D** – main score matrix
  * **P** – vertical gaps
  * **Q** – horizontal gaps
* Support for both **DNA and protein sequences**
* Multiple substitution matrices:

  * `EDNAMAT`
  * `EDNAFULL`
  * `BLOSUM62`
  * `PAM250`
* Outputs:

  * alignment statistics
  * aligned sequences
  * scoring matrices (D, P, Q)

Results are saved automatically into the `outputs/` directory.

## Project Structure

```
.
├── data/          # Example FASTA files for testing
├── outputs/       # Generated alignment results
├── gotoh.py       # Main script
```

The `data/` folder contains several **ready-to-use FASTA sequences** for quick testing.

Example files include:

* `dna_homo_sapiens_insulin.fasta`
* `dna_pan_troglodytes_insulin.fasta`
* `protein_homo_sapiens_hemoglobin.fasta`
* `protein_rat_hemoglobin.fasta`
* `test_dna_*.fasta`
* `test_protein_*.fasta`

## Usage

Run the script from the command line:

```
python gotoh.py -t TYPE -m MATRIX -go GAP_OPEN -ge GAP_EXTEND seq1.fasta seq2.fasta
```

### Arguments

| Argument            | Description                                                     |
| ------------------- |-----------------------------------------------------------------|
| `-t, --type`        | Sequence type: `dna` or `prt`                                   |
| `-m, --matrix`      | Substitution matrix: `EDNAMAT`, `EDNAFULL`, `BLOSUM62`, `PAM250`|
| `-go, --gap_open`   | Gap opening penalty                                             |
| `-ge, --gap_extend` | Gap extension penalty                                           |
| `seq1`              | Path to first FASTA sequence                                    |
| `seq2`              | Path to second FASTA sequence                                   |

### Example (DNA)

```
python gotoh.py -t dna -m EDNAMAT -go 10 -ge 1 data/test_dna_1.fasta data/test_dna_2.fasta
```

### Example (Protein)

```
python gotoh.py -t prt -m BLOSUM62 -go 10 -ge 1 data/test_protein_1.fasta data/test_protein_2.fasta
```

## Output

After execution, results will be stored in:

```
outputs/<sequence1>_<sequence2>_alignment.txt
```

The output file contains:

* alignment score
* identity and similarity statistics
* gap statistics
* formatted alignment
* full scoring matrices

## Requirements

* Python 3
* NumPy

Install dependencies:

```
pip install numpy
```

## Notes

The implementation validates sequence alphabets to ensure compatibility with the selected sequence type and substitution matrix.
