# Mutation Annotation Pipeline

This project provides a pipeline for annotating mutation calls with genomic context and identifying specific mutation signatures, such as those associated with APOBEC3A (A3A).

## Description

The main script, `annotate_context.py`, performs the following steps:
1. Reads mutation call files (in CSV format) from the `mutation_calls` directory.
2. Filters for Single Nucleotide Variants (SNVs).
3. Annotates each SNV with its genomic context (trinucleotide sequence) using a reference FASTA file.
4. Identifies mutations occurring in an A3A-like context.
5. Writes the annotated calls to the `annotated_calls` directory.
6. Consolidates all annotated calls into a master file, removing duplicates.

## Installation

1. **Clone the repository:**
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. **Install required Python packages:**
   ```bash
   pip install pandas BioAid
   ```

## Usage

1. **Place your reference genome** in FASTA format in the root directory of the project (e.g., `AM3422_de_novo.fasta`).
2. **Place your mutation call CSV files** in the `mutation_calls` directory.
3. **Run the script:**
   ```bash
   python annotate_context.py
   ```
4. **Find the results:**
   - Annotated files for each input CSV will be in the `annotated_calls` directory.
   - A consolidated, deduplicated master file will be at `annotated_calls/master_annotated_calls_deduped.csv`.

## Project Structure

```
.
├── AM3422_de_novo.fasta        # Reference genome
├── annotate_context.py         # Main analysis script
├── annotated_calls/            # Output directory for annotated files
│   ├── ...
│   ├── duplicate_rows.csv      # Duplicate mutation calls
│   └── master_annotated_calls_deduped.csv # All mutations without duplicates
├── mutation_calls/             # Input directory for mutation CSVs
│   └── ...
├── .gitignore                  # Git ignore rules
└── README.md                   # This file
```
