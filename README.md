# APOBEC Annotator

This project provides a pipeline for annotating mutation calls with genomic context and identifying specific mutation signatures, such as those associated with APOBEC3A (A3A).

## Description

The main script, `annotate_context.py`, performs the following steps:
1. Reads mutation call files (in CSV format) from the `mutation_calls` directory.
2. Filters for Single Nucleotide Variants (SNVs).
3. Annotates each SNV with its genomic context (trinucleotide sequence) using a reference FASTA file.
4. Identifies mutations occurring in an A3A-like context:
   - Mutation is a C>T or C>G mutation
   - Mutation satisfies one of the following:
      - Mutation occurs in a tCa, tCt, or tCg trinucleotide context (positive strand)
      - Mutation occurs in a tGa, aGa, or cGa trinucleotide context (negative strand)
5. Writes the annotated calls to the `annotated_calls` directory.
6. Consolidates all annotated calls into a master file, removing duplicates.
7. Saves all duplicate mutation calls to a separate file.

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

Script is set up to work with outputs of CLC Genomics Workbench. To use with other variant callers, 
ensure your CSVs have these columns: 
 - `Mapping` - Chromosome or contig name. Ensure naming matches reference FASTA.
 - `Reference Position` - Position of the variant on the reference
 - `Reference` - Reference Allele
 - `Allele` - Variant Allele
 - `Type` - Mutation type (e.g., SNV)

 You can adapt the script for other callers by ensuring these columns are present, or you can modify the script accordingly, for example by changing the default column names in the `annotate_genomic_context()` function.

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
