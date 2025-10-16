import os
from typing import Tuple
import pandas as pd
import BioAid as ba
import logging

#set up logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def annotate_genomic_context(
        df: pd.DataFrame,
        chrom_col: str = "Mapping",
        pos_col: str = "Reference Position",
        fasta_file_path: str = "",
        context_flank: int = 1,
        require_full_context: bool = True
    ) -> pd.DataFrame:
    """Annotates a DataFrame of variant calls with genomic sequence context.

    Uses BioAid.pullGenomicContext to retrieve flanking sequences for each variant.

    Args:
        df: DataFrame containing at least chromosome and position columns.
        chrom_col: Name of the chromosome column in df.
        pos_col: Name of the position column in df.
        fasta_file_path: Path to the reference FASTA file to query.
        context_flank: Number of bases to include on each side (e.g., 1 for trinucleotide).
        require_full_context: If True, rows where the returned sequence length
            is not equal to (2 * context_flank + 1) are dropped.

    Returns:
        A copy of the DataFrame with added columns: ctx_L, ctx_C, ctx_R, ctx_chr, seq, seq_len.
    """

    df = df.copy()
    # ensure integer positions
    df[pos_col] = df[pos_col].astype(int)
    list_of_positions = df[[chrom_col, pos_col]].values.tolist()
    result = ba.pullGenomicContext(list_of_positions, fasta_file_path, context_flank)
    # result expected: list of [left_context, center_base, right_context, chr, pos] or similar
    # code in this repo assigns into four columns; adapt if your BioAid returns different shape
    df[['ctx_L', 'ctx_C', 'ctx_R', 'ctx_chr']] = pd.DataFrame(result, index=df.index)
    df['seq'] = df['ctx_L'] + df['ctx_C'] + df['ctx_R']
    df['seq_len'] = df['seq'].apply(lambda x: len(x))
    if require_full_context:
        expected_len = context_flank * 2 + 1
        df = df[df['seq_len'] == expected_len].copy()
    return df

def annotate_a3a_context(annotated_df: pd.DataFrame) -> pd.DataFrame:
    """Annotates a DataFrame with A3A context.

    Args:
        annotated_df: DataFrame with a 'seq' column containing trinucleotide sequences.

    Returns:
        DataFrame with an added 'A3A_context' boolean column.
    """
    
    annotated_df['A3A_context'] = False #add column called A3A context

    # Ensure case-insensitivity by converting relevant columns to uppercase
    annotated_df['seq'] = annotated_df['seq'].str.upper()
    annotated_df['Reference'] = annotated_df['Reference'].str.upper()
    annotated_df['Allele'] = annotated_df['Allele'].str.upper()

    # Define A3A context trinucleotides
    a3a_context_trinucleotides_top = {"TCA", "TCT", "TCG"}
    a3a_context_trinucleotides_bottom = {"TGA", "AGA", "CGA"}

    # Condition for C>T or C>G mutations in top strand context
    cond1 = (
        (annotated_df['Reference'] == 'C') &
        (annotated_df['Allele'].isin(['T', 'G'])) &
        (annotated_df['seq'].isin(a3a_context_trinucleotides_top))
    )

    # Condition for G>A or G>T mutations in bottom strand context (reverse complement)
    cond2 = (
        (annotated_df['Reference'] == 'G') &
        (annotated_df['Allele'].isin(['A', 'T'])) &
        (annotated_df['seq'].isin(a3a_context_trinucleotides_bottom))
    )

    # Set A3A_context to True where either condition is met
    annotated_df.loc[cond1 | cond2, 'A3A_context'] = True

    return annotated_df

def remove_duplicate_mutations(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Removes duplicate mutations from a DataFrame.

    Args:
        df: DataFrame with mutation data.

    Returns:
        A tuple containing:
            - DataFrame with duplicate mutations removed.
            - DataFrame containing the removed duplicate rows.
    """

    subset_cols = ["Mapping","Reference Position","Reference","Allele","seq","A3A_context"]

    duplicate_rows = df[df.duplicated(subset=subset_cols, keep=False)]

    len_before = len(df)
    df = df.drop_duplicates(subset=subset_cols)
    len_after = len(df)
    print(f"Removed {len_before - len_after} duplicate mutations. {len_after} mutations remain:")

    #print number of mutations per origin file
    for origin_file in df['origin_file'].unique():
        count = df[df['origin_file'] == origin_file].shape[0]
        count_a3a = df[(df['origin_file'] == origin_file) & (df['A3A_context'] == True)].shape[0]
        print(f" - {origin_file}: {count} mutations, {count_a3a} ({count_a3a/count:.2%}) in A3A-like context")

    return df, duplicate_rows

if __name__ == "__main__":

    REFERENCE_GENOME_FASTA = "AM3422_de_novo.fasta"

    #find all csv files in mutation_calls directory
    csv_files = [f for f in os.listdir("mutation_calls") if f.endswith('.csv')]

    master_df = pd.DataFrame()  # Initialize an empty DataFrame to hold all annotated data

    for csv_file in csv_files:
        df = pd.read_csv(os.path.join("mutation_calls", csv_file))

        df = df[df['Type'] == 'SNV']

        # Following line removes the string " mapping" from the 'Mapping' column
        # This is may be necessary as the 'Mapping' column may contain extra 
        # descriptive text that doesn't match the reference FASTA naming.
        df['Mapping'] = df['Mapping'].str.replace(" mapping", "")

        annotated_df = annotate_genomic_context(df, fasta_file_path=REFERENCE_GENOME_FASTA)
        annotated_df = annotate_a3a_context(annotated_df)

        annotated_df.to_csv(os.path.join("annotated_calls", csv_file), index=False)

        annotated_df["origin_file"] = csv_file
        master_df = pd.concat([master_df, annotated_df], ignore_index=True)

    # Remove duplicate mutations across all files
    master_df, duplicate_rows = remove_duplicate_mutations(master_df)

    master_df.to_csv("annotated_calls/master_annotated_calls_deduped.csv", index=False)
    duplicate_rows.to_csv("annotated_calls/duplicate_rows.csv", index=False)
    