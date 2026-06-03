from ViennaRNA import RNA
from Bio import SeqIO
import pandas as pd
import joblib
from tqdm import tqdm

def _compute_mfe(window_seq):
    """
    Helper function to compute MFE for a single window sequence.

    Parameters:
        window_seq (str): RNA sequence

    Returns:
        float: MFE value
    """
    fc = RNA.fold_compound(window_seq)
    (_, mfe) = fc.mfe()
    return mfe



def extract_windows_with_mfe_centered(fasta_file, window_sizes=[30, 60, 90], step=10, n_jobs=-1):
    """
    Extract sliding windows from sequences, calculate MFE (in parallel), and index by central nucleotide position.

    Parameters:
        fasta_file (str): Path to FASTA file
        window_sizes (list): Window sizes to extract
        step (int): Step size for sliding windows
        n_jobs (int): Number of parallel jobs (-1 uses all available cores)

    Returns:
        dict: {seq_id: DataFrame} where rows = window sizes, columns = central positions,
              and values = MFE.
    """
    mfe_tables = {}

    # Count total sequences first for progress bar
    total_seqs = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

    for record in tqdm(SeqIO.parse(fasta_file, "fasta"), total=total_seqs, desc="Computing MFE"):
        seq_id = record.id
        seq = str(record.seq).upper()
        seq_len = len(seq)

        # Initialize a DataFrame with window sizes as rows
        table = pd.DataFrame(index=window_sizes)

        for w in tqdm(window_sizes, desc=f"Windows for {seq_id}", leave=False):
            # Extract all windows upfront
            windows = []
            center_positions = []

            for start in range(0, seq_len - w + 1, step):
                end = start + w
                window_seq = seq[start:end]
                windows.append(window_seq)

                # Central position
                center = start + w // 2
                center_positions.append(center)

            # Compute MFE for all windows in parallel
            mfe_values = joblib.Parallel(n_jobs=n_jobs)(
                joblib.delayed(_compute_mfe)(ws) for ws in windows
            )

            # Add MFE values for this window size as a row, indexed by central position
            table_row = pd.Series(mfe_values, index=center_positions)
            table = table.join(table_row.to_frame(name=w), how='outer')

        # Transpose so rows = window size, columns = center positions
        table = table.T
        mfe_tables[seq_id] = table

    return mfe_tables
