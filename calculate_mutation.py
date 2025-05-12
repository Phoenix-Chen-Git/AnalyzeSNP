import sys
import re

def calculate_mutation_rates(mpileup_file_path):
    """
    Reads an mpileup file and calculates the rate of each base substitution.

    Args:
        mpileup_file_path (str): Path to the input mpileup file.
    """
    # Initialize a dictionary to store counts of each substitution type
    # Keys are in the format 'Ref_Base->Obs_Base'
    mutation_counts = {
        'A->C': 0, 'A->G': 0, 'A->T': 0, 'A->N': 0,
        'C->A': 0, 'C->G': 0, 'C->T': 0, 'C->N': 0,
        'G->A': 0, 'G->C': 0, 'G->T': 0, 'G->N': 0,
        'T->A': 0, 'T->C': 0, 'T->G': 0, 'T->N': 0,
        'N->A': 0, 'N->C': 0, 'N->G': 0, 'N->T': 0, 'N->N': 0 # Include N as reference
    }

    # Counter for the total number of bases observed at reference positions
    # This excludes indels, skips, read starts/ends, etc.
    total_observed_bases = 0

    print(f"Processing mpileup file: {mpileup_file_path}")

    try:
        with open(mpileup_file_path, 'r') as f:
            for line in f:
                # Skip header lines if any (lines starting with @) or empty lines
                if line.startswith('@') or not line.strip():
                    continue

                fields = line.strip().split('\t')

                # Mpileup lines typically have at least 6 columns if quality scores are present
                # We need at least 4 (chrom, pos, ref, depth) and the pileup string (5th column)
                if len(fields) < 5:
                    # print(f"Skipping malformed line: {line.strip()}", file=sys.stderr)
                    continue

                # Extract reference base and pileup string
                # Reference base is the 3rd column (index 2)
                # Pileup string is the 5th column (index 4)
                ref_base = fields[2].upper() # Convert reference to uppercase for consistency
                pileup_string = fields[4]

                # Iterate through the pileup string to parse bases and indels
                i = 0
                while i < len(pileup_string):
                    char = pileup_string[i]

                    if char == '^':
                        # Skip read start marker and the following mapping quality character
                        i += 2
                        continue
                    elif char == '$':
                        # Skip read end marker
                        i += 1
                        continue
                    elif char == '+':
                        # Handle insertion: parse length and skip inserted sequence
                        i += 1 # Skip '+'
                        # Find the integer length
                        indel_len_str = ''
                        while i < len(pileup_string) and pileup_string[i].isdigit():
                            indel_len_str += pileup_string[i]
                            i += 1
                        try:
                            indel_len = int(indel_len_str)
                            # Skip the inserted sequence
                            i += indel_len
                        except ValueError:
                            # Malformed indel length, try to recover by skipping one char
                            # print(f"Warning: Malformed insertion length at position {fields[1]}. Skipping.", file=sys.stderr)
                            i += 1 # Skip the non-digit character
                        continue
                    elif char == '-':
                        # Handle deletion: parse length and skip deleted sequence
                        i += 1 # Skip '-'
                        # Find the integer length
                        indel_len_str = ''
                        while i < len(pileup_string) and pileup_string[i].isdigit():
                            indel_len_str += pileup_string[i]
                            i += 1
                        try:
                            indel_len = int(indel_len_str)
                            # Skip the deleted sequence
                            i += indel_len
                        except ValueError:
                            # Malformed indel length, try to recover by skipping one char
                            # print(f"Warning: Malformed deletion length at position {fields[1]}. Skipping.", file=sys.stderr)
                            i += 1 # Skip the non-digit character
                        continue
                    elif char == '*' or char == '<' or char == '>':
                         # Skip deletion placeholder, reference skip markers
                         i += 1
                         continue
                    elif char == '.' or char == ',':
                        # Match to reference base
                        # We still count this as an observed base, but not a mutation
                        total_observed_bases += 1
                        i += 1
                        continue
                    elif char.upper() in ['A', 'C', 'G', 'T', 'N']:
                        # This is an observed base (potential mismatch)
                        observed_base = char.upper() # Convert observed base to uppercase

                        # Only count as a mutation if it's different from the reference
                        if ref_base in mutation_counts and observed_base != ref_base:
                             mutation_key = f"{ref_base}->{observed_base}"
                             if mutation_key in mutation_counts:
                                 mutation_counts[mutation_key] += 1
                             # else:
                                 # print(f"Warning: Unexpected mutation key {mutation_key}. Skipping.", file=sys.stderr)

                        # Count this as an observed base regardless of match/mismatch
                        total_observed_bases += 1
                        i += 1
                        continue
                    else:
                        # Handle any other unexpected characters by skipping
                        # print(f"Warning: Skipping unexpected character '{char}' at position {fields[1]} in pileup string.", file=sys.stderr)
                        i += 1
                        continue

    except FileNotFoundError:
        print(f"Error: File not found at {mpileup_file_path}", file=sys.stderr)
        return
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        return

    # --- Calculate and Print Rates ---

    print("\n--- Mutation Rates ---")

    if total_observed_bases == 0:
        print("No bases observed with coverage to calculate mutation rates.")
        return

    # Sort mutation types alphabetically for consistent output
    sorted_mutation_types = sorted(mutation_counts.keys())

    for mutation_type in sorted_mutation_types:
        count = mutation_counts[mutation_type]
        # Calculate rate as a percentage
        rate = (count / total_observed_bases) * 100 if total_observed_bases > 0 else 0
        print(f"{mutation_type}: {rate:.4f}%") # Format to 4 decimal places

    print(f"\nTotal bases observed at reference positions: {total_observed_bases}")


# --- Script Entry Point ---
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python your_script_name.py <mpileup_file>")
        sys.exit(1)

    mpileup_file = sys.argv[1]
    calculate_mutation_rates(mpileup_file)
