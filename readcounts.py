import pandas as pd
import argparse

def calculate_mutation_frequencies(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")

    bases = "ACGT"
    sub_cols = [f"{b1}>{b2}" for b1 in bases for b2 in bases if b1 != b2]
    ref_to_cols = [f"ref>{b}" for b in bases]
    all_new_cols = sub_cols + ref_to_cols + ["ref>N"]

    for col in all_new_cols:
        df[col] = 0.0

    for i, row in df.iterrows():
        ref = row["ref"]
        cov = row["cov"]
        if cov == 0 or ref not in bases:
            continue

        total_mut_freq = 0.0
        refN = 0.0

        for b in bases:
            freq = row[b] / cov
            if b != ref:
                total_mut_freq += freq
                refN += freq
                df.at[i, f"{ref}>{b}"] = round(freq * 100, 6)
            df.at[i, f"ref>{b}"] = round(freq * 100, 6)

        df.at[i, "ref>N"] = round(refN * 100, 6)

    df.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Mutation frequencies calculated. Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mutation frequencies from pacbam output")
    parser.add_argument("-i", "--input", required=True, help="Input pacbam .pileup file")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")

    args = parser.parse_args()
    calculate_mutation_frequencies(args.input, args.output)
