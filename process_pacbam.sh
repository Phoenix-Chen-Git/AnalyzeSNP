#!/bin/bash

# Combined script: process_new.sh (auto-detect ref FASTA and analyze mutation frequencies)

ulimit -v unlimited  # unlock mem

# === SETTINGS ===
DATA_DIR="../data"
OUT_DIR="../output"
SAMPLE_LIST="./sample_list.txt"
REF_DIR="../ref"
BT2_DIR="$REF_DIR/bt2_index"
SAM_DIR="$REF_DIR/sam_index"
TMP_DIR="../temp"
SCRIPT_DIR="."  # Directory where calc_mut_freq.py is located

# === Detect reference FASTA ===
fa_files=($(find "$REF_DIR" -maxdepth 1 -type f -name "*.fa"))
if [ ${#fa_files[@]} -eq 0 ]; then
  echo "No FASTA file found in $REF_DIR"
  exit 1
elif [ ${#fa_files[@]} -gt 1 ]; then
  echo "Too many FASTA files found in $REF_DIR. Please keep only one."
  printf "Found files:\n%s\n" "${fa_files[@]}"
  exit 1
fi

ref_file="${fa_files[0]}"
ref_base=$(basename "$ref_file" .fa)
REF_NAME="$ref_base"

# === Step 1: Discover unique sample prefixes ===
samples=""
for f in "$DATA_DIR"/*.fq.gz; do
  [ -e "$f" ] || continue
  fn=${f##*/}
  base=${fn%.fq.gz}
  samp=${base%_[12]}
  case " $samples " in
    *" $samp "*) ;;
    *) samples="$samples $samp" ;;
  esac
done

samples=${samples# }
echo $samples
printf "%s\n" $samples > "$SAMPLE_LIST"
echo "Found $(wc -l < "$SAMPLE_LIST") samples. Written to $SAMPLE_LIST."

# Allocate the computation resources
TOTAL_CPUS=$(nproc)
THREADS=$((TOTAL_CPUS - 2))
[ "$THREADS" -lt 1 ] && THREADS=1

echo "Using $THREADS threads (out of $TOTAL_CPUS total, sparing 2)"

# === Step 2: Run pipeline per sample ===
while read sample; do
  fq1="$DATA_DIR/${sample}_1.fq.gz"
  fq2="$DATA_DIR/${sample}_2.fq.gz"
  out_dir="$OUT_DIR/$sample"
  tmp_dir="$TMP_DIR/${sample}_sort_tmp"

  echo "Checking sample: $sample"

  if [[ -f "$fq1" && -f "$fq2" ]]; then
    echo "  Found: $fq1"
    echo "  Found: $fq2"
    echo "  Creating output dir: $out_dir"
    mkdir -p "$out_dir" "$tmp_dir"
    chmod 777 "$tmp_dir"

    trim_galore \
      -q 35 \
      --clip_R1 15 \
      --three_prime_clip_R1 10 \
      --clip_R2 15 \
      --three_prime_clip_R2 10 \
      --fastqc \
      --paired \
      --cores $THREADS \
      -o "$out_dir" \
      "$fq1" "$fq2"

    fq1_trimmed="$out_dir/${sample}_1_val_1.fq.gz"
    fq2_trimmed="$out_dir/${sample}_2_val_2.fq.gz"
    sam="$out_dir/${sample}.sam"
    bam="$out_dir/${sample}.bam"
    sorted_bam="$out_dir/${sample}.sorted.bam"

    echo 'align using bowtie2'
    bowtie2 \
      -p $THREADS \
      -t \
      -x "$BT2_DIR/$REF_NAME" \
      -1 "$fq1_trimmed" \
      -2 "$fq2_trimmed" \
      -S "$sam"

    echo 'sam2bam'
    samtools view --threads $THREADS -S "$sam" -b > "$bam"
    echo 'sort bam'
    samtools sort --threads $THREADS "$bam" -o "$sorted_bam"
    echo 'index bam'
    samtools index -@ $THREADS "$sorted_bam"

    bed_file="$out_dir/${sample}.regions.bed"
    vcf_file="$out_dir/${sample}.dummy.vcf"
    pacbam_dir="$out_dir/${sample}.pacbam.tsv"
    freq_out="$out_dir/${sample}.mut_freq.tsv"

    echo 'generate bed'
    awk '
      /^>/ {
        if (chr != "") {
          print chr, 0, seq_len
        }
        chr = substr($1, 2)
        seq_len = 0
        next
      }
      {
        seq_len += length($0)
      }
      END {
        if (chr != "") {
          print chr, 0, seq_len
        }
      }
    ' OFS="\t" "$SAM_DIR/${REF_NAME}.fa" > "$bed_file"

    echo 'false vcf'
    echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$vcf_file"
    awk -v OFS="\t" '{ print $1, $2 + 1, ".", "A", "T", ".", ".", "." }' "$bed_file" >> "$vcf_file"

    echo 'pacbam'
    pacbam bam="$sorted_bam" \
        bed="$bed_file" \
        vcf="$vcf_file" \
        fasta="$SAM_DIR/${REF_NAME}.fa" \
        mode=1 \
        out="$pacbam_dir" \
        threads="$THREADS"

    echo 'python'
    pileup_file=$(find "$pacbam_dir" -type f -name "*.pileup" | head -n 1)
    # Clean non-ASCII characters in pileup file before running Python
    clean_pileup="${pileup_file}.cleaned"
    tr -cd '\11\12\15\40-\176' < "$pileup_file" > "$clean_pileup"

    # Confirm the cleanup worked
    if [ ! -s "$clean_pileup" ]; then
      echo "  Error: Cleaned pileup file is empty or failed."
      continue
    fi

    # Run Python script on cleaned pileup
    python "$SCRIPT_DIR/readcounts.py" \
        -i "$clean_pileup" \
        -o "$freq_out"

  else
    echo "  Missing files for sample: $sample"
  fi

  echo "-------------------------"
done < "$SAMPLE_LIST"

echo "所有样本处理完成！"
