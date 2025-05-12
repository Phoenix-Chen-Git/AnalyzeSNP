#!/bin/bash

# Combined script: process_new.sh (auto-detect ref FASTA)

ulimit -v unlimited  # unlock mem

# === SETTINGS ===
DATA_DIR="../data"
OUT_DIR="../output"
SAMPLE_LIST="./sample_list.txt"
REF_DIR="../ref"
BT2_DIR="$REF_DIR/bt2_index"
SAM_DIR="$REF_DIR/sam_index"
TMP_DIR="../temp"



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

# Alocare the computation resourses
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
      -q 30 \
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
    mpileup="$out_dir/${sample}.mpileup"
    tmp_dir="$TMP_DIR/${sample}_sort_tmp"

    bowtie2 \
      -p $THREADS \
      -t \
      -x "$BT2_DIR/$REF_NAME" \
      -1 "$fq1_trimmed" \
      -2 "$fq2_trimmed" \
      -S "$sam"


     # 处理流程
    samtools view --threads $THREADS -S "$sam" -b > "$bam" 
    samtools sort \
      --threads 12 \
      "$bam" -o "$sorted_bam"
    samtools index -@ $THREADS "$sorted_bam"

    # cleanup temp
    echo "清理临时目录: $tmp_dir"
    rm -rf "$tmp_dir"

    samtools mpileup \
      -d 5000000 \
      -f "$SAM_DIR/${REF_NAME}.fa" \
      "$sorted_bam" > "$mpileup"

    varscan readcounts "$mpileup" \
    --min-coverage 1000 \
    --output-file "$out_dir/${sample}.varScan.readcounts"

    #  parse_readcounts
    python parse_readcounts.py \
        "$out_dir/${sample}.varScan.readcounts" \
        "$sample" \
        "$out_dir"

    varscan mpileup2snp "$mpileup" \
    --min-var-freq 0.000 \
    --min-reads2 4 > "$out_dir/${sample}.varScan.snp"
    #varscan mpileup2snp "$mpileup" --min-var-freq 0.000 --min-reads2 4 --output-vcf 1 > "$out_dir/${sample}.vcf"

    # cleanup temp
    rm -rf "$tmp_dir"
    echo "Success: $sample"

  else
    echo "  Missing files for sample: $sample"
  fi

  echo "-------------------------"
done < "$SAMPLE_LIST"

echo "所有样本处理完成！"
