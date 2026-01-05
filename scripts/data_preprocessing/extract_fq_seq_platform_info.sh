#!/usr/bin/env bash

usage() {
  echo "Usage: $0 -s sample_list.txt -f fastq_dir -o output.tsv"
  exit 1
}

# defaults (optional)
FASTQ_DIR="."

while getopts ":s:f:o:" opt; do
  case ${opt} in
    s ) SAMPLE_LIST=$OPTARG ;;
    f ) FASTQ_DIR=$OPTARG ;;
    o ) OUT=$OPTARG ;;
    * ) usage ;;
  esac
done

# required arguments
if [[ -z "${SAMPLE_LIST:-}" || -z "${OUT:-}" ]]; then
  usage
fi

echo -e "${sample}\t${r1}\t${r2}\t${r1_size}\t${r2_size}\t${flowcell}\t${lane}" > "$OUT"

cat $SAMPLE_LIST | while read sample; do 

  r1=$(ls "${FASTQ_DIR}"/"${sample}"* 2>/dev/null | \
       grep -E '(_R1_|_R1\.|_1\.fastq)' | head -n 1 || true)

  r2=$(ls "${FASTQ_DIR}"/"${sample}"* 2>/dev/null | \
       grep -E '(_R2_|_R2\.|_2\.fastq)' | head -n 1 || true)

  if [[ -z "$r1" || -z "$r2" ]]; then
    echo "WARNING: missing FASTQs for sample $sample" >&2
    continue
  fi

  r1_size=$(stat -c %s "$r1")
  r2_size=$(stat -c %s "$r2")

  # Read first FASTQ header
  header=$(zcat "$r1" | head -n 1)

  flowcell=$(echo $header | cut -d: -f3)
  lane=$(echo $header | cut -d: -f4)

  echo -e "${sample}\t${r1}\t${r2}\t${r1_size}\t${r2_size}\t${flowcell}\t${lane}" >> "$OUT" ;

  done
