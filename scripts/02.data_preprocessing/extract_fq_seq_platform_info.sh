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

echo -e "sample\tr1\tr2\tr1_size\tr2_size\tflowcell\tlane" > "$OUT"

while read -r file; do

  fullpath="${FASTQ_DIR}/${file}"

#  if [[ "$file" =~ _1\.f(ast)?q(\.gz)?$ ]]; then
#      r1="$fullpath"
#      r2="${FASTQ_DIR}/$(echo "$file" | sed -E 's/_1(\.f(ast)?q(\.gz)?$)/_2\1/')"
#  elif [[ "$file" =~ _2\.f(ast)?q(\.gz)?$ ]]; then
#      continue   # skip R2 lines, we handle pairs from R1
#  else
#      continue
#  fi

  if [[ "$file" =~ _R1\.f(ast)?q(\.gz)?$ ]]; then
      r1="$fullpath"
      r2="${FASTQ_DIR}/$(echo "$file" | sed -E 's/_R1(\.f(ast)?q(\.gz)?$)/_R2\1/')"

  elif [[ "$file" =~ _1\.f(ast)?q(\.gz)?$ ]]; then
      r1="$fullpath"
      r2="${FASTQ_DIR}/$(echo "$file" | sed -E 's/_1(\.f(ast)?q(\.gz)?$)/_2\1/')"

  elif [[ "$file" =~ _R2\.f(ast)?q(\.gz)?$ || "$file" =~ _2\.f(ast)?q(\.gz)?$ ]]; then
      continue
  else
      continue
  fi

  if [[ ! -f "$r2" ]]; then
      echo "WARNING: missing mate for $r1" >&2
      continue
  fi

  r1_size=$(stat -c %s "$r1")
  r2_size=$(stat -c %s "$r2")

  header=$(zcat "$r1" | head -n 1)
  flowcell=$(echo "$header" | cut -d: -f3)
  lane=$(echo "$header" | cut -d: -f4)

  echo -e "${file}\t${r1}\t${r2}\t${r1_size}\t${r2_size}\t${flowcell}\t${lane}" >> "$OUT"

done < "$SAMPLE_LIST"
