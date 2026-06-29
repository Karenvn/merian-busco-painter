#!/bin/bash

# Master script to generate Merian plots from BUSCO results.
# Run from a project directory containing `tolids` plus an accession table.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${DATA_ROOT:-.}"
BUSCO_DIR="${BUSCO_DIR:-${DATA_ROOT}/busco}"
OUTPUT_DIR="${OUTPUT_DIR:-${DATA_ROOT}/merians}"
MERIAN_REF="${MERIAN_REF:-${SCRIPT_DIR}/Merian_elements_full_table.tsv}"
TOLID_FILE="${TOLID_FILE:-}"

if command -v merian-busco-painter >/dev/null 2>&1; then
  PAINTER_CMD=(merian-busco-painter paint)
  PLOTTER_CMD=(merian-busco-painter plot)
else
  PAINTER_CMD=(python3 "${SCRIPT_DIR}/buscopainter.py")
  PLOTTER_CMD=(python3 "${SCRIPT_DIR}/plot_buscopainter.py")
fi

resolve_accession_file() {
  if [[ -n "${ACCESSION_FILE:-}" ]]; then
    printf '%s\n' "$ACCESSION_FILE"
    return 0
  fi
  if [[ -f "tolids_accessions.tsv" ]]; then
    printf '%s\n' "tolids_accessions.tsv"
    return 0
  fi
  if [[ -f "tolid_accessions.tsv" ]]; then
    printf '%s\n' "tolid_accessions.tsv"
    return 0
  fi
  return 1
}

ACCESSION_FILE="$(resolve_accession_file)" || {
  echo "ERROR: no accession table found. Expected ACCESSION_FILE, tolids_accessions.tsv, or tolid_accessions.tsv in $(pwd)"
  exit 1
}

TOLID_SOURCE=""
TOLID_SOURCE_MODE="list"
if [[ -n "$TOLID_FILE" ]]; then
  if [[ ! -f "$TOLID_FILE" ]]; then
    echo "ERROR: TOLID_FILE not found: $TOLID_FILE"
    exit 1
  fi
  TOLID_SOURCE="$TOLID_FILE"
else
  TOLID_SOURCE="$ACCESSION_FILE"
  TOLID_SOURCE_MODE="accession"
fi

if [[ ! -f "$MERIAN_REF" ]]; then
  echo "ERROR: Merian reference table not found: $MERIAN_REF"
  exit 1
fi

mkdir -p "$OUTPUT_DIR" || {
  echo "ERROR: could not create output directory: $OUTPUT_DIR"
  exit 1
}

get_accession() {
  local tolid="$1"
  awk -F '\t' -v tolid="$tolid" '$1 == tolid {print $2; exit}' "$ACCESSION_FILE"
}

length_table_rows() {
  local lengths_file="$1"
  awk -F '\t' 'NR > 1 && NF >= 2 && $1 != "" && $2 != "" {count++} END {print count + 0}' "$lengths_file"
}

find_local_fai() {
  local tolid="$1"
  find "${BUSCO_DIR}/${tolid}" -maxdepth 1 -type f \
    \( -name "*.fai" -o -name "*.fa.fai" -o -name "*.fasta.fai" -o -name "*.fa.gz.fai" -o -name "*.fasta.gz.fai" \) \
    2>/dev/null | sort | head -n 1
}

process_tolid() {
  local tolid="$1"
  local accession busco_input tolid_output location_file lengths_file local_fai n_lengths

  echo "================================"
  echo "Processing: $tolid"
  echo "================================"

  accession="$(get_accession "$tolid")"
  if [[ -z "$accession" ]]; then
    echo "WARN: No accession found for $tolid in $ACCESSION_FILE. Skipping."
    return 2
  fi

  echo "Accession: $accession"

  busco_input="${BUSCO_DIR}/${tolid}/full_table.tsv"
  if [[ ! -f "$busco_input" ]]; then
    echo "WARN: Missing BUSCO file: $busco_input. Skipping."
    return 2
  fi

  tolid_output="${OUTPUT_DIR}/${tolid}"
  mkdir -p "$tolid_output" || {
    echo "ERROR: Could not create output directory: $tolid_output"
    return 1
  }

  echo "Running BUSCO painter..."
  if ! "${PAINTER_CMD[@]}" \
    --reference_table "$MERIAN_REF" \
    --query_table "$busco_input" \
    --prefix "${tolid_output}/" \
    --accession "$accession"; then
    echo "ERROR: BUSCO painter failed for $tolid"
    return 1
  fi

  location_file="${tolid_output}/all_location.tsv"
  lengths_file="${tolid_output}/chrom_lengths.tsv"

  if [[ ! -f "$location_file" ]]; then
    echo "ERROR: BUSCO painter did not create $location_file"
    return 1
  fi

  if [[ ! -f "$lengths_file" ]]; then
    echo "ERROR: BUSCO painter did not create $lengths_file"
    return 1
  fi

  n_lengths="$(length_table_rows "$lengths_file")"
  if [[ "$n_lengths" -eq 0 ]]; then
    local_fai="$(find_local_fai "$tolid")"
    if [[ -z "$local_fai" ]]; then
      echo "ERROR: NCBI returned no chromosome lengths and no local .fai was found for $tolid in ${BUSCO_DIR}/${tolid}"
      return 1
    fi
    echo "WARN: NCBI returned no chromosome lengths for $accession; plotting with local .fai: $local_fai"
    lengths_file="$local_fai"
  fi

  echo "Plotting Merian elements..."
  rm -f "${tolid_output}/${tolid}.png" "${tolid_output}/${tolid}.svg"
  if ! "${PLOTTER_CMD[@]}" \
    --file "$location_file" \
    --lengths "$lengths_file" \
    --assembly-mode auto \
    --prefix "${tolid_output}/${tolid}" \
    --palette merianbow4 \
    --label-threshold 5; then
    echo "ERROR: Plotting failed for $tolid"
    return 1
  fi

  echo "Finished: $tolid"
  echo ""
  return 0
}

total=0
success=0
failed=0
skipped=0

if [[ "$TOLID_SOURCE_MODE" == "accession" ]]; then
  echo "ToLID source: first column of $TOLID_SOURCE"
else
  echo "ToLID source: $TOLID_SOURCE"
fi

while IFS= read -r line || [[ -n "$line" ]]; do
  if [[ "$TOLID_SOURCE_MODE" == "accession" ]]; then
    tolid="${line%%$'\t'*}"
  else
    tolid="$line"
  fi

  [[ -z "$tolid" ]] && continue
  [[ "$tolid" =~ ^# ]] && continue

  total=$((total + 1))
  process_tolid "$tolid"
  status=$?

  if [[ $status -eq 0 ]]; then
    success=$((success + 1))
  elif [[ $status -eq 2 ]]; then
    skipped=$((skipped + 1))
  else
    failed=$((failed + 1))
  fi
done < "$TOLID_SOURCE"

echo "================================"
echo "Batch complete"
echo "   Total:   $total"
echo "   Success: $success"
echo "   Skipped: $skipped"
echo "   Failed:  $failed"
echo "   Output:  $OUTPUT_DIR"
echo "================================"
