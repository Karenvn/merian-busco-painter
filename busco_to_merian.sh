#!/bin/bash

# Master script to generate Merian plots from BUSCO results.
# Run from a project directory containing `tolids` plus an accession table.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUSCO_DIR="${BUSCO_DIR:-/Users/kh18/server_data/busco}"
OUTPUT_DIR="${OUTPUT_DIR:-/Users/kh18/server_data/merians}"
MERIAN_REF="${MERIAN_REF:-${SCRIPT_DIR}/Merian_elements_full_table.tsv}"
TOLID_FILE="${TOLID_FILE:-tolids}"

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

if [[ ! -f "$MERIAN_REF" ]]; then
  echo "ERROR: Merian reference table not found: $MERIAN_REF"
  exit 1
fi

if [[ ! -f "$TOLID_FILE" ]]; then
  echo "ERROR: $TOLID_FILE not found in current directory"
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

process_tolid() {
  local tolid="$1"
  local accession busco_input tolid_output location_file lengths_file

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
  if ! python3 "${SCRIPT_DIR}/buscopainter.py" \
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

  echo "Plotting Merian elements..."
  if ! python3 "${SCRIPT_DIR}/plot_buscopainter.py" \
    --file "$location_file" \
    --lengths "$lengths_file" \
    --prefix "${tolid_output}/${tolid}" \
    --minimum 1 \
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

while IFS= read -r tolid || [[ -n "$tolid" ]]; do
  [[ -z "$tolid" ]] && continue
  [[ "$tolid" =~ ^# ]] && continue

  total=$((total + 1))
  if process_tolid "$tolid"; then
    success=$((success + 1))
    continue
  fi

  status=$?
  if [[ $status -eq 2 ]]; then
    skipped=$((skipped + 1))
  else
    failed=$((failed + 1))
  fi
done < "$TOLID_FILE"

echo "================================"
echo "Batch complete"
echo "   Total:   $total"
echo "   Success: $success"
echo "   Skipped: $skipped"
echo "   Failed:  $failed"
echo "   Output:  $OUTPUT_DIR"
echo "================================"
