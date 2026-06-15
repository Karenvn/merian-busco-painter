"""Map BUSCO full_table rows to Merian elements and optional chromosome lengths."""

from __future__ import annotations

import csv
import os
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import requests

API_KEY = os.getenv("NCBI_API_KEY")


@dataclass(frozen=True)
class PaintOutputs:
    all_locations: Path
    chrom_lengths: Path
    summary: Path
    wrote_lengths: bool
    wrote_summary: bool


def parse_busco_table(path: Path) -> tuple[list[tuple[str, str, int, int]], list[str]]:
    """Return BUSCO ID, chromosome, start and stop for Complete/Duplicated rows."""
    table: list[tuple[str, str, int, int]] = []
    chromosomes: set[str] = set()
    keep_status = {"Complete", "Duplicated"}

    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 5:
                continue
            busco_id, status, chrom, start, stop = row[:5]
            if status not in keep_status:
                continue
            try:
                start_coord, end_coord = int(start), int(stop)
            except ValueError:
                continue
            table.append((busco_id, chrom, start_coord, end_coord))
            chromosomes.add(chrom)
    return table, sorted(chromosomes)


def build_ref_map(ref_path: Path) -> dict[str, str]:
    """Return BUSCO-ID to Merian element mapping from the reference table."""
    merian_set = {"MZ"} | {f"M{i}" for i in range(1, 32)}
    ref_map: dict[str, str] = {}
    with ref_path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if (
                not row
                or len(row) < 3
                or row[0].startswith("#")
                or row[0].lower().startswith("busco")
            ):
                continue
            merian = row[2].upper().strip()
            if merian in merian_set:
                ref_map[row[0]] = merian
    return ref_map


def build_location_rows(
    ref_map: dict[str, str], query_table: list[tuple[str, str, int, int]]
) -> list[str]:
    rows = ["buscoID\tquery_chr\tposition\tassigned_chr\tstatus"]
    for busco_id, query_chr, start, end in query_table:
        position = (start + end) / 2
        assigned = ref_map.get(busco_id, "NA")
        rows.append(f"{busco_id}\t{query_chr}\t{position}\t{assigned}\t{assigned}")
    return rows


def fetch_sequence_report(accession: str) -> list[dict]:
    url = (
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
        f"{accession}/sequence_reports"
    )
    headers = {"accept": "application/json", "User-Agent": "buscopainter"}
    params = {}
    if API_KEY:
        params["api_key"] = API_KEY
    print(f"[INFO] Fetching chromosome info from NCBI for {accession}...")
    response = requests.get(url, headers=headers, params=params, timeout=60)
    response.raise_for_status()
    payload = response.json()
    return payload.get("sequence_report", {}).get("records") or payload.get(
        "reports", []
    )


def chrom_lengths_with_unloc(records: list[dict]) -> list[tuple[str, int]]:
    """Return main chromosome accession and bp length including unlocalized scaffolds."""
    main_acc: dict[str, str] = {}
    for rec in records:
        if (
            rec.get("role") == "assembled-molecule"
            and rec.get("assigned_molecule_location_type") == "Chromosome"
        ):
            main_acc[rec["chr_name"]] = rec["genbank_accession"]

    bp_tot: dict[str, int] = {acc: 0 for acc in main_acc.values()}
    for rec in records:
        role = rec.get("role")
        loc = rec.get("assigned_molecule_location_type", "")
        if role == "assembled-molecule" and loc == "Chromosome":
            acc = rec["genbank_accession"]
            bp_tot[acc] += int(rec.get("length", 0))
        elif role == "unlocalized-scaffold":
            parent = rec.get("chr_name")
            acc = main_acc.get(parent)
            if acc:
                bp_tot[acc] += int(rec.get("length", 0))

    return sorted(bp_tot.items(), key=lambda item: -item[1])


def write_tsv(lines: list[str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def resolve_output_paths(prefix: str | Path) -> tuple[Path, Path, Path]:
    """Resolve output paths from either a directory-like prefix or file stem."""
    prefix_text = str(prefix)
    prefix_path = Path(prefix)
    if prefix_text.endswith(("/", "\\")) or prefix_path.is_dir():
        out_dir = prefix_path
        return (
            out_dir / "all_location.tsv",
            out_dir / "chrom_lengths.tsv",
            out_dir / "summary.tsv",
        )

    out_dir = prefix_path.parent
    stem = prefix_path.name
    return (
        out_dir / f"{stem}_all_location.tsv",
        out_dir / f"{stem}_chrom_lengths.tsv",
        out_dir / f"{stem}_summary.tsv",
    )


def paint_buscos(
    reference_table: Path,
    query_table: Path,
    prefix: str | Path,
    accession: str | None = None,
    write_summary: bool = False,
) -> PaintOutputs:
    """Run the painter workflow and return generated output paths."""
    out_all, out_len, out_sum = resolve_output_paths(prefix)
    out_all.parent.mkdir(parents=True, exist_ok=True)

    ref_map = build_ref_map(reference_table)
    query_rows, query_chromosomes = parse_busco_table(query_table)
    all_rows = build_location_rows(ref_map, query_rows)
    chrom_order = query_chromosomes.copy()

    wrote_len = False
    if accession:
        pairs = chrom_lengths_with_unloc(fetch_sequence_report(accession))
        length_lines = ["Chrom\tLength_Mb"] + [
            f"{chrom}\t{basepairs / 1e6:.3f}" for chrom, basepairs in pairs
        ]
        write_tsv(length_lines, out_len)
        wrote_len = True
        chrom_order = [chrom for chrom, _ in pairs]

    query_chroms = {chrom for _, chrom, _, _ in query_rows}
    missing = [chrom for chrom in chrom_order if chrom not in query_chroms]
    for chrom in missing:
        all_rows.append(f"NA\t{chrom}\tNA\tNA\tNA")

    write_tsv(all_rows, out_all)

    wrote_sum = False
    if write_summary:
        counts = Counter(chrom for _, chrom, _, _ in query_rows)
        counts.update({chrom: 0 for chrom in missing})
        summary_lines = ["query_chr\tbusco_hits"] + [
            f"{chrom}\t{counts[chrom]}" for chrom in chrom_order
        ]
        write_tsv(summary_lines, out_sum)
        wrote_sum = True

    return PaintOutputs(
        all_locations=out_all,
        chrom_lengths=out_len,
        summary=out_sum,
        wrote_lengths=wrote_len,
        wrote_summary=wrote_sum,
    )
