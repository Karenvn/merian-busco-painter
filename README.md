# merian-plotting

Utilities for plotting Merian element assignments from BUSCO `full_table.tsv`
files, adapted for genome notes workflows.

This code is derived from
[`charlottewright/lep_busco_painter`](https://github.com/charlottewright/lep_busco_painter)
and uses the same Merian reference table. The main changes in this version are:

- a Python plotting script (`plot_buscopainter.py`) instead of the upstream R plotter
- inclusion of duplicated BUSCO hits in the location table and plot
- chromosome length lookup from the NCBI Datasets API `sequence_reports` endpoint
- a batch wrapper for running multiple ToLIDs from genome notes working directories

## Files

- `buscopainter.py`: map BUSCO hits to Merian elements and optionally fetch chromosome lengths
- `plot_buscopainter.py`: generate PNG and SVG Merian plots
- `busco_to_merian.sh`: batch wrapper for ToLID lists
- `Merian_elements_full_table.tsv`: reference BUSCO-to-Merian table
- `LICENSE`: MIT license retained for the adapted codebase

## Requirements

- Python 3.10+
- `pandas`
- `matplotlib`
- `requests`
- `biopython`

Install with:

```bash
python3 -m pip install -r requirements.txt
```

## Single assembly workflow

Prepare a plot from a BUSCO result directory:

```bash
mkdir -p output/ilHelArmi9

python3 buscopainter.py \
  --reference_table Merian_elements_full_table.tsv \
  --query_table /Users/kh18/server_data/busco/ilHelArmi9/full_table.tsv \
  --prefix output/ilHelArmi9/ \
  --accession GCA_963930815.1

python3 plot_buscopainter.py \
  --file output/ilHelArmi9/all_location.tsv \
  --lengths output/ilHelArmi9/chrom_lengths.tsv \
  --prefix output/ilHelArmi9/ilHelArmi9 \
  --minimum 1 \
  --palette merianbow4 \
  --label-threshold 5
```

Outputs:

- `all_location.tsv`
- `chrom_lengths.tsv`
- `*.png`
- `*.svg`

## Batch workflow

The wrapper expects to be run from a directory containing:

- `tolids`: one ToLID per line
- `tolids_accessions.tsv` or `tolid_accessions.tsv`: two-column tab-separated file with `ToLID<TAB>assembly_accession`

Run:

```bash
bash busco_to_merian.sh
```

By default the wrapper reads BUSCO results from `/Users/kh18/server_data/busco`
and writes plots to `/Users/kh18/server_data/merians`. You can override those
paths with environment variables:

```bash
BUSCO_DIR=/path/to/busco \
OUTPUT_DIR=/path/to/output \
bash busco_to_merian.sh
```

## Notes

- `buscopainter.py` keeps both `Complete` and `Duplicated` BUSCO records.
- Chromosome lengths come from NCBI Datasets `sequence_reports`, using the main
  assembled-molecule accession and summing any unlocalized scaffolds assigned to
  that chromosome.
- The plotting script labels each scaffold/chromosome with Merian elements that
  meet the `--label-threshold`.
- If Open Sans is available locally it will be used automatically; otherwise
  matplotlib's default sans-serif font is used.

## Publishing

A minimal GitHub publication workflow from this directory is:

```bash
git init
git add README.md requirements.txt .gitignore LICENSE *.py *.sh Merian_elements_full_table.tsv
git commit -m "Initial import of Merian plotting scripts"
```

Then create a GitHub repository and add a remote:

```bash
git remote add origin git@github.com:YOUR-ORG/merian-plotting.git
git branch -M main
git push -u origin main
```

## Attribution

Please retain attribution to the upstream repository in any public release:

- Charlotte Wright, `lep_busco_painter`
- <https://github.com/charlottewright/lep_busco_painter>
