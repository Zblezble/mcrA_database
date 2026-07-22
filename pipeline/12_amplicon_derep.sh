#!/bin/bash
#PBS -N amplicon_derep
#PBS -l select=1:ncpus=2:mem=16gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 12_amplicon_derep.sh
# Trim the cleaned final alignment to short- and long-read amplicon regions, strip gaps, and dereplicate each region separately with RESCRIPt.
#
# Inputs:
#   alignment_final_cleaned.fasta
#   taxonomy_final_cleaned.tsv
#
# Outputs (unaligned dereplicated amplicon sequences):
#   short: seqs_short_derep.fasta + taxonomy_short_derep.tsv
#   long:  seqs_long_derep.fasta  + taxonomy_long_derep.tsv
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

ALN_IN="${STORAGE}/11_dereplicated_final/alignment_final_derep_guides.fasta"
TAX_IN="${STORAGE}/11_dereplicated_final/taxonomy_derep_final.tsv"

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${ALN_IN}" alignment_final_cleaned.fasta
cp "${TAX_IN}" taxonomy_final_cleaned.tsv

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimming to short and long amplicon regions ..."
python3 - <<'PY'
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

align = AlignIO.read('alignment_final_cleaned.fasta', 'fasta')


def compute_span(prefix):
    """Return (start, end) covering all non-gap positions for IDs with prefix."""
    positions = []
    for rec in align:
        if rec.id.startswith(prefix):
            seq = str(rec.seq)
            for col, ch in enumerate(seq):
                if ch not in '-.':
                    positions.append(col)
    if not positions:
        raise RuntimeError(f'No sequences found with prefix {prefix!r}')
    return min(positions), max(positions)


def trim_and_clean(start, end, out_fasta, out_tax):
    """Trim alignment to [start, end], drop guides/primers and all-gap seqs,
    and strip gap characters so the output is valid QIIME2 DNA FASTA."""
    keep_recs = []
    keep_ids = set()
    for rec in align:
        # Skip guide and primer sequences themselves
        if rec.id.startswith(('mcrA_', 'mlas_', 'guide_')):
            continue
        trimmed_seq = ''.join(rec.seq[i] for i in range(start, end + 1))
        # Skip sequences with no data in the trimmed region
        if all(c in '-.' for c in trimmed_seq):
            continue
        # Strip gaps for QIIME2 import
        ungapped_seq = trimmed_seq.upper().replace('-', '').replace('.', '')
        if not ungapped_seq:
            continue
        new_rec = rec.__class__(rec.seq.__class__(ungapped_seq), id=rec.id, description='')
        keep_recs.append(new_rec)
        keep_ids.add(rec.id)

    SeqIO.write(keep_recs, out_fasta, 'fasta')

    # Filter taxonomy to the kept IDs
    with open('taxonomy_final_cleaned.tsv') as fin, open(out_tax, 'w') as fout:
        header = fin.readline()
        fout.write(header)
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if parts[0] in keep_ids:
                fout.write(line)

    return len(keep_recs)


short_start, short_end = compute_span('guide_short_')
long_start, long_end = compute_span('guide_long_')

n_short = trim_and_clean(short_start, short_end,
                         'seqs_short.fasta', 'taxonomy_short.tsv')
n_long = trim_and_clean(long_start, long_end,
                        'seqs_long.fasta', 'taxonomy_long.tsv')

print(f'Short region: columns {short_start}-{short_end} '
      f'({short_end - short_start + 1} bp), {n_short} sequences')
print(f'Long region:  columns {long_start}-{long_end} '
      f'({long_end - long_start + 1} bp), {n_long} sequences')
PY

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dereplicating short and long amplicon datasets ..."
for region in short long; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] --- ${region} amplicon ---"

    # QIIME2 expects the taxonomy column to be named 'Taxon'
    sed '1s/new_tax/Taxon/' "taxonomy_${region}.tsv" > "taxonomy_${region}_qiime.tsv"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importing ${region} sequences ..."
    qiime tools import \
      --input-path "seqs_${region}.fasta" \
      --output-path "seqs_${region}.qza" \
      --type 'FeatureData[Sequence]'

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importing ${region} taxonomy ..."
    qiime tools import \
      --input-path "taxonomy_${region}_qiime.tsv" \
      --output-path "taxonomy_${region}.qza" \
      --type 'FeatureData[Taxonomy]'

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dereplicating ${region} ..."
    qiime rescript dereplicate --verbose \
      --i-sequences "seqs_${region}.qza" \
      --i-taxa "taxonomy_${region}.qza" \
      --p-mode uniq \
      --p-threads 1 \
      --p-derep-prefix \
      --o-dereplicated-sequences "seqs_${region}_derep.qza" \
      --o-dereplicated-taxa "taxonomy_${region}_derep.qza"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting ${region} dereplicated outputs ..."
    qiime tools export \
      --input-path "seqs_${region}_derep.qza" \
      --output-path "seqs_${region}_derep_exported"
    cp "seqs_${region}_derep_exported/dna-sequences.fasta" "seqs_${region}_derep.fasta"

    qiime tools export \
      --input-path "taxonomy_${region}_derep.qza" \
      --output-path "taxonomy_${region}_derep_exported"
    cp "taxonomy_${region}_derep_exported/taxonomy.tsv" "taxonomy_${region}_derep.tsv"
done

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/12_amplicon_derep"
mkdir -p "${OUT_STORAGE}"

for region in short long; do
    mkdir -p "${OUT_STORAGE}/${region}"
    cp "seqs_${region}.fasta"              "${OUT_STORAGE}/${region}/"
    cp "taxonomy_${region}.tsv"            "${OUT_STORAGE}/${region}/"
    cp "seqs_${region}.qza"                "${OUT_STORAGE}/${region}/"
    cp "taxonomy_${region}.qza"            "${OUT_STORAGE}/${region}/"
    cp "seqs_${region}_derep.qza"          "${OUT_STORAGE}/${region}/"
    cp "taxonomy_${region}_derep.qza"      "${OUT_STORAGE}/${region}/"
    cp "seqs_${region}_derep.fasta"        "${OUT_STORAGE}/${region}/"
    cp "taxonomy_${region}_derep.tsv"      "${OUT_STORAGE}/${region}/"
    cp -r "seqs_${region}_derep_exported"     "${OUT_STORAGE}/${region}/"
    cp -r "taxonomy_${region}_derep_exported" "${OUT_STORAGE}/${region}/"
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 12 complete. Results in ${OUT_STORAGE}"
echo "  short/seqs_short_derep.fasta + short/taxonomy_short_derep.tsv"
echo "  long/seqs_long_derep.fasta   + long/taxonomy_long_derep.tsv"
