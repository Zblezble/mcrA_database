#!/bin/bash
#PBS -N dereplicate_final
#PBS -l select=1:ncpus=2:mem=16gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 11_dereplicate_final.sh
# Dereplicate the manually cleaned final database with RESCRIPt (mode=uniq).
#   1. Strip gap characters from the aligned cleaned FASTA so QIIME2 accepts it.
#   2. Import unaligned sequences + taxonomy to QIIME2.
#   3. Dereplicate with RESCRIPt to obtain representative IDs.
#   4. Use the dereplicated IDs to subset the original aligned FASTA and produce a dereplicated, aligned final database.
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

ALN_IN="${STORAGE}/10_manual_removed/alignment_final_cleaned.fasta"
TAXONOMY="${STORAGE}/10_manual_removed/taxonomy_final_cleaned.tsv"

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${ALN_IN}" alignment_final_cleaned.fasta
cp "${TAXONOMY}" taxonomy_final_cleaned.tsv

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Preparing unaligned FASTA and QIIME2 taxonomy ..."
python3 - <<'PY'
from Bio import SeqIO

# Read taxonomy IDs
with open('taxonomy_final_cleaned.tsv') as f:
    header = f.readline()
    tax_ids = {line.split('\t')[0].strip() for line in f}
print(f'  Taxonomy IDs: {len(tax_ids)}')

# Filter alignment and strip gaps for QIIME2 import.
# Also drop guide/primer sequences, which are not in the taxonomy.
ungapped_recs = []
n_filtered = 0
for rec in SeqIO.parse('alignment_final_cleaned.fasta', 'fasta'):
    if rec.id.startswith(('mcrA_', 'mlas_', 'guide_')):
        n_filtered += 1
        continue
    if rec.id not in tax_ids:
        n_filtered += 1
        continue
    seq = str(rec.seq).upper().replace('-', '').replace('.', '')
    if not seq:
        n_filtered += 1
        continue
    ungapped_recs.append(rec.__class__(rec.seq.__class__(seq), id=rec.id, description=''))

SeqIO.write(ungapped_recs, 'seqs_ungapped.fasta', 'fasta')
print(f'  Ungapped sequences for dereplication: {len(ungapped_recs)}')
print(f'  Filtered out (guides/primers/missing-tax/empty): {n_filtered}')

# Prepare taxonomy for QIIME2 (column must be named 'Taxon')
with open('taxonomy_final_cleaned.tsv') as f, open('taxonomy_qiime.tsv', 'w') as fout:
    h = f.readline()
    fout.write(h.replace('new_tax', 'Taxon'))
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if parts[0] in tax_ids:
            fout.write(line)
PY

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importing unaligned sequences to QIIME2 ..."
qiime tools import \
  --input-path seqs_ungapped.fasta \
  --output-path seqs_ungapped.qza \
  --type 'FeatureData[Sequence]'

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importing taxonomy to QIIME2 ..."
qiime tools import \
  --input-path taxonomy_qiime.tsv \
  --output-path taxonomy_qiime.qza \
  --type 'FeatureData[Taxonomy]'

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dereplicating (mode=uniq) ..."
qiime rescript dereplicate --verbose \
  --i-sequences seqs_ungapped.qza \
  --i-taxa taxonomy_qiime.qza \
  --p-mode uniq \
  --p-threads 1 \
  --p-derep-prefix \
  --o-dereplicated-sequences seqs_derep_final.qza \
  --o-dereplicated-taxa taxonomy_derep_final.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting dereplicated outputs ..."
qiime tools export \
  --input-path seqs_derep_final.qza \
  --output-path seqs_derep_final_exported

cp seqs_derep_final_exported/dna-sequences.fasta seqs_derep_final_ungapped.fasta

qiime tools export \
  --input-path taxonomy_derep_final.qza \
  --output-path taxonomy_derep_final_exported

cp taxonomy_derep_final_exported/taxonomy.tsv taxonomy_derep_final.tsv

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Subsetting aligned FASTA to dereplicated IDs ..."
python3 - <<'PY'
from Bio import SeqIO

# Get representative IDs retained by RESCRIPt
with open('taxonomy_derep_final.tsv') as f:
    next(f)
    kept_ids = {line.split('\t')[0].strip() for line in f}
print(f'  Dereplicated IDs: {len(kept_ids)}')

# Subset original aligned FASTA (drop guides/primers and duplicates)
aligned_kept = []
for rec in SeqIO.parse('alignment_final_cleaned.fasta', 'fasta'):
    if rec.id.startswith(('mcrA_', 'mlas_', 'guide_')):
        continue
    if rec.id in kept_ids:
        aligned_kept.append(rec)

SeqIO.write(aligned_kept, 'alignment_final_derep.fasta', 'fasta')
print(f'  Aligned sequences kept: {len(aligned_kept)}')
PY

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/11_dereplicated_final"
mkdir -p "${OUT_STORAGE}"

cp "${WORKDIR}/seqs_ungapped.fasta"              "${OUT_STORAGE}/"
cp "${WORKDIR}/seqs_ungapped.qza"                "${OUT_STORAGE}/"
cp "${WORKDIR}/taxonomy_qiime.tsv"               "${OUT_STORAGE}/"
cp "${WORKDIR}/taxonomy_qiime.qza"               "${OUT_STORAGE}/"
cp "${WORKDIR}/seqs_derep_final.qza"             "${OUT_STORAGE}/"
cp "${WORKDIR}/taxonomy_derep_final.qza"         "${OUT_STORAGE}/"
cp "${WORKDIR}/seqs_derep_final_ungapped.fasta"  "${OUT_STORAGE}/"
cp "${WORKDIR}/taxonomy_derep_final.tsv"         "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_final_derep.fasta"      "${OUT_STORAGE}/"
cp -r "${WORKDIR}/seqs_derep_final_exported"     "${OUT_STORAGE}/"
cp -r "${WORKDIR}/taxonomy_derep_final_exported" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 11 complete. Results in ${OUT_STORAGE}"
echo "  seqs_derep_final_ungapped.fasta  -- unaligned dereplicated sequences"
echo "  taxonomy_derep_final.tsv         -- dereplicated taxonomy"
echo "  alignment_final_derep.fasta      -- aligned dereplicated sequences"
