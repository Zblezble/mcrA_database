#!/bin/bash
#PBS -N get_ncbi_mcrA
#PBS -l select=1:ncpus=8:mem=32gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail
# ============================================================================
# 01_get_ncbi_mcrA_mrtA.sh
# Fetch mcrA and mrtA nucleotide sequences from NCBI using RESCRIPt
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

declare -A Q
Q['mcrA']='txid2157[ORGN] AND (mcrA[Title] OR coenzyme-b sulfoethylthiotransferase subunit alpha[Title] OR methylcoenzyme M reductase subunit A[Title] OR methylcoenzyme M reductase subunit alpha[Title] OR methylcoenzyme M reductase alpha subunit[Title] OR methylcoenzyme M reductase alpha subunit-like[Title] OR methyl coenzyme M reductase subunit A[Title] OR methyl coenzyme M reductase subunit alpha[Title] OR methyl coenzyme M reductase alpha subunit[Title] OR methyl coenzyme M reductase alpha subunit-like[Title] OR methyl-coenzyme M reductase subunit A[Title] OR methyl-coenzyme M reductase alpha subunit[Title] OR methyl-coenzyme M reductase subunit alpha[Title] OR methyl-coenzyme M reductase alpha subunit-like[Title] OR methylcoenzyme M reductase, subunit A[Title] OR methylcoenzyme M reductase, subunit alpha[Title] OR methylcoenzyme M reductase, alpha subunit[Title] OR methylcoenzyme M reductase, alpha subunit-like[Title] OR methyl coenzyme M reductase, subunit A[Title] OR methyl coenzyme M reductase, subunit alpha[Title] OR methyl coenzyme M reductase, alpha subunit[Title] OR methyl coenzyme M reductase, alpha subunit-like[Title] OR methyl-coenzyme M reductase, subunit A[Title] OR methyl-coenzyme M reductase, alpha subunit[Title] OR methyl-coenzyme M reductase, subunit alpha[Title] OR methyl-coenzyme M reductase, alpha subunit-like[Title]) NOT (euryarchaeota archaeon[Title] OR methanogenic archaeon enrichment culture[Title] OR uncultured[Title])'

Q['mrtA']='txid2157[ORGN] AND (mrtA OR methylcoenzyme M reductase II subunit A[Title] OR methylcoenzyme M reductase II subunit alpha[Title] OR methylcoenzyme M reductase II alpha subunit[Title] OR methylcoenzyme M reductase II alpha subunit-like[Title] OR methyl coenzyme M reductase II subunit A[Title] OR methyl coenzyme M reductase II subunit alpha[Title] OR methyl coenzyme M reductase II alpha subunit[Title] OR methyl coenzyme M reductase II alpha subunit-like[Title] OR methyl-coenzyme M reductase II subunit A[Title] OR methyl-coenzyme M reductase II alpha subunit[Title] OR methyl-coenzyme M reductase II subunit alpha[Title] OR methyl-coenzyme M reductase II alpha subunit-like[Title] OR methylcoenzyme M reductase, II subunit A[Title] OR methylcoenzyme M reductase, II subunit alpha[Title] OR methylcoenzyme M reductase, II alpha subunit[Title] OR methylcoenzyme M reductase, II alpha subunit-like[Title] OR methyl coenzyme M reductase, II subunit A[Title] OR methyl coenzyme M reductase, II subunit alpha[Title] OR methyl coenzyme M reductase, II alpha subunit[Title] OR methyl coenzyme M reductase, II alpha subunit-like[Title] OR methyl-coenzyme M reductase, II subunit A[Title] OR methyl-coenzyme M reductase, II alpha subunit[Title] OR methyl-coenzyme M reductase, II subunit alpha[Title] OR methyl-coenzyme M reductase, II alpha subunit-like[Title]) NOT (euryarchaeota archaeon[Title] OR methanogenic archaeon enrichment culture[Title] OR uncultured[Title] OR genome[Title])'

# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

# Activate QIIME2
module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fetching mcrA sequences from NCBI ..."
qiime rescript get-ncbi-data \
  --p-query "${Q['mcrA']}" \
  --p-ranks domain kingdom phylum class order family genus species \
  --verbose \
  --p-n-jobs 8 \
  --o-sequences mcrA_ncbi_seqs.qza \
  --o-taxonomy mcrA_ncbi_taxonomy.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fetching mrtA sequences from NCBI ..."
qiime rescript get-ncbi-data \
  --p-query "${Q['mrtA']}" \
  --p-ranks domain kingdom phylum class order family genus species \
  --verbose \
  --p-n-jobs 8 \
  --o-sequences mrtA_ncbi_seqs.qza \
  --o-taxonomy mrtA_ncbi_taxonomy.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Tagging mcrA and mrtA feature IDs before merging ..."
for gene in mcrA mrtA; do
    echo "  Tagging ${gene} ..."

    qiime tools export \
      --input-path "${gene}_ncbi_seqs.qza" \
      --output-path "${gene}_exported"

    qiime tools export \
      --input-path "${gene}_ncbi_taxonomy.qza" \
      --output-path "${gene}_tax_exported"

    # Append _mcrA / _mrtA to sequence headers and taxonomy feature IDs
    sed -E "s/^>(.+)$/>\\1_${gene}/" "${gene}_exported/dna-sequences.fasta" > "${gene}_tagged.fasta"

    awk -F'\t' -v suffix="${gene}" 'BEGIN {OFS="\t"} NR==1 {print; next} {$1=$1"_"suffix; print}' \
      "${gene}_tax_exported/taxonomy.tsv" > "${gene}_tagged_taxonomy.tsv"

    qiime tools import \
      --input-path "${gene}_tagged.fasta" \
      --output-path "${gene}_ncbi_seqs_tagged.qza" \
      --type 'FeatureData[Sequence]'

    qiime tools import \
      --input-path "${gene}_tagged_taxonomy.tsv" \
      --output-path "${gene}_ncbi_taxonomy_tagged.qza" \
      --type 'FeatureData[Taxonomy]'
done

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging tagged mcrA + mrtA ..."
qiime feature-table merge-seqs \
  --i-data mcrA_ncbi_seqs_tagged.qza \
  --i-data mrtA_ncbi_seqs_tagged.qza \
  --o-merged-data ncbi_mcrA_mrtA_seqs.qza

qiime feature-table merge-taxa \
  --i-data mcrA_ncbi_taxonomy_tagged.qza \
  --i-data mrtA_ncbi_taxonomy_tagged.qza \
  --o-merged-data ncbi_mcrA_mrtA_taxonomy.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting merged FASTA and taxonomy ..."
qiime tools export \
  --input-path ncbi_mcrA_mrtA_seqs.qza \
  --output-path ncbi_exported

qiime tools export \
  --input-path ncbi_mcrA_mrtA_taxonomy.qza \
  --output-path ncbi_exported

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/01_ncbi_rescript"
mkdir -p "${OUT_STORAGE}"
cp -r "${WORKDIR}/"*.qza "${OUT_STORAGE}/"
cp -r "${WORKDIR}/ncbi_exported" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 01 complete. Results in ${OUT_STORAGE}"
