#!/bin/bash
#PBS -N dereplicate
#PBS -l select=1:ncpus=2:mem=16gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 05_dereplicate.sh
# Dereplicate merged (RESCRIPt + GeneScoop) sequences with RESCRIPt (uniq).
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

SEQ_IN="${STORAGE}/04_merged_filtered/seqs_filt.qza"
TAXONOMY="${STORAGE}/04_merged_filtered/all_taxonomy.qza"

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${SEQ_IN}" .
cp "${TAXONOMY}" .

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Dereplicating (mode=uniq) ..."
qiime rescript dereplicate --verbose \
  --i-sequences seqs_filt.qza \
  --i-taxa all_taxonomy.qza \
  --p-mode uniq \
  --p-threads 1 \
  --p-derep-prefix \
  --o-dereplicated-sequences seqs_derep.qza \
  --o-dereplicated-taxa taxonomy_derep.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting dereplicated FASTA ..."
qiime tools export \
  --input-path seqs_derep.qza \
  --output-path derep_exported

cp derep_exported/dna-sequences.fasta seqs_for_alignment.fasta

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting dereplicated taxonomy ..."
qiime tools export \
  --input-path taxonomy_derep.qza \
  --output-path taxonomy_derep_exported

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/05_dereplicated"
mkdir -p "${OUT_STORAGE}"
cp "${WORKDIR}/seqs_derep.qza" "${OUT_STORAGE}/"
cp "${WORKDIR}/taxonomy_derep.qza" "${OUT_STORAGE}/"
cp -r "${WORKDIR}/derep_exported" "${OUT_STORAGE}/"
cp -r "${WORKDIR}/taxonomy_derep_exported" "${OUT_STORAGE}/"
cp "${WORKDIR}/seqs_for_alignment.fasta" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 05 complete. Results in ${OUT_STORAGE}"
