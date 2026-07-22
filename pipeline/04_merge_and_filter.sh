#!/bin/bash
#PBS -N merge_filter
#PBS -l select=1:ncpus=8:mem=32gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 04_merge_and_filter.sh
# Merge pre-fixed RESCRIPt + GeneScoop taxonomies, import to QIIME2, cull, filter length. Requires fix_rescript_taxonomy.py and fix_genescoop_taxonomy.py to have been run beforehand.
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

# RESCRIPt outputs
RESCRIPT_SEQ="${STORAGE}/01_ncbi_rescript/ncbi_mcrA_mrtA_seqs.qza"
RESCRIPT_TAX="${STORAGE}/01_ncbi_rescript/ncbi_mcrA_mrtA_taxonomy_fixed.tsv"
# GeneScoop outputs for mcrA and mrtA
GENESCOOP_SEQ_MCRA="${STORAGE}/03_genescoop/genescoop_clean_mcrA.fasta"
GENESCOOP_TAX_MCRA="${STORAGE}/03_genescoop/genescoop_taxonomy_mcrA.tsv"
GENESCOOP_SEQ_MRTA="${STORAGE}/03_genescoop/genescoop_clean_mrtA.fasta"
GENESCOOP_TAX_MRTA="${STORAGE}/03_genescoop/genescoop_taxonomy_mrtA.tsv"

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${RESCRIPT_SEQ}" .
cp "${RESCRIPT_TAX}" .
cp "${GENESCOOP_SEQ_MCRA}" .
cp "${GENESCOOP_TAX_MCRA}" .
cp "${GENESCOOP_SEQ_MRTA}" .
cp "${GENESCOOP_TAX_MRTA}" .

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building unified GeneScoop files ..."
# Concatenate mcrA and mrtA clean FASTAs into a single unified file
cat genescoop_clean_mcrA.fasta genescoop_clean_mrtA.fasta > genescoop_clean.fasta

# Build unified taxonomy with one header line and all data rows from both files
head -1 genescoop_taxonomy_mcrA.tsv > genescoop_taxonomy.tsv
tail -n +2 genescoop_taxonomy_mcrA.tsv >> genescoop_taxonomy.tsv
tail -n +2 genescoop_taxonomy_mrtA.tsv >> genescoop_taxonomy.tsv

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting RESCRIPt FASTA ..."
qiime tools export \
  --input-path ncbi_mcrA_mrtA_seqs.qza \
  --output-path rescript_fasta

RESCRIPT_FASTA="rescript_fasta/dna-sequences.fasta"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Using pre-generated GeneScoop cleaned files ..."

# Merge FASTAs
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging sequences ..."
cat "${RESCRIPT_FASTA}" genescoop_clean.fasta > all_mcrA_mrtA_raw.fasta

# Build unified taxonomy from pre-fixed RESCRIPt and GeneScoop taxonomies.
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building unified taxonomy ..."
cat ncbi_mcrA_mrtA_taxonomy_fixed.tsv > all_taxonomy.tsv
tail -n +2 genescoop_taxonomy.tsv >> all_taxonomy.tsv

# QIIME2 import
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Importing to QIIME2 ..."
qiime tools import \
  --input-path all_mcrA_mrtA_raw.fasta \
  --output-path all_seqs.qza \
  --type 'FeatureData[Sequence]'

qiime tools import \
  --input-path all_taxonomy.tsv \
  --output-path all_taxonomy.qza \
  --type 'FeatureData[Taxonomy]'

# Cull and filter
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Culling sequences (degenerates <=1, homopolymer <=12) ..."
qiime rescript cull-seqs \
  --i-sequences all_seqs.qza \
  --p-num-degenerates 1 \
  --p-homopolymer-length 12 \
  --p-n-jobs 8 \
  --o-clean-sequences seqs_clean.qza

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Filtering by length (172 and 2000 bp) ..."
qiime rescript filter-seqs-length \
  --i-sequences seqs_clean.qza \
  --p-global-min 172 \
  --p-global-max 2000 \
  --p-threads 8 \
  --o-filtered-seqs seqs_filt.qza \
  --o-discarded-seqs seqs_discarded.qza

# Export filtered FASTA for MAFFT
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Exporting filtered FASTA ..."
qiime tools export \
  --input-path seqs_filt.qza \
  --output-path seqs_filt_exported

cp seqs_filt_exported/dna-sequences.fasta seqs_for_alignment.fasta

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/04_merged_filtered"
mkdir -p "${OUT_STORAGE}"

# Only copy newly generated outputs (avoid copying input QZAs we brought in)
cp "${WORKDIR}/all_seqs.qza"          "${OUT_STORAGE}/" || true
cp "${WORKDIR}/all_taxonomy.qza"      "${OUT_STORAGE}/" || true
cp "${WORKDIR}/seqs_clean.qza"        "${OUT_STORAGE}/" || true
cp "${WORKDIR}/seqs_filt.qza"         "${OUT_STORAGE}/" || true
cp "${WORKDIR}/seqs_discarded.qza"    "${OUT_STORAGE}/" || true
cp -r "${WORKDIR}/seqs_filt_exported" "${OUT_STORAGE}/" || true
cp "${WORKDIR}/seqs_for_alignment.fasta" "${OUT_STORAGE}/" || true
cp "${WORKDIR}/all_taxonomy.tsv"         "${OUT_STORAGE}/" || true
cp "${WORKDIR}/genescoop_taxonomy.tsv"   "${OUT_STORAGE}/" || true
cp "${WORKDIR}/genescoop_clean.fasta"    "${OUT_STORAGE}/" || true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 04 complete. Results in ${OUT_STORAGE}"
