#!/bin/bash
#PBS -N mafft_align
#PBS -l select=1:ncpus=16:mem=64gb:scratch_ssd=500gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 06_mafft_align.sh
# 1) Align dereplicated sequences with MAFFT --auto
# 2) Add guide sequences using --addfragments --multipair --keeplength --mapout
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

SEQ_IN="${STORAGE}/05_dereplicated/seqs_for_alignment.fasta"
GUIDE_SEQS="${STORAGE}/guide_seqs.fas"
MAFFT_THREADS=16

module add mafft
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${SEQ_IN}" seqs_for_alignment.fasta
cp "${GUIDE_SEQS}" guide_seqs.fas

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running MAFFT --auto ..."
mafft \
  --auto \
  --thread "${MAFFT_THREADS}" \
  --reorder \
  seqs_for_alignment.fasta > alignment_raw.fasta

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Adding guide sequences ..."
mafft \
  --multipair \
  --addfragments guide_seqs.fas \
  --keeplength \
  --mapout \
  --thread "${MAFFT_THREADS}" \
  alignment_raw.fasta > alignment_with_guides.fasta

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/06_alignment"
mkdir -p "${OUT_STORAGE}"
cp "${WORKDIR}/alignment_raw.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_with_guides.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_raw.fasta.map" "${OUT_STORAGE}/" 2>/dev/null || true

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 06 complete. Results in ${OUT_STORAGE}"
