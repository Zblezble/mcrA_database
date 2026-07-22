#!/bin/bash
#PBS -N remove_flagged
#PBS -l select=1:ncpus=4:mem=16gb:scratch_ssd=10gb
#PBS -l walltime=1:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 09_remove_flagged.sh
# Three-step removal from the trimmed alignment:
#   1. Remove low-quality sequences listed in remove_quality.tsv
#      -> outputs alignment_after_quality.fasta
#   2. Remove short-length sequences listed in remove_length.tsv
#      -> outputs alignment_after_length.fasta
#   3. Remove sequences with missing/problematic taxonomy listed in remove_taxonomy.tsv
#      -> outputs alignment_final.fasta
#
# All steps also collapse all-gap columns after removal.
# ============================================================================

set -euo pipefail

STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

ALN_IN="${STORAGE}/07_trimmed/alignment_trimmed.fasta"
QUALITY_IDS="${STORAGE}/08_qc/remove_quality.tsv"
LENGTH_IDS="${STORAGE}/08_qc/remove_length.tsv"
TAXONOMY_IDS="${STORAGE}/08_qc/remove_taxonomy.tsv"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

cp "${ALN_IN}" .

# ---- Step 1: remove quality-flagged sequences -------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 1: removing quality-flagged sequences ..."
if [[ -f "${QUALITY_IDS}" ]]; then
    QUAL_COUNT=$(tail -n +2 "${QUALITY_IDS}" 2>/dev/null | wc -l || echo "0")
else
    QUAL_COUNT=0
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${QUAL_COUNT} IDs in remove_quality.tsv"

if [[ "${QUAL_COUNT}" -gt 0 ]]; then
    python3 - "${QUALITY_IDS}" <<'PY'
import sys
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

tsv_path = sys.argv[1]
remove_ids = set()
with open(tsv_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if parts:
            remove_ids.add(parts[0])

records = [rec for rec in SeqIO.parse('alignment_trimmed.fasta', 'fasta')
           if rec.id not in remove_ids]
print(f'After quality removal: {len(records)} sequences')

# collapse all-gap columns
align = MultipleSeqAlignment(records)
keep_cols = [i for i in range(align.get_alignment_length())
             if any(c not in '-.' for c in align[:, i])]
clean = MultipleSeqAlignment([])
for rec in align:
    rec.seq = rec.seq.__class__(''.join(rec.seq[i] for i in keep_cols))
    clean.append(rec)
AlignIO.write(clean, 'alignment_after_quality.fasta', 'fasta')
print(f'Collapsed {align.get_alignment_length()} -> {len(keep_cols)} columns')
PY
else
    cp alignment_trimmed.fasta alignment_after_quality.fasta
fi

# ---- Step 2: remove length-flagged sequences --------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 2: removing length-flagged sequences ..."
if [[ -f "${LENGTH_IDS}" ]]; then
    LEN_COUNT=$(tail -n +2 "${LENGTH_IDS}" 2>/dev/null | wc -l || echo "0")
else
    LEN_COUNT=0
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${LEN_COUNT} IDs in remove_length.tsv"

if [[ "${LEN_COUNT}" -gt 0 ]]; then
    python3 - "${LENGTH_IDS}" <<'PY'
import sys
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

tsv_path = sys.argv[1]
remove_ids = set()
with open(tsv_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if parts:
            remove_ids.add(parts[0])

records = [rec for rec in SeqIO.parse('alignment_after_quality.fasta', 'fasta')
           if rec.id not in remove_ids]
print(f'After length removal: {len(records)} sequences')

# collapse all-gap columns again
align = MultipleSeqAlignment(records)
keep_cols = [i for i in range(align.get_alignment_length())
             if any(c not in '-.' for c in align[:, i])]
clean = MultipleSeqAlignment([])
for rec in align:
    rec.seq = rec.seq.__class__(''.join(rec.seq[i] for i in keep_cols))
    clean.append(rec)
AlignIO.write(clean, 'alignment_after_length.fasta', 'fasta')
print(f'Collapsed {align.get_alignment_length()} -> {len(keep_cols)} columns')
PY
else
    cp alignment_after_quality.fasta alignment_after_length.fasta
fi

# ---- Step 3: remove taxonomy-flagged sequences ------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 3: removing taxonomy-flagged sequences ..."
if [[ -f "${TAXONOMY_IDS}" ]]; then
    TAX_COUNT=$(tail -n +2 "${TAXONOMY_IDS}" 2>/dev/null | wc -l || echo "0")
else
    TAX_COUNT=0
fi
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${TAX_COUNT} IDs in remove_taxonomy.tsv"

if [[ "${TAX_COUNT}" -gt 0 ]]; then
    python3 - "${TAXONOMY_IDS}" <<'PY'
import sys
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

tsv_path = sys.argv[1]
remove_ids = set()
with open(tsv_path) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if parts:
            remove_ids.add(parts[0])

records = [rec for rec in SeqIO.parse('alignment_after_length.fasta', 'fasta')
           if rec.id not in remove_ids]
print(f'After taxonomy removal: {len(records)} sequences')

# collapse all-gap columns again
align = MultipleSeqAlignment(records)
keep_cols = [i for i in range(align.get_alignment_length())
             if any(c not in '-.' for c in align[:, i])]
clean = MultipleSeqAlignment([])
for rec in align:
    rec.seq = rec.seq.__class__(''.join(rec.seq[i] for i in keep_cols))
    clean.append(rec)
AlignIO.write(clean, 'alignment_final.fasta', 'fasta')
print(f'Collapsed {align.get_alignment_length()} -> {len(keep_cols)} columns')
PY
else
    cp alignment_after_length.fasta alignment_final.fasta
fi

# ---- Copy results back to storage -------------------------------------------
OUT_STORAGE="${STORAGE}/09_clean"
mkdir -p "${OUT_STORAGE}"
cp "${WORKDIR}/alignment_after_quality.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_after_length.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_final.fasta" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 09 complete."
echo "  alignment_after_quality.fasta  — after removing ${QUAL_COUNT} quality-flagged"
echo "  alignment_after_length.fasta   — after removing ${LEN_COUNT} length-flagged"
echo "  alignment_final.fasta          — after removing ${TAX_COUNT} taxonomy-flagged"
echo "  Results in ${OUT_STORAGE}"
