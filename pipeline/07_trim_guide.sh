#!/bin/bash
#PBS -N trim_guide
#PBS -l select=1:ncpus=4:mem=16gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 07_trim_guide.sh
# Trim the MAFFT alignment to the span of the reference guide
# (guide_Y00158.1_Methanosarcina_barkeri) plus optional radius.
# Also computes per-sequence pairwise identity vs the guide within the span.
#
# Output:
#   alignment_trimmed.fasta  — alignment trimmed to guide span
#   guide_identity_scores.tsv — identity vs P07962 for each sequence
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

SEQ_ALIGNED="${STORAGE}/06_alignment/alignment_with_guides.fasta"
GUIDE_NAME="guide_Y00158.1_Methanosarcina_barkeri"
GUIDE_RADIUS=10            # bp to keep outside guide span (0 = exact trim)

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying alignment to scratch ..."
cp "${SEQ_ALIGNED}" .

python3 - "${GUIDE_NAME}" "${GUIDE_RADIUS}" <<'PY'
import sys
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

GUIDE_NAME = sys.argv[1]
GUIDE_RADIUS = int(sys.argv[2])

align = AlignIO.read('alignment_with_guides.fasta', 'fasta')

# Find guide index
guide_idx = None
for i, rec in enumerate(align):
    if GUIDE_NAME in rec.id:
        guide_idx = i
        break

if guide_idx is None:
    raise RuntimeError(f'{GUIDE_NAME} not found in alignment. Cannot trim.')

guide_seq = str(align[guide_idx].seq)

# Find first and last non-gap positions in guide
first = next((i for i, ch in enumerate(guide_seq) if ch not in '-.'), 0)
last = next((i for i in range(len(guide_seq) - 1, -1, -1) if guide_seq[i] not in '-.'), len(guide_seq) - 1)

# Apply radius
radius = GUIDE_RADIUS
start = max(0, first - radius)
end = min(len(guide_seq) - 1, last + radius)

print(f'Guide span: columns {first}-{last} ({last - first + 1} bp)')
print(f'Trim range: columns {start}-{end} ({end - start + 1} bp, radius ±{radius})')

# Trim alignment
span = end - start + 1
trimmed = MultipleSeqAlignment([])
for rec in align:
    rec.seq = rec.seq[start:end + 1]
    trimmed.append(rec)

AlignIO.write(trimmed, 'alignment_trimmed.fasta', 'fasta')

# ---- Compute identity scores vs guide ----
guide_trimmed = str(trimmed[guide_idx].seq)
n_cols = len(guide_trimmed)

with open('guide_identity_scores.tsv', 'w') as f:
    f.write('seq_id\tguide_match\tguide_mismatch\tguide_gap\tidentity_pct\tnogap_identity_pct\tungapped_len\n')
    for rec in trimmed:
        if GUIDE_NAME in rec.id:
            continue
        seq = str(rec.seq)
        match = mismatch = gap = 0
        guide_cov = 0
        for g, q in zip(guide_trimmed, seq):
            if g in '-.':
                continue
            guide_cov += 1
            if q in '-.':
                gap += 1
            elif q.upper() == g.upper():
                match += 1
            else:
                mismatch += 1
        identity = (match / guide_cov * 100) if guide_cov > 0 else 0.0
        ungapped = guide_cov - gap
        # no-gap identity: (1 - mismatch / ungapped_len) * 100
        nogap_identity = ((1 - mismatch / ungapped) * 100) if ungapped > 0 else 0.0
        f.write(f'{rec.id}\t{match}\t{mismatch}\t{gap}\t{identity:.2f}\t{nogap_identity:.2f}\t{ungapped}\n')

print(f'Scored {len(trimmed) - 1} sequences against {GUIDE_NAME}')
PY

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/07_trimmed"
mkdir -p "${OUT_STORAGE}"
cp "${WORKDIR}/alignment_trimmed.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/guide_identity_scores.tsv" "${OUT_STORAGE}/"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 07 complete. Results in ${OUT_STORAGE}"
