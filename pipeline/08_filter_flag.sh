#!/bin/bash
#PBS -N filter_flag
#PBS -l select=1:ncpus=8:mem=32gb:scratch_ssd=10gb
#PBS -l walltime=2:00:00
#PBS -m ae
#PBS -M mail@mail.mail

# ============================================================================
# 08_filter_flag.sh
# Post-trim QC: guide identity, short-region length, tree outliers, etc.
# ============================================================================

set -euo pipefail

# ---- CONFIGURATION ---------------------------------------------------------
STORAGE="/path/to/storage"
WORKDIR="${SCRATCHDIR}/mcrA_work"

SEQ_ALIGNED="${STORAGE}/07_trimmed/alignment_trimmed.fasta"
IDENTITY_SCORES="${STORAGE}/07_trimmed/guide_identity_scores.tsv"
MIN_GUIDE_IDENTITY=0.46
OUTLIER_MULTIPLIER=3.0

module add mambaforge
conda activate /path/to/.conda/envs/qiime2-amplicon-2024.2
module add fasttree
# ----------------------------------------------------------------------------

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying inputs to scratch ..."
cp "${SEQ_ALIGNED}" .
cp "${IDENTITY_SCORES}" .

# ---- 1. Flag low guide identity ---------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Flagging low guide identity ..."
awk -v min="${MIN_GUIDE_IDENTITY}" '
NR > 1 {
    id = $6 / 100.0;
    if (id < min) print $1 "\t" $6
}' guide_identity_scores.tsv | sort -t$'\t' -k2 -n > guide_low_identity_flag.txt
LOW=$(wc -l < guide_low_identity_flag.txt)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${LOW} sequences flagged (guide identity < ${MIN_GUIDE_IDENTITY}%)"

# ---- 2. Short-read region length filter (using mcrA_forward/reverse guides) --
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Checking amplicon region lengths ..."

python3 - <<'PY'
from Bio import SeqIO

recs = list(SeqIO.parse('alignment_trimmed.fasta', 'fasta'))
ids = [r.id for r in recs]
seqs = [list(str(r.seq)) for r in recs]
n_cols = len(seqs[0])

def get_span(indices):
    """Return (start, end) covering the first/last non-gap across given rows."""
    positions = []
    for idx in indices:
        for col, ch in enumerate(seqs[idx]):
            if ch not in '-.':
                positions.append(col)
    return (min(positions), max(positions)) if positions else (0, n_cols - 1)

# ---- Short-read region = span of guide_short_* only ----
# Primers are NOT included because newer sequencing platforms strip them.
short_indices = [i for i, name in enumerate(ids) if name.startswith('guide_short_')]

short_start, short_end = get_span(short_indices)
print(f'Short-read span: columns {short_start}-{short_end} ({short_end - short_start + 1} bp)')

short_flag = []
for qi, sid in enumerate(ids):
    if sid.startswith('mcrA_') or sid.startswith('guide_'):
        continue
    length = sum(1 for col in range(short_start, short_end + 1) if seqs[qi][col] not in '-.')
    if length < 172:
        short_flag.append((sid, length))

with open('guide_short_flag.txt', 'w') as f:
    f.write('seq_id\tungapped_length_in_short_region\n')
    for sid, L in short_flag:
        f.write(f'{sid}\t{L}\n')
print(f'Short-region (<172 bp) flagged: {len(short_flag)}')

# ---- Long-read region = span of guide_long_* ----
long_indices = [i for i, name in enumerate(ids) if name.startswith('guide_long_')]
long_start, long_end = get_span(long_indices)
print(f'Long-read span: columns {long_start}-{long_end} ({long_end - long_start + 1} bp)')

long_flag = []
for qi, sid in enumerate(ids):
    if sid.startswith('mcrA_') or sid.startswith('guide_'):
        continue
    length = sum(1 for col in range(long_start, long_end + 1) if seqs[qi][col] not in '-.')
    if length < 392:
        long_flag.append((sid, length))

with open('guide_long_flag.txt', 'w') as f:
    f.write('seq_id\tungapped_length_in_long_region\n')
    for sid, L in long_flag:
        f.write(f'{sid}\t{L}\n')
print(f'Long-region (<392 bp) flagged: {len(long_flag)}')

# ---- Write out trimmed alignments for sanity check ----
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

# short trim
short_align = MultipleSeqAlignment([])
for rec in recs:
    rec.seq = rec.seq.__class__(''.join(rec.seq[i] for i in range(short_start, short_end + 1)))
    short_align.append(rec)
AlignIO.write(short_align, 'alignment_trimmed_short.fasta', 'fasta')
print(f'Wrote alignment_trimmed_short.fasta ({short_end - short_start + 1} columns)')

# long trim
long_align = MultipleSeqAlignment([])
# reload because we mutated recs above
recs = list(SeqIO.parse('alignment_trimmed.fasta', 'fasta'))
for rec in recs:
    rec.seq = rec.seq.__class__(''.join(rec.seq[i] for i in range(long_start, long_end + 1)))
    long_align.append(rec)
AlignIO.write(long_align, 'alignment_trimmed_long.fasta', 'fasta')
print(f'Wrote alignment_trimmed_long.fasta ({long_end - long_start + 1} columns)')
PY

# ---- 3. Build tree + flag outliers ------------------------------------------
# Trim gappy columns from trimmed alignment for tree speed
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Trimming gappy columns for tree ..."
python3 - <<'PY'
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment

align = AlignIO.read('alignment_trimmed.fasta', 'fasta')
keep_cols = [i for i in range(align.get_alignment_length())
             if any(c not in '-.Nn?' for c in align[:, i])]
new_align = MultipleSeqAlignment([])
for rec in align:
    new_seq = ''.join(rec.seq[i] for i in keep_cols)
    rec.seq = rec.seq.__class__(new_seq)
    new_align.append(rec)
AlignIO.write(new_align, 'alignment_gappy_trimmed.fasta', 'fasta')
print(f'Trimmed from {align.get_alignment_length()} to {len(keep_cols)} columns')
PY

if command -v FastTree &> /dev/null; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building tree with FastTree ..."
    FastTree -gtr -nt alignment_gappy_trimmed.fasta > tree.nwk 2> fasttree.log
else
    echo "WARNING: FastTree not found. Skipping tree."
    touch tree.nwk
fi

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Flagging long-branch outliers ..."
python3 - <<'PY'
import statistics, os
from Bio import Phylo

if not os.path.exists('tree.nwk') or os.path.getsize('tree.nwk') == 0:
    with open('tree_outliers_flag.txt', 'w') as f:
        f.write('seq_id\tbranch_length\tcutoff\n')
    exit()

tree = Phylo.read('tree.nwk', 'newick')
term_bl = {tip.name: tree.distance(tip.name) for tip in tree.get_terminals()}
values = list(term_bl.values())
med = statistics.median(values)
mad = statistics.median([abs(v - med) for v in values])
cutoff = med + 3.0 * mad
if cutoff <= med:
    cutoff = med * 2.0

outliers = []
for name, bl in term_bl.items():
    if bl > cutoff:
        outliers.append((name, bl, cutoff))
outliers.sort(key=lambda x: x[1], reverse=True)

with open('tree_outliers_flag.txt', 'w') as f:
    f.write('seq_id\tbranch_length\tcutoff (med + 3*MAD)\n')
    for name, bl, co in outliers:
        f.write(f'{name}\t{bl:.6f}\t{co:.6f}\n')
print(f'Tree outliers flagged: {len(outliers)}')
PY

# ---- 4. Merge all flags -----------------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Merging flag lists ..."

# Start with guide low identity
cp guide_low_identity_flag.txt outliers_flag.txt 2>/dev/null || touch outliers_flag.txt

# Append short-region
if [[ -s guide_short_flag.txt ]]; then
    echo "" >> outliers_flag.txt
    echo "# Short-region <172 bp (guide-based)" >> outliers_flag.txt
    tail -n +2 guide_short_flag.txt >> outliers_flag.txt
fi

# Append long-region
if [[ -s guide_long_flag.txt ]]; then
    echo "" >> outliers_flag.txt
    echo "# Long-region <392 bp (guide-based, 80% of ~490 bp)" >> outliers_flag.txt
    tail -n +2 guide_long_flag.txt >> outliers_flag.txt
fi

# Append tree outliers
if [[ -s tree_outliers_flag.txt ]]; then
    echo "" >> outliers_flag.txt
    echo "# Long-branch outliers (tree-based)" >> outliers_flag.txt
    tail -n +2 tree_outliers_flag.txt >> outliers_flag.txt
fi

# ---- Copy results back to storage -------------------------------------------
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Copying results to storage ..."
OUT_STORAGE="${STORAGE}/08_qc"
mkdir -p "${OUT_STORAGE}"
cp "${WORKDIR}/guide_low_identity_flag.txt" "${OUT_STORAGE}/"
cp "${WORKDIR}/guide_short_flag.txt" "${OUT_STORAGE}/"
cp "${WORKDIR}/guide_long_flag.txt" "${OUT_STORAGE}/"
cp "${WORKDIR}/tree_outliers_flag.txt" "${OUT_STORAGE}/"
cp "${WORKDIR}/outliers_flag.txt" "${OUT_STORAGE}/"
cp "${WORKDIR}/tree.nwk" "${OUT_STORAGE}/"
cp "${WORKDIR}/fasttree.log" "${OUT_STORAGE}/" 2>/dev/null || true
cp "${WORKDIR}/alignment_trimmed_short.fasta" "${OUT_STORAGE}/"
cp "${WORKDIR}/alignment_trimmed_long.fasta" "${OUT_STORAGE}/"

TOT=$(grep -vc '^#' outliers_flag.txt || wc -l < outliers_flag.txt)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Step 08 complete. Results in ${OUT_STORAGE}"
echo ""
echo "TOTAL FLAGGED SEQUENCES FOR MANUAL REVIEW: ${TOT}"
echo "  - guide_low_identity_flag.txt  (identity < ${MIN_GUIDE_IDENTITY}% vs guide)"
echo "  - guide_short_flag.txt         (<172 bp in short-read amplicon region, 80% of ~215 bp)"
echo "  - guide_long_flag.txt          (<392 bp in long-read amplicon region, 80% of ~490 bp)"
echo "  - tree_outliers_flag.txt       (unexpectedly long branches)"
echo ""
echo "SANITY-CHECK ALIGNMENTS:"
echo "  - alignment_trimmed_short.fasta  (trimmed to guide_short span)"
echo "  - alignment_trimmed_long.fasta   (trimmed to guide_long span)"
echo ""
echo "NEXT STEPS (manual):"
echo "  1. Inspect tree (tree.nwk) and flagged IDs"
echo "  2. Remove confirmed bad sequences by ID"
echo "  3. Trim to amplicon boundaries (LRD/SRD)"
echo "  4. Dereplicate each trimmed dataset separately"
echo "  5. Reassign mrtA misannotations"
