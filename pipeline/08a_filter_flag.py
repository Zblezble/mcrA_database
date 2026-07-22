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
