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
