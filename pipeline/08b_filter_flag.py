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
