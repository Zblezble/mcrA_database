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
