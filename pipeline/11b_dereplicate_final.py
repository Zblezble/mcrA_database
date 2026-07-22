from Bio import SeqIO

# Get representative IDs retained by RESCRIPt
with open('taxonomy_derep_final.tsv') as f:
    next(f)
    kept_ids = {line.split('\t')[0].strip() for line in f}
print(f'  Dereplicated IDs: {len(kept_ids)}')

# Subset original aligned FASTA (drop guides/primers and duplicates)
aligned_kept = []
for rec in SeqIO.parse('alignment_final_cleaned.fasta', 'fasta'):
    if rec.id.startswith(('mcrA_', 'mlas_', 'guide_')):
        continue
    if rec.id in kept_ids:
        aligned_kept.append(rec)

SeqIO.write(aligned_kept, 'alignment_final_derep.fasta', 'fasta')
print(f'  Aligned sequences kept: {len(aligned_kept)}')
