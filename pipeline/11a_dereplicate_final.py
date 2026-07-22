from Bio import SeqIO

# Read taxonomy IDs
with open('taxonomy_final_cleaned.tsv') as f:
    header = f.readline()
    tax_ids = {line.split('\t')[0].strip() for line in f}
print(f'  Taxonomy IDs: {len(tax_ids)}')

# Filter alignment and strip gaps for QIIME2 import.
# Also drop guide/primer sequences, which are not in the taxonomy.
ungapped_recs = []
n_filtered = 0
for rec in SeqIO.parse('alignment_final_cleaned.fasta', 'fasta'):
    if rec.id.startswith(('mcrA_', 'mlas_', 'guide_')):
        n_filtered += 1
        continue
    if rec.id not in tax_ids:
        n_filtered += 1
        continue
    seq = str(rec.seq).upper().replace('-', '').replace('.', '')
    if not seq:
        n_filtered += 1
        continue
    ungapped_recs.append(rec.__class__(rec.seq.__class__(seq), id=rec.id, description=''))

SeqIO.write(ungapped_recs, 'seqs_ungapped.fasta', 'fasta')
print(f'  Ungapped sequences for dereplication: {len(ungapped_recs)}')
print(f'  Filtered out (guides/primers/missing-tax/empty): {n_filtered}')

# Prepare taxonomy for QIIME2 (column must be named 'Taxon')
with open('taxonomy_final_cleaned.tsv') as f, open('taxonomy_qiime.tsv', 'w') as fout:
    h = f.readline()
    fout.write(h.replace('new_tax', 'Taxon'))
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if parts[0] in tax_ids:
            fout.write(line)
