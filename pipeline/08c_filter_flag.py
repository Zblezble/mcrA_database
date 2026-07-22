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
