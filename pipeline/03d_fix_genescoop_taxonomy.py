#!/usr/bin/env python3
"""
03d_fix_genescoop_taxonomy.py
Extract clean FASTA + QIIME2 taxonomy TSV from GeneScoop mrtA headers.

Input:
    03_genescoop/genescoop_mrtA.fasta

Outputs:
    03_genescoop/genescoop_clean_mrtA.fasta
    03_genescoop/genescoop_taxonomy_mrtA.tsv
"""

import re

INPUT_FASTA = '03_genescoop/genescoop_mrtA.fasta'
OUTPUT_FASTA = '03_genescoop/genescoop_clean_mrtA.fasta'
OUTPUT_TAX = '03_genescoop/genescoop_taxonomy_mrtA.tsv'


def extract_species(definition):
    """Pull the first binomial name out of the GenBank description."""
    blacklist = {
        'strain', 'sp', 'sp.', 'contig', 'scaffold', 'complete',
        'partial', 'sequence', 'gene', 'product', 'note', 'genome',
        'shotgun', 'chromosome', 'plasmid', 'draft', 'clone', 'cell',
        'line', 'assembly', 'type', 'and', 'the'
    }
    for m in re.finditer(r'\b(Candidatus\s+)?([A-Z][a-z]+)\s+([a-z][a-z]+)\b', definition):
        cand, genus, epithet = m.group(1), m.group(2), m.group(3)
        if epithet.lower() in blacklist or any(ch.isdigit() for ch in epithet):
            continue
        return f'{cand.strip()} {genus} {epithet}' if cand else f'{genus} {epithet}'
    return None


with open(INPUT_FASTA) as fin, open(OUTPUT_FASTA, 'w') as fseq, open(OUTPUT_TAX, 'w') as ftax:
    ftax.write('Feature ID\tTaxon\n')
    for line in fin:
        if not line.startswith('>'):
            fseq.write(line)
            continue

        header = line[1:].rstrip('\n')
        sp = header.find(' ')
        if sp == -1:
            fseq.write('>' + header + '\n')
            continue
        seq_id = header[:sp]
        rest = header[sp + 1:]

        # split by `, ` into blocks
        blocks = [b.strip() for b in rest.split(', ')]
        block1 = blocks[0]          # match_type + definition
        block3 = blocks[-1] if len(blocks) > 1 else ''

        # species from definition (block 1, after match_type)
        mt_end = block1.find(' ')
        definition = block1[mt_end + 1:] if mt_end != -1 else block1
        species = extract_species(definition)

        # taxonomy from block 3
        raw_tokens = [t.strip() for t in block3.split(';') if t.strip()]
        # drop non-standard ranks (anything containing ' group')
        clean = [t for t in raw_tokens if ' group' not in t.lower()]

        # assign ranks left-to-right
        rank_names = ['domain', 'kingdom', 'phylum', 'class',
                      'order', 'family', 'genus', 'species']
        assigned = {r: '' for r in rank_names}
        for i, token in enumerate(clean):
            if i < len(rank_names) - 1:
                assigned[rank_names[i]] = token

        assigned['domain'] = 'Archaea'

        # If no species was extracted, fall back to "<genus> sp."
        if species:
            assigned['species'] = species
        elif assigned['genus']:
            assigned['species'] = f"{assigned['genus']} sp."
        else:
            assigned['species'] = ''

        tax_str = '; '.join(f'{r[0]}__{assigned[r]}' for r in rank_names)
        fseq.write('>' + seq_id + '\n')
        ftax.write(seq_id + '\t' + tax_str + '\n')
