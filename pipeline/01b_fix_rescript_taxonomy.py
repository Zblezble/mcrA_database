#!/usr/bin/env python
# coding: utf-8

"""
01b_fix_rescript_taxonomy.py
Fix a RESCRIPt-exported QIIME2 taxonomy TSV so it can be merged cleanly with
the GeneScoop taxonomy.

Changes applied:
  - Remove superkingdom rank(s) such as sk__Archaea.
  - Convert short species names to long format: s__horonobensis ->
    s__Methanosarcina horonobensis.
  - Missing species becomes s__<genus> sp.
  - Strip strain/accession suffixes from species (e.g. DSM 5219, Kol 5).

The input is the TSV exported by step 01 (``01_get_ncbi_mcrA_mrtA.sh``).
"""

import re

DEFAULT_INPUT_TSV = '01_ncbi_rescript/ncbi_exported/taxonomy.tsv'
DEFAULT_OUTPUT_TSV = '01_ncbi_rescript/ncbi_mcrA_mrtA_taxonomy_fixed.tsv'

# Tokens that, when they appear after the species epithet, mark the start of a strain/accession/designation suffix that should be discarded.
STRAIN_PREFIXES = {
    'atcc', 'bcm', 'bp', 'ccug', 'cell', 'cfbp', 'chromosome', 'cip', 'clone',
    'complete', 'contig', 'ctg', 'dna', 'draft', 'dsm', 'genome', 'gp', 'jf',
    'jcm', 'ka', 'kin', 'kctc', 'kol', 'kole', 'lam', 'lmg', 'mag', 'mc',
    'mc-s', 'mcs', 'mo', 'nbi', 'nbrc', 'ncimb', 'nctc', 'no', 'nobi', 'os',
    'partial', 'plasmid', 'sael', 'sarl', 'scaffold', 'sequence', 'shotgun',
    'smsp', 'sp', 'spp', 'strain', 'tm', 'type', 'wgs', 'whole',
}

def looks_like_strain(word):
    """True if word looks like a strain/accession designation."""
    if not word:
        return False
    w = word.lower().rstrip('.').rstrip('-')
    if w in STRAIN_PREFIXES:
        return True
    # All-caps / alphanumeric IDs such as NOBI-1, DSM 5219, Mc-S-70
    if re.match(r'^[A-Z0-9][A-Z0-9\-/]*$', word):
        return True
    # Starts with a digit
    if re.match(r'^\d', word):
        return True
    return False

def is_valid_epithet(word):
    """A real species epithet: lowercase letters, possibly hyphenated."""
    return bool(word) and bool(re.match(r'^[a-z][a-z\-]+$', word))

def normalize_species(genus, species):
    """Return a long-format species string given genus and raw species."""
    if not species or species.lower() in {'', 'unknown', 'unidentified'}:
        return f'{genus} sp.' if genus else ''

    s = species.strip()

    # Pull out a leading "Candidatus" and propagate it to the genus.
    if s.lower().startswith('candidatus '):
        s = s[11:].strip()
        if genus and not genus.lower().startswith('candidatus '):
            genus = f'Candidatus {genus}'

    words = s.split()
    if not words:
        return f'{genus} sp.' if genus else ''

    first = words[0]

    # "sp.", "sp", "spp." -> genus sp.
    if first.lower().rstrip('.') in {'sp', 'spp'}:
        return f'{genus} sp.' if genus else 'sp.'

    # Determine the bare genus name for comparison.
    genus_base = ''
    if genus:
        genus_base = genus.lower().replace('candidatus ', '').strip()

    epithet = None

    # Species already includes the genus name: "Methanosarcina mazei ..."
    if genus_base and first.lower() == genus_base:
        if len(words) >= 2 and is_valid_epithet(words[1]) and not looks_like_strain(words[1]):
            epithet = words[1]
        else:
            return f'{genus} sp.' if genus else 'sp.'
    elif is_valid_epithet(first) and not looks_like_strain(first):
        epithet = first

    if not epithet:
        return f'{genus} sp.' if genus else 'sp.'

    return f'{genus} {epithet}' if genus else epithet

def fix_taxonomy_line(line):
    """Apply all fixes to one taxonomy TSV line."""
    line = line.rstrip('\n')
    if not line or line.startswith('#') or line.startswith('Feature ID'):
        return line

    cols = line.split('\t')
    if len(cols) < 2:
        return line

    taxon = cols[1]

    # Split taxonomy into ranks and drop superkingdom (sk__).
    parts = [p.strip() for p in taxon.split(';')]
    filtered = [p for p in parts if not p.startswith('sk__')]

    # Build rank lookup.
    ranks = {}
    for p in filtered:
        if '__' in p:
            prefix, value = p.split('__', 1)
            ranks[prefix.strip()] = value.strip()

    genus = ranks.get('g', '')
    species = ranks.get('s', '')
    ranks['s'] = normalize_species(genus, species)

    # Rebuild taxonomy string preserving original order (minus sk__ ranks).
    new_parts = []
    for p in filtered:
        if '__' in p:
            prefix, _ = p.split('__', 1)
            prefix = prefix.strip()
            if prefix in ranks:
                new_parts.append(f'{prefix}__{ranks[prefix]}')
            else:
                new_parts.append(p)
        else:
            new_parts.append(p)

    cols[1] = '; '.join(new_parts)
    return '\t'.join(cols)

with open(DEFAULT_INPUT_TSV, encoding='utf-8') as fin,      open(DEFAULT_OUTPUT_TSV, 'w', encoding='utf-8') as fout:
    for line in fin:
        fout.write(fix_taxonomy_line(line) + '\n')

print(f'Wrote fixed taxonomy to {DEFAULT_OUTPUT_TSV}')

