#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# ============================================================
# USER-CONFIGURABLE SETTINGS
# ============================================================

# The gene_target parameter is used in the output ID format:
#   {accession}_{gene_target}_{match_type}_{n}
# Example IDs:
#   AP019780.1_mcrA_gene_1
#   BAAJPY010000024.1_mrtA_product_or_note_1

THREADS = 8                                                    # Number of threads for parallel processing
input_dir = Path('02a_genomes/genbank_files')                  # Update this path
output_file = Path('03_genescoop/genescoop_mcrA.fasta')        # Update this path
gene_target = 'mcrA'                                           # Set to 'mcrA' or 'mrtA'
genes_of_interest = ['mcrA']                                   # Gene names (e.g. 'mcrA', 'mrtA')

# Keywords for matching in 'product' and 'note' fields
product_and_note_keywords = [
    'coenzyme-b sulfoethylthiotransferase subunit alpha',
    'methylcoenzyme M reductase subunit A', 'methyl coenzyme M reductase subunit A',
    'methyl coenzyme M reductase subunit alpha', 'methyl coenzyme M reductase alpha subunit',
    'methyl-coenzyme M reductase subunit A', 'methyl-coenzyme M reductase alpha subunit',
    'methyl coenzyme M reductase, subunit A', 'methyl coenzyme M reductase, subunit alpha',
    'methyl coenzyme M reductase, alpha subunit', 'methyl-coenzyme M reductase, subunit A',
    'methyl-coenzyme M reductase, alpha subunit', 'methyl-coenzyme M reductase, subunit alpha'
]

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# Precompute lowercase versions for comparison
genes_of_interest_lower = [g.lower() for g in genes_of_interest]
product_and_note_keywords_lower = [k.lower() for k in product_and_note_keywords]

def extract_sequences(genbank_file, gene_names, keyword_list, gene_target):
    """Extract sequences matching gene name or keywords from a GenBank file."""
    extracted = []
    processed_coords = set()

    try:
        gb_obj = SeqIO.read(genbank_file, 'genbank')
    except Exception as e:
        print(f"Error reading {genbank_file.name}: {e}")
        return []

    cds_features = [f for f in gb_obj.features if f.type == 'CDS']
    print(f" {genbank_file.name}: {len(cds_features)} CDS features")

    gene_hits = []

    for feature in cds_features:
        coords = (feature.location.start, feature.location.end)
        if coords in processed_coords:
            continue

        qualifiers = feature.qualifiers
        gene_qual = ' '.join(qualifiers.get('gene', [])).lower()
        product_note = ' '.join(
            qualifiers.get('product', []) + qualifiers.get('note', [])
        ).lower()

        if any(gene in gene_qual for gene in gene_names):
            gene_hits.append((feature, 'gene'))
            processed_coords.add(coords)
            continue

        if any(keyword in product_note for keyword in keyword_list):
            gene_hits.append((feature, 'product_or_note'))
            processed_coords.add(coords)

    print(f" Matches found: {len(gene_hits)}")

    for idx, (hit, match_type) in enumerate(gene_hits):
        try:
            seq = hit.extract(gb_obj)
            definition = gb_obj.description
            taxonomy = '; '.join(gb_obj.annotations.get('taxonomy', [])) or 'unknown taxonomy'
            desc = f"{match_type} {definition}, {taxonomy}"

            # New ID format: {accession}_{gene_target}_{match_type}_{n}
            # Example: AP019780.1_mcrA_gene_1 or BAAJPY010000024.1_mrtA_product_or_note_1
            record = SeqRecord(
                seq.seq,
                id=f"{gb_obj.id}_{gene_target}_{match_type}_{idx+1}",
                description=desc
            )
            extracted.append(record)
        except Exception as e:
            print(f" Error extracting feature: {e}")

    return extracted

# ============================================================
# DISCOVER GENBANK FILES
# ============================================================

gbk_files = sorted(input_dir.glob("*.gbk"))
print(f" Found {len(gbk_files)} GenBank files in {input_dir}")
gbk_files[:5]  # preview first few paths

# ============================================================
# MAIN EXTRACTION (PARALLEL)
# ============================================================

all_sequences = []

with ThreadPoolExecutor(max_workers=THREADS) as executor:
    future_to_file = {
        executor.submit(
            extract_sequences,
            file,
            genes_of_interest_lower,
            product_and_note_keywords_lower,
            gene_target  # Pass gene_target for ID construction
        ): file
        for file in gbk_files
    }

    for future in tqdm(
        as_completed(future_to_file),
        total=len(future_to_file),
        desc="Processing files (parallel)"
    ):
        try:
            sequences = future.result()
            all_sequences.extend(sequences)
        except Exception as e:
            file = future_to_file[future]
            print(f" Error processing {file.name}: {e}")

print(f"\n Total extracted sequences: {len(all_sequences)}")

# ============================================================
# VALIDATION / QUICK INSPECTION
# ============================================================

# Fix missing IDs / descriptions
for seq in all_sequences:
    if not seq.id or seq.id == "<unknown id>":
        seq.id = "unknown_id"
    if not seq.description or seq.description == "<unknown description>":
        seq.description = "unknown description"

# Quick look at first few records (optional)
# Note: IDs now include gene_target, e.g.: AP019780.1_mcrA_gene_1
for s in all_sequences[:5]:
    print(s.id, "|", s.description[:120], "...")

# Save to FASTA
try:
    SeqIO.write(all_sequences, output_file, "fasta")
    print(f" Sequences saved to {output_file}")
except Exception as e:
    print(f" Error writing FASTA: {e}")

