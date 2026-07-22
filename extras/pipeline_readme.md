# mcrA/mrtA reference database pipeline

This repository contains the reproducible pipeline used to build the curated short-read (SRD) and long-read (LRD) *mcrA* / *mrtA* reference databases.

All scripts are numbered in the order they are intended to run.

## Repository contents

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_get_ncbi_mcrA_mrtA.sh` | Download *mcrA* and *mrtA* nucleotide sequences from NCBI with RESCRIPt |
| 01b | `01b_fix_rescript_taxonomy.py` | Clean RESCRIPt taxonomy for later merging |
| 02a | `02a_download_genomes.sh` | Download completed genomes / WGS assemblies containing *mcrA* |
| 02b | `02b_download_genomes.sh` | Download completed genomes / WGS assemblies containing *mrtA* |
| 03a | `03a_GeneScoop.py` | Extract *mcrA* sequences from GenBank files downloaded in 02a |
| 03b | `03b_GeneScoop.py` | Extract *mrtA* sequences from GenBank files downloaded in 02b |
| 03c | `03c_fix_genescoop_taxonomy.py` | Build clean FASTA + taxonomy TSV from 03a *mcrA* output |
| 03d | `03d_fix_genescoop_taxonomy.py` | Build clean FASTA + taxonomy TSV from 03b *mrtA* output |
| 04 | `04_merge_and_filter.sh` | Merge RESCRIPt + GeneScoop data, filter by length and quality |
| 05 | `05_dereplicate.sh` | Dereplicate merged sequences with RESCRIPt |
| 06 | `06_mafft_align.sh` | Align dereplicated sequences with MAFFT and add guide sequences |
| 07 | `07_trim_guide.sh` | Trim alignment to the *Methanosarcina barkeri* guide span |
| 08 | `08_filter_flag.sh` | Flag low-quality, short, and phylogenetically outlier sequences |
| 09 | `09_remove_flagged.sh` | Remove flagged sequences from the alignment |
| 10 | `10_manual_remove.R` | Remove manually curated sequences listed in `manual_removal.txt` |
| 11 | `11_dereplicate_final.sh` | Dereplicate the manually cleaned database |
| 12 | `12_amplicon_derep.sh` | Trim to short- and long-read amplicon regions and dereplicate |
| 13 | `13_taxonomy_test.R` | Validate SRD/LRD as DADA2 taxonomy references |

## User-provided files

The following files are expected to exist in the working directory (i.e. the directory you set as `STORAGE`):

- `guide_seqs.fas` — guide sequences including primers, representative short/long amplicons, and the full-length *Methanosarcina barkeri* reference (`guide_Y00158.1_Methanosarcina_barkeri`). This file is provided in the repository.
- `remove_quality.tsv`, `remove_length.tsv`, `remove_taxonomy.tsv` — lists of sequence IDs to remove, produced during manual review of step 08 outputs. Create empty files (with only a header) if nothing needs to be removed.
- `manual_removal.txt` — list of sequence IDs to remove in step 10, created during manual tree inspection.
- `new_taxonomy.tsv` — taxonomy file matching `alignment_final.fasta` from step 09. This is the curated taxonomy that accompanies the alignment going into step 10. It must have the sequence ID in the first column and the taxonomy string in the second column.

## Important manual dependency: step 11 → step 12

Step 11 (`11_dereplicate_final.sh`) produces an aligned, dereplicated database **without** guide sequences:

```
11_dereplicated_final/alignment_final_derep.fasta
```

Step 12 (`12_amplicon_derep.sh`) needs an alignment **with** guide sequences to locate the short- and long-read amplicon boundaries:

```
11_dereplicated_final/alignment_final_derep_guides.fasta
```

Before running step 12, you must re-add `guide_seqs.fas` to `alignment_final_derep.fasta`. 


This manual step is required because guide sequences are removed during dereplication but are needed again to define the amplicon trimming boundaries.

