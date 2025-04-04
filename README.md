## 🧬 Curated _mcrA_ database

A curated collection of _mcrA_ and _mrtA_ sequences assembled for taxonomical classification of _mcrA_ reads from amplicon sequencing.

---

## 📂 Contents

- mcrA_SRD.fas – suitable for 2×150 bp sequencing (Short Read Database, e.g. MiniSeq, qmcrAF and mcrA-rev primers)
- mcrA_LRD.fas – suitable for 2×300 or 2×250 bp sequencing (Long Read Database, e.g. MiSeq, Mlas-mod-F and mcrA-rev primers)
- LICENSE – Licensing information  
- README.md – This file  

---

## 🧫 Description

This database contains curated DNA _mcrA_ and _mrtA_ sequences sourced from the NCBI database using [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) and sequences extracted from genomes sourced from NCBI and extracted using [GeneScoop](https://github.com/Zblezble/GeneScoop) and followed up by a manual review.
The goal is to provide a high-quality, reusable dataset for the community. Only sequences featuring species annotations are present.

Featururing:
- LRD features 478 unique sequences with median length of 425 bp encompassing 188 different methanogen species 
- SRD features 411 unique sequences with median length of 211 bp encompassing 182 different methanogen species 

---

## 🛠️ Usage

The database is in DADA2 comaptible format.

This database is released under a CC BY 4.0 license. See LICENSE for more details.
