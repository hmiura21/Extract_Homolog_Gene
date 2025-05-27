# 🔍 Extract_Homolog_Gene

This repository contains an enhanced version of the `find_homologs.py` script. It identifies putative homologous genes from BLAST results, maps them to gene annotations in a BED file, extracts their sequences from a reference assembly, and writes the results to a FASTA file. If a gene is encoded on the minus strand, the script will return the reverse complement of the sequence.

---

## 📌 Assignment Description

This project builds on a previous homolog identification task and emphasizes working with and modifying existing code. The final script accepts four command-line arguments:

1. **BLAST output file**: Tabular BLAST results
2. **BED file**: Gene annotations
3. **FASTA file**: Assembly from which to extract homologous gene sequences
4. **Output file**: FASTA-format output of homologous gene sequences

---

## ✅ Features

- Filters BLAST hits by:
  - ≥30% identity
  - ≥90% of query length
- Maps filtered hits to genes based on BED file coordinates
- Extracts matched gene sequences from the assembly
- Handles reverse complementing for minus-strand genes
- Outputs homologous sequences in FASTA format

---

## 🔧 Usage

```bash
python extract_homolog_gene.py <querydata.faa> <annotations.bed> <assembly.fna> <output_genes.fasta>
```
---

## 📁 Repository Structure
- `Vibrio_cholerae_N16961.fna`: example bacteria assembly for Vibrio cholerae (example query data)
- `Vibrio_cholerae_N16961.bed`: example bacteria bed file for Vibrio cholerae containing gene names (example bed file)
- `HK_domain.faa`: sequence of histidine kinase (HK) domains from the organism Escherichia coli strain K-12 (example subject file)
- `extract_homolog_gene.sh`: script that contains all commands


