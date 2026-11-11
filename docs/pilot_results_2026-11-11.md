# Three-genome annotated small-protein pilot — 11 November 2026

**Source workflow:** [successful corrected run #43071848450](https://github.com/shixiaoxuan8540/bsubtilis-small-protein-pangenome/actions/runs/43071848450). Source artifact: `bsubtilis-annotated-small-protein-pilot` (includes `annotated_small_proteins.csv`, `assembly_qc.csv`, and `summary.json`).

## Genome and annotation inventory

| Assembly | Strain / source label | Whole-genome FASTA bp | FASTA sequences | Annotated proteins 1–100 aa | Distinct short-protein AA sequences |
|---|---|---:|---:|---:|---:|
| GCF_000009045.1 | 168, reference | 4,215,606 | 1 | 636 | 635 |
| GCA_000321395.1 | BSP1, chicken feces | 4,043,754 | 1 | 430 | 430 |
| GCA_000209795.2 | BEST195, natto | 4,111,218 | 2 | 954 | 951 |

The initial run #43071193793 mistakenly selected `cds_from_genomic.fna` instead of the whole-genome FASTA when reporting assembly size and contig number for two assemblies. This was fixed in commit `c6c5bc9`, and the corrected run's QC numbers above supersede the original QC output. The short-protein inventories and exact-sequence counts were unchanged.

## Exact amino-acid sequence overlaps (exploratory only)

Each distinct sequence is counted once per strain. Presence-bit order: **168, BSP1, BEST195**.

| Present in | Count |
|---|---:|
| 168 only | 340 |
| BSP1 only | 199 |
| BEST195 only | 748 |
| 168 + BSP1, not BEST195 | 107 |
| 168 + BEST195, not BSP1 | 79 |
| BSP1 + BEST195, not 168 | 15 |
| All three | 109 |
| **Union of distinct exact protein sequences** | **1,597** |

These results describe **perfect amino-acid sequence identity**, **not orthologous gene families**. Even a single amino-acid change makes two genuine orthologs count as different. Consequently, the 109 identical sequences are NOT a core small-protein pangenome estimate, and the 748 BEST195-only sequences are NOT proof of natto-associated gene acquisition.

## Major biological caveat

BEST195 has 954 annotated proteins <=100 aa, versus 636 in 168 and 430 in BSP1. A substantial portion are especially short (BEST195: 388 below 50 aa; 168: 118; BSP1: 47). Differences in automated gene calling and short-ORF annotation may dominate the apparent strain-specific content. Repeat all annotation using one method, then perform sensitive nucleotide/protein homology searches and curate ortholog groups before making ecological interpretations.

## Next analysis

1. Confirm that assemblies represent distinct, correctly classified B. subtilis strains and review original isolation metadata.
2. Standardize genome annotations, explicitly retaining short ORFs.
3. Search all genome nucleotide sequences for known B. subtilis small-protein homologs, including unannotated loci.
4. Cluster and manually verify ambiguous ortholog families, gene lengths, and fragmented loci.
5. Construct a gene-family presence/absence matrix and calculate core/accessory/niche prevalence with a phylogeny.
6. Add more independent isolates before any niche association inference.

**No phenotype, survival advantage, or ecological causality is established by this pilot.**
