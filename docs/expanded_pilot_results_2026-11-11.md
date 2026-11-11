# Expanded annotated short-protein pilot (2026-11-11)

Source: [GitHub Actions run 43073213469](https://github.com/shixiaoxuan8540/bsubtilis-small-protein-pangenome/actions/runs/43073213469), artifact `bsubtilis-expanded-annotated-small-protein-pilot`.

## Results

| Strain | Accession | Annotated proteins <=100 aa | Genome FASTA bp | FASTA sequences |
|---|---|---:|---:|---:|
| 168 | GCF_000009045.1 | 636 | 4,215,606 | 1 |
| NCIB 3610 chromosome only | CP020102 | 643 | 4,215,607 | 1 |
| BSP1 | GCA_000321395.1 | 430 | 4,043,754 | 1 |
| BEST195 | GCA_000209795.2 | 954 | 4,111,218 | 2 |
| MARUCo01 | GCA_025837035.1 | 583 | 4,061,384 | 2 |
| NARUSE | GCA_030297075.1 | 743 | 4,211,563 | 2 |
| Eramic25 | GCA_031932165.1 | not available | not available | not available |

The Eramic25 assembly was included in the download request but lacks `protein.faa` in the returned dataset. It was excluded from sequence comparisons; investigate assembly annotation availability before interpreting absence.

The NCIB 3610 sequence was downloaded from chromosome accession **CP020102**. The plasmid **CP020103** was explicitly excluded. Other strains were not chromosome-only filtered.

## Exact sequence analysis

- Total annotated short-protein records among six processed genomes: **3,989**.
- Distinct exact amino-acid sequences (union): **2,134**.
- Exact amino-acid sequences shared by all six: **90**.
- Exact amino-acid sequences shared between 168 and NCIB 3610: **556**.

These counts do **not** represent orthologous gene families or a biologically validated pangenome. A single substitution prevents an exact match; differences in short-ORF annotation can dominate counts. Genome FASTA sequence counts are not proof of assembly completeness. The next phase must harmonize annotations, recover unannotated homologs, and verify taxonomy and source metadata.
