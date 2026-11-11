# Curated 589 ecological comparison — 2026-11-11

GitHub Actions: https://github.com/shixiaoxuan8540/bsubtilis-small-protein-pangenome/actions/runs/43109215171

## Methods
The supplied 589-entry manually curated B. subtilis small-protein reference contains 586 unique full-length amino-acid sequences. We scored a curated entry as present only if its full reference amino-acid sequence exactly matched a protein in the downloaded native annotation. The analysis used 11 genomes, including five added soil/plant-associated chromosomes; plasmids were excluded from those five and NCIB 3610, but assembly-level annotations in some older isolates may still include plasmid proteins. These are **annotation-level exact matches, not gene presence/absence, orthology, or sequence-similarity conservation**.

## Counts
| Genome | Curated matches / 589 |
|---|---:|
| Reference 168 | 575 |
| NCIB 3610 | 554 |
| Chicken BSP1 | 216 |
| Natto BEST195 | 182 |
| Aquatic sediment MARUCo01 | 330 |
| Natto NARUSE | 252 |
| Soil 73 | 326 |
| Soil 75 | 332 |
| Soil MB8_B1 | 329 |
| Plant root UD1022 | 259 |
| Plant endophyte BSn5 | 298 |

66 curated entries have exact annotated-protein matches in all 11 genomes; 11 have none. The original 589 curated entries include three duplicated amino-acid sequences, so entry counts and unique peptide counts are not identical. Eramic25 GCA_031932165.1 lacked a downloadable protein.faa and was excluded.

## Figures and files
Download the workflow artifact `curated-589-ecological-sequence-comparison` from the run above. It contains:
- `curated_589_proportional_venn.png` — 168, NCIB 3610, soil 73; area-proportional approximation
- `curated_589_upset.png` — top 30 nonempty exact match patterns across 11 genomes
- `curated_589_presence_matrix.csv` — 589 rows x 11 genomes
- `curated_589_genome_summary.csv`, `curated_589_intersections.csv`, `assembly_qc.csv`

## Limitations / next analysis
A failure to match exactly may reflect amino-acid substitutions, missing CDS annotation, gene loss, or assembly errors. Next perform a genome-level translated search and standardized reannotation, reconcile the three duplicate reference sequences, normalize chromosome-only comparisons across every isolate, and account for phylogenetic relatedness before claiming niche association.
