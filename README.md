# Bacillus subtilis small-protein pangenome (pilot)

Goal: compare small-protein gene families across 20 carefully curated B. subtilis genomes and explore associations with isolation niches. This is a **pilot**, not a complete species-wide pangenome.

## Sampling plan
Two reference strains (168 and NCIB 3610), plus three aquatic, three soil, three human gut, three chicken gut, three fish gut, and three food isolates. References are not ecological replicates.

## Workflow
1. Curate and independently verify NCBI assembly and BioSample metadata, including original isolation source.
2. Confirm taxonomic identity and assembly quality; remove redundant isolates.
3. Download genome FASTA, GFF3, CDS, and protein sequences using NCBI Datasets.
4. Extract annotated proteins <=100 aa as a **preliminary baseline only**.
5. Recover known short-ORF homologs from all genomes, including unannotated ORFs.
6. Build curated ortholog-family presence/absence matrix.
7. Calculate strict core, soft core, accessory, singleton, and niche-specific sets.
8. Visualize Venn/UpSet, phylogeny-aligned heatmap, and family identities.

**Important:** No environmental isolate set has been validated, no genomes downloaded, and no pangenome results calculated. Do not interpret annotation absence as biological gene absence.

The genome-selection CSV and remaining scripts are included in the downloadable starter ZIP shared in chat.
