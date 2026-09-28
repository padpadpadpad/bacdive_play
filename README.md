# BacDive

Integrating bacterial phenotypic data from [BacDive](https://bacdive.dsmz.de/) (morphology, physiology, metabolism) with genomic and taxonomic data from NCBI and [GTDB](https://gtdb.ecogenomic.org/), to explore relationships between bacterial traits, genome features, and phylogeny.

## Workflow

1. **Query BacDive** for bacterial traits and strain records using species names or accessions (`scripts/bacdive_harvester-main/`, `scripts/api_use_new.R`, `scripts/using_the_api.R`).
2. **Link** BacDive records to NCBI genome assemblies via accession mappings.
3. **Download** genome assemblies/stats from NCBI (`scripts/ncbi_download.sh`).
4. **Retrieve genomic and taxonomic context** from GTDB (`data/gtdb/`).
5. **Wrangle and visualize** trait data across the bacterial phylogeny (`scripts/wrangle_phenotype_data.R`), producing figures in `figures/`.
6. **Analyze** links between phenotype and genome features (`scripts/check_sequencing_info.R`).

## Repository structure

```
BacDive/
├── scripts/
│   ├── bacdive_harvester-main/    # Toolkit for extracting BacDive traits (see its own README)
│   ├── api_use_new.R              # BacDive API usage / experimentation
│   ├── using_the_api.R            # Original BacDive API usage tutorial
│   ├── check_sequencing_info.R    # Genome accession analysis for BacDive records
│   ├── wrangle_phenotype_data.R   # Phenotype wrangling + phylogeny visualization
│   ├── ncbi_download.sh           # Downloads genome assemblies from NCBI
│   └── sparkql_query.txt.R        # Example SPARQL query against BacDive
├── data/
│   ├── bacdive_all_info.rds           # Cached full BacDive dataset (all records)
│   ├── bacdive_hasLength.csv          # BacDive records with cell length measurements
│   ├── bacdive_ncbi_accession.csv/txt # BacDive ID -> NCBI genome accession mapping
│   ├── bacdive_ncbi_genome_info.csv   # Genome metadata (length, GC%, accession)
│   ├── bacdive_patric_accession.csv   # BacDive ID -> PATRIC accession mapping
│   ├── bacdive_img_accession.csv      # BacDive ID -> IMG accession mapping
│   ├── ncbi_download/                 # Downloaded NCBI assembly stats files
│   └── gtdb/                          # GTDB taxonomy, metadata, and phylogenetic tree
│       ├── bac120_metadata.tsv        # ⚠️ 909 MB — too large for GitHub, not tracked
│       └── bac120_metadata.tsv.gz     # ⚠️ 225 MB — too large for GitHub, not tracked
├── figures/                       # Generated plots (body size, genome size, Topt vs. phylogeny)
├── data.tsv                       # NCBI assembly metadata table
└── BacDive.Rproj                  # RStudio project file
```

## Data access

Querying BacDive requires an API account: register at https://api.bacdive.dsmz.de and set credentials as expected by `scripts/bacdive_harvester-main/` and the API usage scripts.

Genome downloads use [`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download) via conda (see `scripts/ncbi_download.sh`).

## License

Scripts are marked with a modified MIT non-AI license (see individual file headers); this restricts use of the code for training machine learning / AI models.
