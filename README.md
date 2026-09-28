# BacDive Play

BacDive is a database containing standardised phenotypic information of cultured bacterial isolates. This allows us to test interesting questions about bacterial ecology and evolution by integrating bacterial phenotypic data from [BacDive](https://bacdive.dsmz.de/) (morphology, physiology, metabolism) with genomic and taxonomic data from NCBI and [GTDB](https://gtdb.ecogenomic.org/).

For example, does genome size scale with cell size and temperature.

## Using this folder

This folder is currently highly messy and imperfect. The main script is **wrangle_phenotype_data.R** which does some introductory analyses of the compiled dataset.

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

## License

This code is licensed under an MIT license and has been developed with the help of GitHub copilot, and potentially other LLMs.