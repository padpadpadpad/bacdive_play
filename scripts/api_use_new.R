#----------------------------------
# Purpose of the script:
#
# What this script does:
# 1.
# 2.
# 3.
#
# Author: Dr. Daniel Padfield
# Date Created: 28 September 2026
#
# Daniel Padfield, 2026
# This code is licensed under an MIT license and has been developed with the help of GitHub copilot, and potentially other LLMs. However, all output has been checked and modified by a human and I take full responsibility for any errors.
#
#----------------------------------
#
# Notes (potentially about software used:
#
#
#
#----------------------------------
# if librarian is not installed, install it
if (!requireNamespace('librarian', quietly = TRUE)) {
  install.packages('librarian')
}
# if BiocManager is not installed, install it
if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}
# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace('Biobase', quietly = TRUE)) {
  BiocManager::install('Biobase')
}
# load packages
librarian::shelf(tidyverse, BacDive)
## ---------------------------

# load in helper functions
source('scripts/bacdive_harvester-main/helper_functions.R')

# open API
bacdive <- open_bacdive(
  username = 'd.padfieldscfc@gmail.com',
  password = 'REDACTED'
)

# read in records that have length as defined by sparkql
d_spark <- read.csv('data/bacdive_hasLength.csv')
head(d_spark)

test <- d_spark$bacdiveid[2]

id <- fetch(bacdive, test) %>%
  data.frame()

id2 <- id[[1]]

# extract oxygen tolerance
extract_oxygen_tolerance(id)

# extract positive enzymes
extract_positive_enzymes(id)

# extract positive metabolites
extract_positive_metabolites(id)

# get preferred record
get_preferred_record(id)

# get type strain index
get_type_strain_index(id)

extract_morphology_trait(id, "cell shape")
extract_morphology_trait(id, "cell shape")

sp <- 'Escherichia coli'
id <- retrieve(bacdive, query = sp, sleep = 0.1)

# Select most authoritative record (preferring type strains)
idx <- get_type_strain_index(id)
preferred <- get_preferred_record(id)

extract_morphology_trait(preferred, "cell size")
