#----------------------------------
# Purpose of the script: Play with BacDive scripts
#
# What this script does:
# 1.
# 2.
# 3.
#
# Author: Dr. Daniel Padfield
# Date Created: 19 September 2025
#
# Daniel Padfield, 2025
# This code is licensed under a modified MIT non-AI license. The code and any modifications made to it may not be used for the purpose of training or improving machine learning algorithms, including but not limited to artificial intelligence, natural language processing, or data mining. This condition applies to any derivatives, modifications, or updates based on the Software code. Any usage of the Software in an AI-training dataset is considered a breach of this License.
# The full license can be found here: https://github.com/padpadpadpad/non-ai-licenses/blob/main/NON-AI-MIT
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
