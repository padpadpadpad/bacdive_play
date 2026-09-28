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
librarian::shelf(tidyverse)
## ---------------------------

# load data
d <- readRDS('data/bacdive_all_info.rds')
length(d$bacdiveid)

# extract sequence information
d_seq_info <- select(d, bacdiveid, sequence_info) %>%
  unnest(sequence_info) %>%
  filter(!is.na(description))
length(unique(d_seq_info$bacdiveid))

# count the number of times ncbi is there
d_ncbi <- d_seq_info %>%
  filter(database == 'ncbi')

unique(d_seq_info$bacdiveid) %>%
  length()
unique(d_ncbi$bacdiveid) %>%
  length()

# find the sequences that are not in ncbi
d_not_ncbi <- d_seq_info %>%
  filter(!bacdiveid %in% d_ncbi$bacdiveid)
# these need to be updated by the looks of it!

# save out the sequences for NCBI, JGI, and patric separately
write.csv(
  select(d_ncbi, bacdiveid, accession),
  'data/bacdive_ncbi_accession.csv',
  row.names = FALSE
)
pull(d_ncbi, accession) %>%
  paste(., '.1', sep = '') %>%
  write.table(
    'data/bacdive_ncbi_accession.txt',
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

d_not_ncbi %>%
  filter(database == 'patric') %>%
  select(bacdiveid, accession) %>%
  write.csv('data/bacdive_patric_accession.csv', row.names = FALSE)
d_not_ncbi %>%
  filter(database == 'img') %>%
  select(bacdiveid, accession) %>%
  write.csv('data/bacdive_img_accession.csv', row.names = FALSE)


# download and format the NCBI data

# list files
files <- list.files(
  "data/ncbi_download/genbank/bacteria",
  full.names = TRUE,
  recursive = TRUE,
  pattern = '.txt'
)
length(files)

file <- files[84]

get_assembly_data <- function(file) {
  temp <- read.table(file, comment.char = '#', nrows = 13, fill = TRUE)

  # select key columns
  temp <- temp %>%
    select(V5, V6) %>%
    filter(V5 %in% c('total-length', 'gc-perc')) %>%
    pivot_wider(
      names_from = V5,
      values_from = V6,
      values_fn = list(V6 = as.character)
    ) %>%
    janitor::clean_names() %>%
    mutate(accession = basename(dirname(file)))

  return(temp)
}

genome_output <- files %>%
  map(., get_assembly_data) %>%
  list_rbind()

# save out genome data
write.csv(
  genome_output,
  'data/bacdive_ncbi_genome_info.csv',
  row.names = FALSE
)
