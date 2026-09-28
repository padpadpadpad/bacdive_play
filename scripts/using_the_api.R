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
if (!requireNamespace("librarian", quietly = TRUE)) {
  install.packages("librarian")
}

# if BiocManager is not installed, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# if Biobase is not installed, install it from Bioconductor
if (!requireNamespace("Biobase", quietly = TRUE)) {
  BiocManager::install("Biobase")
}

# load packages
librarian::shelf(tidyverse, BacDive)

## ---------------------------

# read in records that have length as defined by sparkql
d_spark <- read.csv('data/bacdive_hasLength.csv')
head(d_spark)

# open API
source('scripts/bacdive_password.R')

bacdive <- open_bacdive(
  username = 'd.padfieldscfc@gmail.com',
  password = bacdive_password
)

test <- d_spark$bacdiveid[1]

id <- fetch(bacdive, test) %>%
  as.data.frame()

list_columns <- id %>% purrr::keep(is.list) %>% names()

# functions to fetch data from individual sections of the BacDive database

# general
get_general <- function(x) {
  temp <- x$General %>%
    unlist() %>%
    tibble(names = names(.), values = .) %>%
    filter(names %in% c('BacDive-ID')) %>%
    pivot_wider(id_cols = NULL, names_from = names, values_from = values) %>%
    janitor::clean_names()

  return(temp)
}

get_general(x)

# taxonomy
get_taxonomy <- function(x) {
  # if NULL then return NA
  if (is.null(x$`Name and taxonomic classification`)) {
    return(NA)
  }
  temp <- x$`Name and taxonomic classification` %>%
    unlist() %>%
    tibble(names = names(.), values = .) %>%
    filter(
      names %in%
        c(
          'domain',
          'phylum',
          'class',
          'order',
          'family',
          'genus',
          'species',
          'full scientific name'
        )
    ) %>%
    pivot_wider(id_cols = NULL, names_from = names, values_from = values) %>%
    janitor::clean_names()

  return(temp)
}

get_taxonomy(x)

# get morphology
get_morphology <- function(x) {
  # if NA then return NA
  if (any(is.na(x$Morphology[[1]]))) {
    return(NA)
  }

  # if NULL then return NA
  if (is.null(x$Morphology[[1]])) {
    return(NA)
  }
  temp <- x$Morphology[[1]] %>%
    unlist() %>%
    tibble(names = names(.), values = .) %>%
    filter(
      str_detect(
        names,
        paste(
          c(
            'cell length',
            'cell width',
            'motility',
            'gram stain',
            'cell shape'
          ),
          collapse = '|'
        )
      )
    ) %>%
    # only keep things after the last .
    mutate(names = sub(".*\\.", "", names)) %>%
    distinct() %>%
    pivot_wider(id_cols = NULL, names_from = names, values_from = values) %>%
    janitor::clean_names()

  return(temp)
}

get_morphology(x)

# get temperature conditions
get_temps <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Culture and growth conditions`[[1]]))) {
    return(NA)
  }

  # if NULL then return NA
  if (is.null(x$`Culture and growth conditions`[[1]]$`culture temp`)) {
    return(NA)
  }
  temp <- x$`Culture and growth conditions`[[1]]$`culture temp` %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    select(-ref)

  return(temp)
}

get_temps(x)

# get pH conditions
get_ph <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Culture and growth conditions`[[1]]))) {
    return(NA)
  }

  # if NULL then get pH conditions
  if (is.null(x$`Culture and growth conditions`[[1]]$`culture pH`)) {
    return(NA)
  }
  temp <- x$`Culture and growth conditions`[[1]]$`culture pH` %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    select(-ref)

  return(temp)
}

get_ph(x)

# get halophily information
get_halophily <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Physiology and metabolism`[[1]]))) {
    return(NA)
  }

  # if NULL then return NA
  if (is.null(x$`Physiology and metabolism`[[1]]$halophily)) {
    return(NA)
  }
  temp <- x$`Physiology and metabolism`[[1]]$halophily %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    select(-ref)

  return(temp)
}

get_halophily(x)

# create safe version of this function
safe_get_halophily <- safely(get_halophily, otherwise = NA)

# get isolation information
get_isolation <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Isolation, sampling and environmental information`))) {
    return(NA)
  }

  # if  NULL then return NA
  if (
    is.null(
      x$`Isolation, sampling and environmental information`[[
        1
      ]]$isolation
    )
  ) {
    return(NA)
  }
  temp <- x$`Isolation, sampling and environmental information`[[
    1
  ]] %>%
    unlist() %>%
    tibble(names = names(.), values = .) %>%
    group_by(names) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    filter(
      str_detect(
        names,
        paste(
          c(
            'continent',
            'country',
            'geographic location',
            'sample type',
            'isolation date'
          ),
          collapse = '|'
        )
      )
    ) %>%
    filter(names != 'origin.country') %>%
    # only keep things after the last .
    mutate(names = sub(".*\\.", "", names)) %>%
    pivot_wider(id_cols = NULL, names_from = names, values_from = values) %>%
    janitor::clean_names()

  return(temp)
}

get_isolation(x)

# get 16s abundance info
get_16s_abundance <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Isolation, sampling and environmental information`))) {
    return(NA)
  }

  # if NULL then return NA
  if (
    is.null(
      x$`Isolation, sampling and environmental information`[[
        1
      ]]$`taxonmaps`
    )
  ) {
    return(NA)
  }
  temp <- x$`Isolation, sampling and environmental information`[[
    1
  ]]$`taxonmaps` %>%
    unlist() %>%
    tibble(sample = names(.), count = .) %>%
    filter(
      sample %in%
        c('animal counts', 'aquatic counts', 'plant counts', 'soil counts')
    ) %>%
    mutate(sample = sub(" counts", "", sample)) %>%
    janitor::clean_names()

  return(temp)
}

get_16s_abundance(x)

# get aerobic information
get_aerobe <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Physiology and metabolism`[[1]]))) {
    return(NA)
  }

  # if NULL then return NA
  if (is.null(x$`Physiology and metabolism`[[1]]$`oxygen tolerance`)) {
    return(NA)
  }
  temp <- x$`Physiology and metabolism`[[1]]$`oxygen tolerance` %>%
    unlist() %>%
    unlist() %>%
    tibble(names = names(.), values = .)

  # if names does not exist, create it
  if (!('names' %in% names(temp))) {
    temp <- mutate(temp, names = 'oxygen tolerance')
  }

  temp <- group_by(temp, names) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    filter(
      names %in%
        c(
          'oxygen tolerance'
        )
    ) %>%
    pivot_wider(id_cols = NULL, names_from = names, values_from = values) %>%
    janitor::clean_names()

  return(temp)
}

get_aerobe(x)

# create safe function
safe_get_aerobe <- safely(get_aerobe, otherwise = NA)

get_sequence_info <- function(x) {
  # if NA then return NA
  if (any(is.na(x$`Sequence information`))) {
    return(NA)
  }

  # if NULL then return NA
  if (is.null(x$`Sequence information`[[1]]$`Genome sequences`)) {
    return(NA)
  }
  temp <- x$`Sequence information`[[1]]$`Genome sequences` %>%
    bind_rows() %>%
    janitor::clean_names() %>%
    select(-ref)

  return(temp)
}

get_sequence_info(x)

# create empty dataset for all info
all_ids <- data.frame(
  bacdiveid = d_spark$bacdiveid,
  stringsAsFactors = FALSE
) %>%
  mutate(
    taxonomy = list(NA),
    morphology = list(NA),
    temps = list(NA),
    ph = list(NA),
    halophily = list(NA),
    isolation = list(NA),
    abund_16s = list(NA),
    aerobe = list(NA),
    sequence_info = list(NA)
  )

# open API
source('scripts/bacdive_password.R')

bacdive <- open_bacdive(
  username = 'd.padfieldscfc@gmail.com',
  password = bacdive_password
)

for (i in 3954:nrow(all_ids)) {
  print(all_ids$bacdiveid[i])

  # fetch data for each ID
  x <- fetch(bacdive, all_ids$bacdiveid[i]) %>%
    as.data.frame()

  # get taxonomy
  all_ids$taxonomy[[i]] <- get_taxonomy(x)

  # get morphology
  all_ids$morphology[[i]] <- get_morphology(x)

  # get temperature conditions
  all_ids$temps[[i]] <- get_temps(x)

  # get pH conditions
  all_ids$ph[[i]] <- get_ph(x)

  # get halophily information
  all_ids$halophily[[i]] <- get_halophily(x)

  # get isolation information
  all_ids$isolation[[i]] <- get_isolation(x)

  # get 16s abundance info
  all_ids$abund_16s[[i]] <- get_16s_abundance(x)

  # get aerobic information
  all_ids$aerobe[[i]] <- get_aerobe(x)

  # get sequence info
  all_ids$sequence_info[[i]] <- get_sequence_info(x)
}

# save this out
saveRDS(all_ids, 'data/bacdive_all_info.rds')
