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
librarian::shelf(tidyverse, ggtree, ape, RColorBrewer, ggnewscale)
## ---------------------------

# load in data
d <- readRDS(file = "data/bacdive_all_info.rds")

# first wrangle the size data
d_size <- d %>%
  select(bacdiveid, morphology) %>%
  mutate(
    first_element_class = map_chr(morphology, function(x) class(x[[1]][1]))
  )

# when first element is a character, unnest
d_size2 <- d_size %>%
  filter(first_element_class == "character") %>%
  unnest(morphology)

# when first element is a list, try to unnest
d_size3 <- d_size %>%
  filter(first_element_class == "list")


# function to concatenate any elements that have multiple entries
concat_multiple <- function(x) {
  return(
    unlist(x) %>%
      paste(collapse = ", ")
  )
}

# for every column in the list column morphology, apply concat_multiple
d_size3 <- d_size3 %>%
  unnest(morphology)

d_size4 <- d_size3 %>%
  distinct() %>%
  mutate(,
    gram_stain = map_vec(gram_stain, concat_multiple),
    cell_length = map_vec(cell_length, concat_multiple),
    cell_width = map_vec(cell_width, concat_multiple),
    cell_shape = map_vec(cell_shape, concat_multiple),
    motility = map_vec(motility, concat_multiple)
  )

# check when there are multiple width or cell lengths
# keep when there is a , in either morphology_cell_width or morphology_cell_length
d_size5 <- d_size4 %>%
  filter(
    str_detect(cell_length, ",") |
      str_detect(cell_width, ",")
  )

# have a look at the easy ones
unique(d_size2$cell_shape) %>% clipr::write_clip()

d_size2 %>%
  group_by(cell_shape) %>%
  tally() %>%
  arrange(desc(n))

different_shapes <- c(
  "rod-shaped",
  "coccus-shaped",
  "NA",
  "star-shaped",
  "oval-shaped",
  "vibrio-shaped",
  "filament-shaped",
  "ovoid-shaped",
  "spiral-shaped",
  "curved-shaped",
  "pleomorphic-shaped",
  "helical-shaped",
  "ring-shaped",
  "other",
  "ellipsoidal",
  "sphere-shaped",
  "spore-shaped",
  "diplococcus-shaped"
)

# if there is a - calculate the mean of the two numbers
d_size2 <- d_size2 %>%
  mutate(
    cell_length = ifelse(
      str_detect(cell_length, "-"),
      map_chr(str_split(cell_length, "-"), ~ mean(parse_number(.x))),
      parse_number(cell_length)
    ),
    cell_width = ifelse(
      str_detect(cell_width, "-"),
      map_chr(str_split(cell_width, "-"), ~ mean(parse_number(.x))),
      parse_number(cell_length)
    )
  ) %>%
  mutate(across(
    c(cell_length, cell_width),
    as.numeric
  ))

# calculate biovolume
d_rod <- d_size2 %>%
  filter(cell_shape %in% c("rod-shaped", "coccus-shaped")) %>%
  mutate(
    biovolume = (pi / 4) * cell_width^2 * (cell_length - (cell_width / 3))
  )

# ok now we can look at the size data
filter(d_rod) %>%
  ggplot(aes(biovolume)) +
  geom_histogram(fill = 'white', col = 'black') +
  theme_bw(base_size = 16) +
  scale_x_log10() +
  labs(x = "Biovolume (µm³)")

# extract temperature data
d_temp <- d %>%
  select(bacdiveid, temps) %>%
  unnest(temps)

# get optimum temperature
d_temp_opt <- d_temp %>%
  filter(type == "optimum") %>%
  mutate(
    temperature = ifelse(
      str_detect(temperature, "-"),
      map_chr(str_split(temperature, "-"), ~ mean(parse_number(.x))),
      parse_number(temperature)
    )
  ) %>%
  mutate(temperature = as.numeric(temperature)) %>%
  select(bacdiveid, temperature) %>%
  distinct() %>%
  group_by(bacdiveid) %>%
  summarise(
    topt = mean(temperature, na.rm = TRUE),
    range = max(temperature, na.rm = TRUE) - min(temperature, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(range <= 3)

ggplot(d_temp_opt, aes(topt)) +
  geom_histogram(fill = 'white', col = 'black') +
  geom_vline(aes(xintercept = 37), linetype = "dashed", col = "red") +
  theme_bw(base_size = 16) +
  labs(x = "Optimum temperature (°C)")

# join these two datasets together
d_joined <- d_rod %>%
  left_join(., d_temp_opt, by = "bacdiveid")
d_joined <- filter(d_joined, range <= 3)


d_joined <- filter(d_joined, biovolume > 0) %>%
  filter(!is.na(topt))
d_joined <- mutate(d_joined, log_biovolume = log10(biovolume))

lm(log10(biovolume) ~ topt, data = d_joined, na.action = "na.omit") %>%
  summary()

# need to do phylogenetic mixed model on this, need to see if there is a better proxy for isolated temperature, but then there is body size was not measured for each species at each temperature, so maybe not a good idea to do this.

ggplot(d_joined, aes(topt, biovolume)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 16) +
  labs(
    x = "Optimum temperature (°C)",
    y = "Biovolume (µm³)"
  ) +
  scale_y_log10() +
  stat_smooth(method = 'lm', se = FALSE)

# read in genome data
d_genome <- read.csv(file = "data/bacdive_ncbi_genome_info.csv")

# read in d_ncbi genome info
d_ncbi_genome <- read.csv(file = "data/bacdive_ncbi_accession.csv") %>%
  mutate(
    accession = paste(accession, ".1", sep = "")
  )

d_ncbi_genome <- left_join(d_ncbi_genome, d_genome, by = "accession") %>%
  group_by(bacdiveid) %>%
  summarise(
    mean_size = mean(total_length, na.rm = TRUE),
    se_size = sd(total_length, na.rm = TRUE) / sqrt(n()),
    mean_gc = mean(gc_perc, na.rm = TRUE),
    se_gc = sd(gc_perc, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

d_joined <- left_join(d_joined, d_ncbi_genome, by = "bacdiveid")

d_joined <- filter(d_joined, !is.na(mean_size))

ggplot(d_joined, aes(topt, mean_size)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 16) +
  labs(
    x = "Optimum temperature (°C)",
    y = "Genome size (bp)"
  ) +
  stat_smooth(method = 'lm', se = FALSE) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(10^5.5, 10^7.5)
  )

ggplot(d_joined, aes(topt, mean_gc)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 16) +
  labs(
    x = "Optimum temperature (°C)",
    y = "Genome Composition (%)"
  ) +
  stat_smooth(method = 'lm', se = FALSE)

ggplot(d_joined, aes(biovolume, mean_size)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 16) +
  labs(
    x = "Biovolume (µm³)",
    y = "Genome size (bp)"
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(10^5.5, 10^7.5)
  ) +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(10^-2.75, 10^4.5)
  ) +
  stat_smooth(method = 'lm', se = FALSE)

lm(log10(mean_size) ~ topt, data = d_joined, na.action = "na.omit") %>%
  summary()

1 - 10^-0.0058

# wrangle phylogenetic tree
tree <- read.tree("data/gtdb/bac120.tree")
tree

# read in gtdb taxonomy
taxonomy <- read_tsv("data/gtdb/bac120_metadata.tsv", col_names = TRUE) %>%
  # rename columns
  select(accession, gtdb_taxonomy, ncbi_tax_id = ncbi_taxid) %>%
  filter(., accession %in% tree$tip.label)
nrow(taxonomy)

# get taxid from bacdive data
# load data
d_taxid <- readRDS('data/bacdive_all_info.rds') %>%
  select(bacdiveid, sequence_info) %>%
  unnest(sequence_info) %>%
  filter(!is.na(description)) %>%
  select(bacdiveid, ncbi_tax_id) %>%
  filter(bacdiveid %in% d_rod$bacdiveid) %>%
  distinct()
nrow(d_taxid)

d_taxid <- left_join(d_taxid, taxonomy, by = "ncbi_tax_id") %>%
  select(bacdiveid, ncbi_tax_id, gtdb_taxonomy, accession) %>%
  filter(!is.na(gtdb_taxonomy)) %>%
  distinct()

d_taxid <- left_join(
  d_taxid,
  select(d_joined, bacdiveid, biovolume, log_biovolume, mean_size, topt)
) %>%
  filter(!is.na(topt) & !is.na(mean_size) & !is.na(log_biovolume)) %>%
  distinct()

head(d_taxid)

tree_sub <- keep.tip(tree, d_taxid$accession)
tree_sub

# create d_meta
d_meta <- tibble(tip_label = tree_sub$tip.label) %>%
  left_join(., d_taxid, by = c("tip_label" = "accession")) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  separate(
    gtdb_taxonomy,
    c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";"
  ) %>%
  # get rid of p__ c__ etc
  mutate(across(kingdom:species, function(x) {
    gsub(".*__", "", x)
  }))
row.names(d_meta) <- d_meta$tip_label

d_meta %>%
  group_by(genus) %>%
  tally() %>%
  filter(n > 1) %>%
  arrange(desc(n))

filter(d_meta, genus == 'Prevotella') %>%
  ggplot(aes(log(biovolume), log(mean_size))) +
  geom_point()


group_by(d_meta, phylum) %>%
  tally() %>%
  arrange(desc(n))

d_meta <- d_meta %>%
  mutate(
    phylum2 = ifelse(
      phylum %in%
        c('Pseudomonadota', 'Bacteroidota', 'Bacillota', 'Actinomycetota'),
      phylum,
      "Other"
    )
  )

# group tip labels together in terms of their order
to_group <- split(d_meta$tip_label, d_meta$phylum2)
tree2 <- groupOTU(tree_sub, to_group)

# add columns for phyla
cols <- c(colorRampPalette(brewer.pal(4, "Spectral"))(4), 'grey')
names(cols) <- c(
  sort(c('Pseudomonadota', 'Bacteroidota', 'Bacillota', 'Actinomycetota')),
  'Other'
)


tree_base <- ggtree(
  tree2,
  aes(col = group),
  layout = "circular",
  branch.length = "none"
) +
  scale_color_manual('Phylum (branch colours)', values = cols) +
  guides(color = guide_legend(override.aes = list(linewidth = 3)))

tree_base %<+%
  d_meta +
  new_scale_color() +
  geom_tippoint(
    aes(x = x + x * 0.05, col = log10(mean_size)),
    stroke = NA,
    position = position_jitter(width = 1.5),
    show.legend = FALSE
  ) +
  theme_void(base_size = 20) +
  scale_color_viridis_c(
    name = "Genome size (bp)"
  )

ggsave(
  filename = "figures/phylogeny_genome_size.png",
  last_plot(),
  width = 12,
  height = 10
)

ggplot(d_meta, aes(log10(mean_size))) +
  geom_histogram(
    aes(fill = after_stat(x)),
    show.legend = FALSE,
    col = 'white'
  ) +
  theme_bw(base_size = 16) +
  labs(x = "Genome size (bp)") +
  scale_fill_viridis_c(name = "Genome size (bp)")


tree_base %<+%
  d_meta +
  new_scale_color() +
  geom_tippoint(
    aes(x = x + x * 0.05, col = log_biovolume),
    stroke = NA,
    position = position_jitter(width = 1.5),
    show.legend = FALSE
  ) +
  theme_void(base_size = 20) +
  scale_color_viridis_c(
    name = "Log biovolume (µm³)"
  )

ggsave(
  filename = "figures/phylogeny_body_size.png",
  last_plot(),
  width = 12,
  height = 10
)

ggplot(d_meta, aes(log_biovolume)) +
  geom_histogram(
    aes(fill = after_stat(x)),
    show.legend = FALSE,
    col = 'white'
  ) +
  theme_bw(base_size = 16) +
  labs(x = "Log biovolume (µm³)") +
  scale_fill_viridis_c(name = "Log biovolume ((µm³)")

tree_base %<+%
  d_meta +
  new_scale_color() +
  geom_tippoint(
    aes(x = x + x * 0.05, col = topt),
    stroke = NA,
    position = position_jitter(width = 1.5),
    show.legend = FALSE
  ) +
  theme_void(base_size = 20) +
  scale_color_viridis_c(
    name = "Optimum temperature (°C)"
  )

ggsave(
  filename = "figures/phylogeny_topt.png",
  last_plot(),
  width = 12,
  height = 10
)

ggplot(d_meta, aes(topt)) +
  geom_histogram(
    aes(fill = after_stat(x)),
    show.legend = FALSE,
    col = 'white',
    bins = 10
  ) +
  theme_bw(base_size = 16) +
  labs(x = "Optimum temperature (°C)") +
  scale_fill_viridis_c(name = "Optimum temperature (°C)")


# make a plot about random stuff
ggplot(d_meta, aes(topt, log_biovolume)) +
  geom_point(aes(col = phylum2), alpha = 0.5) +
  theme_bw(base_size = 16) +
  labs(
    x = "Optimum temperature (°C)",
    y = "Log biovolume (µm³)"
  ) +
  scale_color_manual(values = cols) +
  stat_smooth(method = 'lm', se = FALSE)
