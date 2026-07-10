# library ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

library('readr', warn.conflicts = F)
library('dplyr', warn.conflicts = F)
library('stringr', warn.conflicts = F)
library('dtplyr', warn.conflicts = F)
snakemake@source('./functions.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# setting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

# snakemake

## input
bedpe_file <- snakemake@input[['bedpe']]

## output
output_file <- snakemake@output[['counts']]

## parameters
reads_renamed <- snakemake@params[['reads_renamed']]
strand_as_R1 <- snakemake@params[['strand_as_R1']]
sample <- snakemake@params[['sample']]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

reads <- read_bedpe(bedpe_file)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

fragments <- form_fragments(reads, strand_as_R1)

counts <- count_frag(fragments, reads_renamed)

counts <- counts %>%
    mutate(sample = sample)

write_tsv(counts, output_file, col_names = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# vim:fdm=marker
