# library ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

library('arrow')
library('tidyverse')
library('dtplyr')
library('readxl')
library('furrr')
plan(multicore, workers = 6)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

source('./function.r')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

ratio_offset_max <- 100

relative_path <- "../../../processing"

library_method_recode2 <- c(L1 = 'ligation_oligo1', L2 = 'ligation_oligo2', TA = 'tailing_polyA', TC = 'tailing_polyC')

library_method_recode <- names(library_method_recode2) %>% setNames(library_method_recode2)


sample_info <- read_xlsx("../data/sample_info.xlsx") |> filter(!(!is.na(exclude) & exclude == 'Y'))


sample_info_recode <- sample_info |>
    mutate(
        sample_rep = case_match(sample_rep, 'rep1' ~ 'r1', 'rep2' ~ 'r2', 'rep3' ~ 'r3'),
        library_method = recode(library_method, !!!library_method_recode),
        sample_recode = str_c(species, library_method, sample_rep, sample_type, sep = '|')
    ) |>
    select(species, sample, sample_recode, library_method, sample_rep, sample_type, dataset)

sample_recode_vector <- sample_info_recode$sample_recode |> setNames(sample_info_recode$sample)


sample_info_path <- sample_info_recode %>%
    mutate(path = file.path(relative_path, str_replace(dataset, '#.+', ''), 'results', '04.bed', paste0(sample, '_genome.bed')))

sample_info_path_spikeins <- sample_info_recode %>%
    mutate(path = file.path(relative_path, str_replace(dataset, '#.+', ''), 'results', '04.bed', paste0(sample, '_spikeins.bed')))


# chr recode {{{ #

genome_path <- '../../../biodata/genomes'

chr_recode_GRCm39 <- read_tsv(file.path(genome_path, "GRCm39/mm39.chromAlias.txt"), comment = '#', col_names = c('UCSC', 'ensembl', 'ENA', 'genbank'))


chr_recode_GRCz11 <- read_tsv(file.path(genome_path, "GRCz11/GRCz11_no_alt/danRer11.chromAlias.txt"), comment = '#', col_names = c('UCSC', 'ensembl', 'ENA', 'genbank')) %>%
    mutate(genbank = ifelse(is.na(genbank), ENA, genbank)) %>%
    filter(!str_detect(UCSC, 'alt$'))


chr_recode_IRGSP <- read_tsv(file.path(genome_path, "IRGSP-1.0/GCF_001433935.1_IRGSP-1.0_genomic_chr_main.tsv"), col_names = c('genbank', 'UCSC'))


chr_recode <- chr_recode_GRCz11 %>%
    bind_rows(chr_recode_GRCm39) %>%
    select(UCSC, genbank) %>%
    bind_rows(chr_recode_IRGSP)

chr_recode_vector <- chr_recode$UCSC %>% setNames(chr_recode$genbank)

chr_recode_vector2_GRCm39 <- chr_recode_GRCm39$genbank %>% setNames(chr_recode_GRCm39$UCSC)

chr_recode_vector2_GRCz11 <- chr_recode_GRCz11$genbank %>% setNames(chr_recode_GRCz11$UCSC)

chr_recode_vector2_IRGSP <- chr_recode_IRGSP$genbank %>% setNames(chr_recode_IRGSP$UCSC)


chr_recode_main_chr <- chr_recode_GRCz11 %>%
    mutate(species = 'zebrafish') |>
    bind_rows(chr_recode_GRCm39 |> mutate(species = 'mouse')) %>%
    select(species, UCSC, genbank) %>%
    bind_rows(chr_recode_IRGSP |> mutate(species = 'rice')) |>
    filter(str_detect(UCSC, 'chr\\d+$') | str_detect(UCSC, 'chr[XY]$'))


# }}} chr recode #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{


# pre process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

# counts frag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

counts_frag <- map(
    sample_info_path$path,
    read_counts_frag_bed
)

counts_frag_spikeins <- map(
    sample_info_path_spikeins$path,
    read_counts_frag_bed
)

# counts frag ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# sample counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

sample_counts <- map(
    counts_frag,
    \(x) summarise(x, count_frag_sample = sum(count_frag), .by = sample)
) |> list_rbind()


sample_counts_spikeins <- map(
    counts_frag_spikeins,
    \(x) summarise(x, count_frag_sample = sum(count_frag), .by = sample)
) |> list_rbind()

sample_counts_all <- bind_rows(sample_counts, sample_counts_spikeins) %>%
    summarise(count_frag_sample = sum(count_frag_sample), .by = sample) |>
    mutate(sample = recode(sample, !!!sample_recode_vector))

write_tsv(sample_counts_all, "../data/sample_counts_all.tsv")

sample_counts_all <- read_tsv("../data/sample_counts_all.tsv")

# sample counts ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# counts end5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

counts_end5_uncapped_info <- map(
    counts_frag,
    \(x, sample_info) {
        x |>
            dtplyr::lazy_dt() %>%
            filter(type != 'capped', chr %in% chr_recode_main_chr$genbank) |>
            change_coor_from_0_base_to_1_base() |>
            count_end5() |>
            left_join(sample_info) |>
            collect()
    },
    sample_info = sample_info_recode |> select(-dataset) |> distinct()
) |> list_rbind()


counts_end5_uncapped_info_merged <- counts_end5_uncapped_info |>
    dtplyr::lazy_dt() |>
    summarise(count_end5 = sum(count_end5), .by = c(species, sample_recode, library_method, sample_rep, sample_type, chr, strand, end5)) |>
    dplyr::rename(sample = sample_recode) |>
    collect()


write_feather(counts_end5_uncapped_info_merged, "../results/counts_end5_uncapped_info_merged.feather")


counts_end5_uncapped_info_merged |>
    group_by(species, library_method, chr, strand) |>
    write_dataset(
        "../results/counts_end5_uncapped_info",
        format = 'feather'
    )


counts_end5_uncapped_info <- open_dataset("../results/counts_end5_uncapped_info/", format = 'feather')

list_str_df <- chr_recode_main_chr |>
    select(-UCSC) |>
    mutate(library_method = list(library_method_recode)) |>
    mutate(strand = list(c('+', '-'))) |>
    unnest(cols = c(strand)) |>
    unnest(cols = c(library_method)) |>
    mutate(name = str_c(species, library_method, genbank, strand, sep = '%'))

counts_end5_uncapped_info_list_o_flatten <- pmap(
    list_str_df,
    \(species, genbank, library_method, strand, name) {
        xspecies <- species
        xchr <- genbank
        xlibrary_method <- library_method
        xstrand <- strand
        counts_end5_uncapped_info %>%
            filter(species == xspecies, library_method == xlibrary_method, chr == xchr, strand == xstrand)
    }
)

names(counts_end5_uncapped_info_list_o_flatten) <- list_str_df$name

# counts end5 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# pre process ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# peaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

sample_rep_num <- sample_counts_all |>
    separate(sample, into = c('species', 'library_method', 'sample_rep', 'sample_type'), sep = '\\|') |> 
    select(species, library_method, sample_rep) |>
    distinct() |>
    count(species, library_method)


sample_type_num <- sample_counts_all |>
    separate(sample, into = c('species', 'library_method', 'sample_rep', 'sample_type'), sep = '\\|') |> 
    select(species, library_method, sample_type) |>
    distinct() |>
    count(species, library_method)


list_names <- names(counts_end5_uncapped_info_list_o_flatten)

all_species <- c('mouse', 'zebrafish', 'rice')

walk(
    c('mouse', 'zebrafish', 'rice'),
    \(x) {
        counts_end5_uncapped_info_list_o_flatten3 <- counts_end5_uncapped_info_list_o_flatten[str_detect(list_names, x)]
        counts_end5_uncapped_info_hit_all <- future_imap(
            counts_end5_uncapped_info_list_o_flatten3,
            function(.x, .y, sample_rep_num, sample_counts) {
                print(.y)
                info <- str_split_1(.y, '%')
                n <- pull(sample_rep_num |> filter(species == info[1], library_method == info[2]), n)
                .x |>
                    lazy_dt() |>
                    left_join(sample_counts) |>
                    mutate(cpm_end5 = calculate_cpm(count_end5, count_frag_sample, 0)) |>
                    group_by(chr, end5, sample_rep) |>
                    summarise(cpm_end5_sum = sum(cpm_end5)) |>
                    filter(n() >= n) |>
                    collect()
            },
            sample_rep_num = sample_rep_num,
            sample_counts = sample_counts_all
        )
        counts_end5_uncapped_info_hit_all2 <- list_rbind(counts_end5_uncapped_info_hit_all, names_to = 'info') |>
            separate(info, into = c('species', 'library_method', 'chr', 'strand'), sep = '%', remove = TRUE)
        write_feather(counts_end5_uncapped_info_hit_all2, glue::glue("../analysis/counts_end5_uncapped_info_hit_all_{x}.feather"))
    }
)


counts_end5_uncapped_info_hit_all_all <- list.files(
        "../analysis/",
        'counts_end5_uncapped_info_hit_all',
        full.names = TRUE
    ) |>
    map(read_feather) |> list_rbind()


peaks_all_bed <- counts_end5_uncapped_info_hit_all_all |>
    summarise(cpm_end5_sum = sum(cpm_end5_sum), .by = c(species, library_method, chr, strand, end5)) |>
    left_join(sample_type_num) |>
    dplyr::rename(n_type = n) |>
    left_join(sample_rep_num) |>
    mutate(cpm_end5_mean = cpm_end5_sum / (n_type * n)) |>
    mutate(chrom = str_c(species, library_method, chr, sep = '|')) |>
    mutate(start = end5 - 1, end = end5) |>
    mutate(id = str_c(species, chr, strand, end5, library_method, sep = '|')) |>
    select(chrom, start, end, id, cpm_end5_mean, strand) |>
    arrange(chrom, start)


write_tsv(peaks_all_bed, "../analysis/peaks_all.bed", col_names = FALSE)


peaks_all_cpm_mean <- peaks_all_bed |>
    select(id, cpm_end5_mean)


system("bedtools merge -s -d 5 -c 6,4 -o distinct -i ../analysis/peaks_all.bed > ../analysis/peaks_all_merged.bed")

peaks_all_merged <- read_tsv(
    "../analysis/peaks_all_merged.bed",
    col_names = c('chr', 'start', 'end', 'strand', 'id'),
    col_types = 'ciicc'
)

peaks_all_region_max <- peaks_all_merged |>
    mutate(id = str_split(id, ',')) |>
    unnest(id) |>
    lazy_dt() |>
    left_join(peaks_all_cpm_mean) |>
    group_by(chr, start, end, strand) |>
    slice_max(cpm_end5_mean, with_ties = FALSE) |>
    select(chr, start, end, id, cpm_end5_mean, strand) |>
    ungroup() |>
    collect()

write_tsv(peaks_all_region_max, "../analysis/peaks_all_region_max.bed", col_names = FALSE)


peaks_lm <- peaks_all_region_max |>
    select(id, cpm_end5_mean) %>%
    separate(id, into = c('species', 'chr', 'strand', 'end5', 'library_method'), remove = FALSE, sep = '\\|', convert = TRUE) %>%
    split(.$library_method) |>
    future_map(
        ~ .x |>
            mutate(start = end5 - 1, end = end5) |>
            select(chr, start, end, id, cpm_end5_mean, strand) |>
            arrange(chr, start)
    )


iwalk(
    peaks_lm,
    ~ write_tsv(.x, glue::glue("../analysis/peaks_{.y}.bed"), col_names = FALSE)
)

system("bash ./get_proceeding_seq_for_peaks.sh")

system("bash ./calculate_seq_identity_for_peaks.sh")

peak_seq_iden <- list.files("../analysis", 'iden', full.names = TRUE) |>
    read_tsv(col_names = c('oligo', 'id', 'identical_base_num')) |>
    mutate(id = str_replace(id, '\\([+-]\\)', '')) |>
    mutate(oligo_len = str_split_i(oligo, '_', 3))


peak_seq_iden_list <- peak_seq_iden %>%
    split(.$oligo_len)


peak_seq_iden_filtered_10 <- peak_seq_iden_list[['10']] |>
    filter(identical_base_num < 8)


peak_seq_iden_filtered_4 <- peak_seq_iden_list[['4']] |>
    filter(identical_base_num < 4)


peaks_seq_iden_filtered <- peaks_all_region_max |>
    filter(id %in% intersect(peak_seq_iden_filtered_4$id, peak_seq_iden_filtered_10$id))


peaks_seq_iden_filtered_lm <- peaks_seq_iden_filtered |>
    select(id, cpm_end5_mean) %>%
    separate(id, into = c('species', 'chr', 'strand', 'end5', 'library_method'), remove = FALSE, sep = '\\|', convert = TRUE) %>%
    split(.$library_method) |>
    map(
        ~ .x |>
            mutate(start = end5 - 1, end = end5) |>
            select(chr, start, end, id, cpm_end5_mean, strand) |>
            arrange(chr, start)
    )


iwalk(
    peaks_seq_iden_filtered_lm,
    ~ write_tsv(.x, glue::glue("../analysis/peaks_seq_iden_filtered_{.y}.bed"), col_names = FALSE)
)


# peaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# merge peaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

system("cat ../analysis/peaks_seq_iden_filtered_TA.bed ../analysis/peaks_seq_iden_filtered_TC.bed | bedtools sort -i | bedtools merge -s -d 5 -c 6,4 -o distinct -i -  > ../analysis/peaks_tailing_merged.bed")

peaks_tailing_merged <- read_tsv(
    "../analysis/peaks_tailing_merged.bed",
    col_names = c('chr', 'start', 'end', 'strand', 'id'),
    col_types = 'ciicc'
)

peaks_tailing_merged_list <- peaks_tailing_merged %>%
    split((str_detect(.$id, 'TA') & str_detect(.$id, 'TC')))


peaks_tailing_combined <- peaks_tailing_merged_list[['TRUE']] %>%
    mutate(id = str_split(id, ',')) %>%
    unnest(id) %>%
    mutate(library_method = str_extract(id, '[^\\|]+$')) %>%
    left_join(peaks_all_cpm_mean) %>%
    group_by(chr, start, end, strand, library_method) %>%
    slice_max(cpm_end5_mean, with_ties = FALSE) %>%
    ungroup()


peaks_tailing_combined_bed <- peaks_tailing_combined %>%
    mutate(end5 = as.numeric(str_split_i(id, '\\|', 4))) %>%
    group_by(chr, start, end, strand) %>%
    mutate(start = min(end5) - 1, end = max(end5)) %>%
    ungroup() %>%
    select(chr, start, end, id, library_method, strand) %>%
    pivot_wider(names_from = library_method, values_from = id) %>%
    mutate(id = str_c(TA, TC, sep = ','), score = '.') %>%
    select(chr, start, end, id, score, strand) %>%
    arrange(chr, start)


write_tsv(peaks_tailing_combined_bed, "../analysis/peaks_tailing_combined.bed", col_names = FALSE)


peaks_tailing_single_bed <- peaks_tailing_merged_list[['FALSE']] |>
    mutate(score = '.') |>
    select(chr, start, end, id, score, strand) %>%
    arrange(chr, start)

write_tsv(peaks_tailing_single_bed, "../analysis/peaks_tailing_single.bed", col_names = FALSE)

system("cat ../analysis/peaks_seq_iden_filtered_L1.bed ../analysis/peaks_seq_iden_filtered_L2.bed ../analysis/peaks_tailing_combined.bed ../analysis/peaks_tailing_single.bed | bedtools sort -i | bedtools merge -s -d 0 -c 6,4 -o distinct -i - > ../analysis/peaks_merged.bed")


peaks_merged <- read_tsv(
    "../analysis/peaks_merged.bed",
    col_names = c('chr', 'start', 'end', 'strand', 'id'),
    col_types = 'ciicc'
)


peaks_merged_long <- peaks_merged %>%
    mutate(id = str_split(id, ',')) %>%
    unnest(id) %>%
    mutate(library_method = str_extract(id, '[^\\|]+$')) %>%
    mutate(end5 = as.numeric(str_split_i(id, '\\|', 4))) %>%
    lazy_dt() %>%
    left_join(peaks_all_cpm_mean) %>%
    group_by(chr, start, end, strand, library_method) %>%
    slice_max(cpm_end5_mean, with_ties = FALSE) %>%
    ungroup() %>%
    collect()

write_tsv(peaks_merged_long, "../analysis/peaks_merged_long.tsv")

peaks_merged_long <- read_tsv("../analysis/peaks_merged_long.tsv")


# merge peaks ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# tier 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

peaks_t0 <- peaks_merged_long %>%
    group_by(chr, start, end, strand) %>%
    filter(n()>3) %>%
    mutate(start = min(end5 - 1), end = max(end5)) %>%
    ungroup() %>%
    dplyr::rename(peak = end5) %>%
    mutate(species = str_split_i(id, "\\|", 1)) %>%
    mutate(peak_id = str_c(species, chr, start + 1, strand, sep = '|'))

peaks_t0_site <- peaks_t0 %>%
    dplyr::rename(subset = species)

peaks_t0_ratio <- calculate_ratio_genome(counts_end5_uncapped_info_list_o_flatten, peaks_t0_site, ratio_offset_max)


peaks_t0_ratio_all <- peaks_t0_ratio %>%
    list_rbind(names_to = 'ratio_type')

write_tsv(peaks_t0_ratio_all, "../analysis/peaks_t0_ratio_all.tsv")

peaks_t0_ratio_all <- read_tsv("../analysis/peaks_t0_ratio_all.tsv")


# tier 0 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# tier 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{


peaks_t1 <- peaks_merged_long %>%
    group_by(chr, start, end, strand) %>%
    filter(n() == 3) %>%
    mutate(start = min(end5 - 1), end = max(end5)) %>%
    ungroup() %>%
    dplyr::rename(peak = end5) %>%
    mutate(species = str_split_i(id, "\\|", 1)) %>%
    mutate(peak_id = str_c(species, chr, start + 1, strand, sep = '|'))


peaks_t1_site <- peaks_t1 %>%
    dplyr::rename(subset = species)


peaks_t1_ratio <- calculate_ratio_genome(counts_end5_uncapped_info_list_o_flatten, peaks_t1_site, ratio_offset_max)


peaks_t1_ratio_all <- peaks_t1_ratio %>%
    list_rbind(names_to = 'ratio_type')


write_tsv(peaks_t1_ratio_all, "../analysis/peaks_t1_ratio_all.tsv")


# tier 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# tier 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{


peaks_t2 <- peaks_merged_long %>%
    group_by(chr, start, end, strand) %>%
    filter(n() == 2) %>%
    mutate(start = min(end5 - 1), end = max(end5)) %>%
    ungroup() %>%
    dplyr::rename(peak = end5) %>%
    mutate(species = str_split_i(id, "\\|", 1)) %>%
    mutate(peak_id = str_c(species, chr, start + 1, strand, sep = '|'))


peaks_t2_site <- peaks_t2 %>%
    dplyr::rename(subset = species)

peaks_t2_ratio <- calculate_ratio_genome(counts_end5_uncapped_info_list_o_flatten, peaks_t2_site, ratio_offset_max)

peaks_t2_ratio_all <- peaks_t2_ratio %>%
    list_rbind(names_to = 'ratio_type')

write_tsv(peaks_t2_ratio_all, "../analysis/peaks_t2_ratio_all.tsv")


# tier 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# calculate measures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{


peaks_tier_all_list <- lst(
    't0' = peaks_t0_ratio_all,
    't1' = peaks_t1_ratio_all,
    't2' = peaks_t2_ratio_all
)


peaks_tier_all <- list_rbind(peaks_tier_all_list, names_to = 'tier')


# peaks distribution on chr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

genome_len <- read_tsv("../../../biodata/genomes/combined_GRCm39_GRCz11_IRGSP-1.0/combined_GRCm39_GRCz11_IRGSP-1.0.genome", col_names = c('chr', 'length'))


peaks_ids <- peaks_tier_all |>
    select(tier, peak_id, sample) |>
    distinct() |>
    separate(sample, into = c('species', 'library_method', 'sample_rep', 'sample_type'), remove = FALSE, sep = '\\|') |>
    select(tier, peak_id, library_method) |>
    distinct()

peaks_ids_coor <- peaks_ids |>
    separate(peak_id, into = c('species', 'chr', 'end5', 'strand'), convert = TRUE, remove = FALSE, sep = '\\|') |>
    left_join(genome_len) |>
    mutate(coor_per = end5 / length)


p_peaks_distribution_mouse <- (
    ggplot(peaks_ids_coor |> filter(species == 'mouse'), aes(coor_per, color = library_method))
    + stat_ecdf()
    + facet_wrap(vars(chr))
)

ggsave("../analysis/plots/peaks_distribution_mouse.pdf", p_peaks_distribution_mouse, width = 525, height = 300, units = 'mm')


peaks_ids_coor_mouse_chr11 <- peaks_ids_coor |>
    filter(species == 'mouse') |>
    filter(chr == 'NC_000077.7') |>
    arrange(coor_per)

xtabs(~ library_method + coor_per, peaks_ids_coor_mouse_chr11)
# < 0.03

peaks_ids_coor_mouse_chr17 <- peaks_ids_coor |>
    filter(species == 'mouse') |>
    filter(chr == 'NC_000083.7') |>
    arrange(coor_per)

xtabs(~ library_method + coor_per, peaks_ids_coor_mouse_chr17)
# > 0.42, < 0.43


p_peaks_distribution_zebrafish <- (
    ggplot(peaks_ids_coor |> filter(species == 'zebrafish'), aes(coor_per, color = library_method))
    + stat_ecdf()
    + facet_wrap(vars(chr))
)

ggsave("../analysis/plots/peaks_distribution_zebrafish.pdf", p_peaks_distribution_zebrafish, width = 525, height = 300, units = 'mm')


p_peaks_distribution_rice <- (
    ggplot(peaks_ids_coor |> filter(species == 'rice'), aes(coor_per, color = library_method))
    + stat_ecdf()
    + facet_wrap(vars(chr))
)

ggsave("../analysis/plots/peaks_distribution_rice.pdf", p_peaks_distribution_rice, width = 525, height = 300, units = 'mm')


peaks_mouse_chr11_to_exclude <- peaks_ids_coor |>
    filter(species == 'mouse', chr == 'NC_000077.7') |>
    filter(coor_per < 0.03)

peaks_mouse_chr17_to_exclude <- peaks_ids_coor |>
    filter(species == 'mouse', chr == 'NC_000083.7') |>
    filter(between(coor_per, 0.42, 0.43))


peaks_chr_to_exclude <- bind_rows(
    peaks_mouse_chr11_to_exclude,
    peaks_mouse_chr17_to_exclude
) |> select(peak_id) |> distinct()


# peaks distribution on chr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


peaks_tier_all_filtered <- peaks_tier_all |>
    filter(!peak_id %in% peaks_chr_to_exclude$peak_id) |>
    separate(sample, into = c('species', 'library_method', 'sample_rep', 'sample_type'), remove = FALSE) |>
    left_join(sample_type_num %>% dplyr::rename(n_type = n)) |>
    left_join(sample_rep_num %>% dplyr::rename(n_rep = n))

write_tsv(peaks_tier_all_filtered, "../analysis/peaks_tier_all_filtered.tsv")

peaks_tier_ratio <- peaks_tier_all_filtered |>
    select(tier, peak_id, sample, chr, strand, species, library_method, sample_rep, n_type, n_rep, ratio, ratio_type) |>
    distinct() |>
    group_by(tier, peak_id, ratio_type, library_method, sample_rep, n_rep) |>
    summarise(ratio_mean = sum(ratio) / unique(n_type), .groups = 'drop') |>
    group_by(tier, peak_id, ratio_type, library_method) |>
    summarise(ratio_mean_mean = sum(ratio_mean) / unique(n_rep)) |>
    group_by(tier, peak_id, ratio_type) |>
    summarise(ratio_mean_mean_mean = sum(ratio_mean_mean) / 4, .groups = 'drop') |>
    mutate(ratio_type = str_c('ratio', ratio_type, sep = '_')) |>
    pivot_wider(id_cols = tier:peak_id, names_from = ratio_type, values_from = c(ratio_mean_mean_mean))


peaks_tier_peak <- peaks_tier_all_filtered |>
    select(tier, peak_id, species, chr, strand, peak, library_method) |>
    distinct() |>
    pivot_wider(id_cols = tier:strand, names_from = library_method, values_from = peak, names_prefix = 'peak_') |>
    mutate(chr_UCSC = recode(chr, !!!chr_recode_vector)) |>
    relocate(chr_UCSC, .after = chr) |>
    select(tier:strand, peak_L1, peak_L2, peak_TA, peak_TC)


peaks_tier_info_all <- peaks_tier_peak |>
    left_join(peaks_tier_ratio |> select(peak_id, starts_with('ratio')))

write_tsv(peaks_tier_info_all, "../analysis/peaks_tier_info_all.tsv")


# calculate measures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# spikeins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{{{

counts_end5_spikeins_uncapped_info <- map(
    counts_frag_spikeins,
    \(x, sample_info) {
        x |>
            dtplyr::lazy_dt() %>%
            filter(type != 'capped') |>
            change_coor_from_0_base_to_1_base() |>
            count_end5() |>
            left_join(sample_info) |>
            collect()
    },
    sample_info = sample_info_recode |> select(-dataset) |> distinct()
) |> list_rbind()



counts_end5_spikeins_uncapped_info_merged <- counts_end5_spikeins_uncapped_info |>
    dtplyr::lazy_dt() |>
    summarise(count_end5 = sum(count_end5), .by = c(species, sample_recode, library_method, sample_rep, sample_type, chr, strand, end5)) |>
    dplyr::rename(sample = sample_recode) |>
    collect()

writexl::write_xlsx(counts_end5_spikeins_uncapped_info_merged, "../analysis/spikeins_counts.xlsx")


# spikeins ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}}}

# vim:fdm=marker
