# read_counts_frag_bed {{{ #
read_counts_frag_bed <- function (x) {
    arrow::read_tsv_arrow(
        x,
        col_names = c('chr', 'start', 'end', 'type',
                      'count_frag', 'strand', 'sample'),
        as_data_frame = TRUE
    )
}
# }}} read_counts_frag_bed #

# change_coor_from_0_base_to_1_base {{{ #
change_coor_from_0_base_to_1_base <- function (x) {
    x %>%
        mutate(end5 = ifelse(strand == '+', start + 1, end)) %>%
        mutate(end3 = ifelse(strand == '+', end, start + 1)) %>%
        select(sample, chr, end5, end3, strand, type, count_frag)
}
# }}} change_coor_from_0_base_to_1_base #

# count_end5 {{{ #
count_end5 <- function (x) {
    x %>%
        group_by(sample, chr, strand, type, end5) %>%
        summarise(count_end5 = sum(count_frag), .groups = 'drop')
}
# }}} count_end5 #

# calculate_cpm {{{ #
calculate_cpm <- function (count, count_sample, fake_count) {
    cpm <- (count + fake_count) * 1e6 / count_sample
    return(cpm)
}
# }}} calculate_cpm #

# filter_region_genome {{{ #
filter_region_genome <- function (counts_list, region) {
    pmap(
        region |> select(starts_with('region'), chr, subset, strand, library_method),
        function(region_id, region_start, region_end, chr, subset, strand, library_method, counts_list) {
            print(region_id)
            xregion_id <- region_id
            xregion_start <- region_start
            xregion_end <- region_end
            xchr <- chr
            xsubset <- subset
            xstrand <- strand
            xlibrary_method <- library_method
            list_name <- str_c(xsubset, xlibrary_method, xchr, xstrand, sep = '%')
            counts_list[[list_name]] |>
                filter(between(end5, xregion_start, xregion_end)) |>
                mutate(region_id = xregion_id)
        },
        counts_list = counts_list
    )
}
# }}} filter_region_genome #

# calculate_ratio {{{ #

calculate_ratio <- function (counts_region, offset, upstream_only) {
   n <- counts_region |>
        mutate(rel_coor = ifelse(strand == '+', end5 - peak, peak - end5))
    if (upstream_only) {
        n <- n |>
            filter(rel_coor <= 0)
    }
    n <- n %>%
        filter(between(abs(end5 - peak), 0, offset + 2)) %>%
        filter(!between(abs(end5 - peak), 1, 2)) %>%
        mutate(tag = ifelse(end5 == peak, 'cs', 'other')) %>%
        group_by(peak_id, sample, chr, strand, peak, tag) %>%
        summarise(count_end5_sum = sum(count_end5), .groups = 'drop') %>%
        pivot_wider(names_from = tag, values_from = count_end5_sum, values_fill = 0)
    if (!('cs' %in% names(n))) {
        n <- n %>% mutate(cs = 0)
    }
    if (!('other' %in% names(n))) {
        n <- n %>% mutate(other = 0)
    }
    n %>%
        mutate(ratio = cs / (1 + other))
}

# }}} calculate_ratio_co #

# calculate_ratio_genome {{{ #
calculate_ratio_genome <- function (counts_list, peaks, offset) {
    n <- lst()
    regions <- peaks %>%
        mutate(
            region_start = peak - offset - 2,
            region_end = peak + offset + 2,
            region_id = str_c(peak_id, peak, sep = '%')
        )
    counts_region <- filter_region_genome(counts_list, regions)
    m <- list_rbind(counts_region) |>
        separate(region_id, into = c('peak_id', 'peak'), sep = '%', convert = TRUE)
    list_str <- lst(
        offset = c(100),
        upstream_only = c(TRUE)
    )
    list_str_df <- expand_grid(!!!list_str) |>
        mutate(tag = paste(ifelse(upstream_only, 'u', 'ud'), offset, sep = '_'))
    n <- pmap(
        list_str_df,
        function(offset, upstream_only, tag) {
            calculate_ratio(m, offset = offset, upstream_only = upstream_only)
        }
    )
    n <- set_names(n, list_str_df$tag)
    return(n)
}
# }}} calculate_ratio_genome #

# vim:fdm=marker
