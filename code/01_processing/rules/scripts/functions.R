# read_bedpe {{{ #
read_bedpe <- function (x) {
    cols <- c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'name', 'score', 'strand1', 'strand2')
    readr::read_tsv(x, col_names = cols, show_col_types = F)
}
# }}} read_bedpe #

# form_fragments {{{ #
form_fragments <- function (x, strand_as_R1) {
    n1 <- x %>%
        mutate(
            chr = chrom1,
            #strand = if_else(strand_as_R1 == 1, strand1, strand2),
            start = pmin(start1, start2),
            end = pmax(end1, end2)
        )
    if (strand_as_R1) {
        n2 <- n1 %>% mutate(strand = strand1)
    } else {
        n2 <- n1 %>% mutate(strand = strand2)
    }
    n3 <- n2 %>%
        select(chr, start, end, name, score, strand)
    return(n3)
}
# }}} form_fragments #

# count_frag {{{ #
count_frag <- function (fragments, reads_renamed) {
    n <- fragments %>%
        dtplyr::lazy_dt()
    if (reads_renamed) {
        m <- n %>%
           mutate(name = str_sub(str_extract(name, '_.+'), 2, -1)) %>%
           count(chr, start, end, strand, name)
    } else {
        m <- n %>%
           count(chr, start, end, strand) %>%
           mutate(name = '.')
    }
    counts_frag <- m %>%
        arrange(chr, start, end) %>%
        as_tibble() %>%
        rename(score = n) %>%
        select(chr, start, end, name, score, strand)
}
# }}} count_frag #

# vim:fdm=marker
