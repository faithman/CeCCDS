#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
# Read data from Rscript input
args <- commandArgs(trailingOnly=TRUE)

#df_raw <- readr::read_tsv("out_gt.tsv", col_names = c("CHROM", "POS", "A", "B"))
df_raw <- readr::read_tsv(args[3], col_names = TRUE)

df_raw_1 <- df_raw %>%
  dplyr::select(CHROM, POS, args[1], args[2]) #BRC20113_JU1530#AB4-ECA251

#df_raw_1 <- df_raw %>%
#  dplyr::select(CHROM, POS, BRC20113, JU1530)


names(df_raw_1) <- c("CHROM", "POS", "A", "B")
# Calculate discordant count
df <- df_raw_1 %>%
  dplyr::mutate(bin = round(POS/1E6)) %>%
  dplyr::filter(A != "./.", B != "./.", A != "0/1", B != "0/1") %>%
  dplyr::mutate(concordant = (A == B), discordant = (A != B)) %>%
  dplyr::group_by(CHROM,bin) %>%
  dplyr::summarize(count=n(),
                   discordant_count = sum(discordant),
                   concordant = sum(concordant)/count,
                   discordant = sum(discordant)/count) %>%
  dplyr::filter(CHROM != "MtDNA") %>%
  dplyr::mutate(type = ifelse(concordant > 0.99, "FALSE", "TRUE")) %>%
  dplyr::mutate(dis_type = ifelse(concordant < 0.98, "TRUE", "FALSE")) #%>%
  #dplyr::mutate(less_0.58 = ifelse(concordant > 0.02 & concordant < 0.58, "TRUE", "FALSE"))

# calculate bin counts
bin_count <- df %>% 
  dplyr::group_by(type) %>% 
  dplyr::summarise(count = n()) %>%
  dplyr::filter(type == "FALSE")


# accordant bin count
concordant_bin <- df %>%
  dplyr::ungroup() %>%
  dplyr::select(concordant, dis_type) %>%
  dplyr::filter(dis_type == "TRUE") %>%
  dplyr::mutate(bin = round(concordant*1e2)) %>%
  dplyr::group_by(bin) %>%
  dplyr::summarise(count = n())

# output dataframe for cutoff, accordant bin counts > 70; maximum disconcordant bin count <= 3, mean disconcordant bin count less than 2.
condition_results <- data.frame(pairwise_group = glue::glue("{args[1]}-{args[2]}"),
           accordant_bin_count = isTRUE(bin_count$count > 70),
           max_discordant_bin_count_lt_3 = isTRUE(max(concordant_bin$count) <= 5), # cutoff should be 3, 5 for test 
           mean_discordant_bin_count_lt_2.5 = isTRUE(mean(concordant_bin$count) <= 2),
           no_bin_lt_0.9 = isTRUE(min(concordant_bin$bin) >= 90))


if(nrow(bin_count) == 0) {
  print("No concordant bin")
  for_distribution <- data.frame(pairwise_group=character(),
                 accordant_bin_count=integer(), 
                 max_discordant_bin_count_lt_3=integer(), 
                 mean_discordant_bin_count_lt_2.5=integer(),
                 no_bin_lt_0.9=integer()) 
  readr::write_tsv(for_distribution, "for_distribution.tsv", col_names = FALSE)
} else if(nrow(concordant_bin) == 0) {
  print("no cordant bin")
  for_distribution <- data.frame(pairwise_group=character(),
                 accordant_bin_count=integer(), 
                 max_discordant_bin_count_lt_3=integer(), 
                 mean_discordant_bin_count_lt_2.5=integer(),
                 no_bin_lt_0.9=integer()) 
  readr::write_tsv(for_distribution, "for_distribution.tsv", col_names = FALSE)
} else if(is.infinite(max(concordant_bin$count))) {
  print("no cordant bin")
  for_distribution <- data.frame(pairwise_group=character(),
                 accordant_bin_count=integer(), 
                 max_discordant_bin_count_lt_3=integer(), 
                 mean_discordant_bin_count_lt_2.5=integer(),
                 no_bin_lt_0.9=integer()) 
  readr::write_tsv(for_distribution, "for_distribution.tsv", col_names = FALSE)
} else {
  for_distribution <- data.frame(pairwise_group = glue::glue("{args[1]}-{args[2]}"),
           accordant_bin_count = bin_count$count,
           max_discordant_bin_count_lt_3 = max(concordant_bin$count, na.rm= TRUE),
           mean_discordant_bin_count_lt_2.5 = mean(concordant_bin$count, na.rm= TRUE),
           no_bin_lt_0.9 = min(concordant_bin$bin, na.rm= TRUE))
  readr::write_tsv(for_distribution, "for_distribution.tsv", col_names = FALSE)
}


condition_results <- condition_results %>%
  dplyr::mutate(suspected_introgress = ifelse(accordant_bin_count == "TRUE" & max_discordant_bin_count_lt_3 == "TRUE" & mean_discordant_bin_count_lt_2.5 == "TRUE",
    "YES", "NO"))

readr::write_tsv(condition_results, "condition_results.tsv", col_names = FALSE)


## Plot pairwise disconcordance across Chromosomes
if(nrow(bin_count) == 0) {
  print("No concordant bin")
} else if(condition_results$suspected_introgress == "NO") {
  print("concordant bins smaller than cutoff")
} else if(condition_results$suspected_introgress == "YES"){
  ggplot(df, aes(x = bin, y = discordant)) + 
    geom_bar(aes(fill = discordant >= 0.02), stat='identity') +
    scale_fill_manual(values = c("gray", "red")) +
    facet_grid(. ~ CHROM, scales = "free_x", space = "free_x") +
    theme_bw() + 
    #scale_y_continuous(limits = c(0, 0.10)) +
    labs(x = "Position", y = "Discordant SNPs (%)", title = glue::glue("{args[1]} VS {args[2]}"))
  ggsave(glue::glue("{args[1]}-{args[2]}.disconcordance.png"), width = 12, height = 3)#
  
  ## Plot histogram for concordance distribution
  ggplot(df) + 
    aes(x = concordant) +
    geom_histogram(aes(fill = discordant >= 0.02), binwidth = 0.01) +
    scale_fill_manual(values = c("#808080", "#0080FF")) +
    labs(x = "Condordance", y = "count", title = glue::glue("{args[1]} VS {args[2]}"))
  ggsave(glue::glue("{args[1]}-{args[2]}.hist.png"), width = 12, height = 3)
}
