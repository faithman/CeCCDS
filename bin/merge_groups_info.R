#!/usr/bin/env Rscript
library(roperators)
library(tidyverse)
# Read data from Rscript input
args <- commandArgs(trailingOnly=TRUE)

df_isotype_group <- readr::read_tsv(args[1], col_names = TRUE)
#df_isotype_group <- readr::read_tsv("isotyep_groups.tsv", col_names = TRUE)

df_pairwise <- readr::read_tsv(args[2], col_names = TRUE)
#df_pairwise <- readr::read_tsv("merge_betweengroup_pairwise_output.tsv", col_names = c("pairwise", "strain_pairwise"))

df_npr1 <- readr::read_tsv(args[3], col_names = c("strain", "npr1_allele"))
#df_npr1 <- readr::read_tsv("npr1_allele_strain.tsv", col_names = c("strain", "npr1_allele"))

df_combine_1 <- df_isotype_group %>%
   dplyr::mutate(npr1_allele = ifelse(strain %in% df_npr1$strain, "TRUE", "FALSE"))


df_combine_2 <- df_pairwise %>%
   tidyr::separate(pairwise, into = c("a", "b"), sep = "-") %>%
   dplyr::filter(suspected_introgress == "YES") %>%
   dplyr::select(a,b,suspected_introgress) %>%
   tidyr::gather(type, strain, -suspected_introgress) %>%
   dplyr::group_by(strain) %>%
   dplyr::summarise(fail_pairwise = n())

df <- dplyr::left_join(df_combine_1, df_combine_2, by = "strain")

df$fail_pairwise %na<-% "FALSE"

# write out

readr::write_tsv(df, "new_isotype_groups.tsv", col_names = TRUE)