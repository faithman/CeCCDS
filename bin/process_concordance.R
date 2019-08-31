#!/usr/bin/env Rscript
library(ggplot2)
library(tidyverse)
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
# Used in calculating isotypes
stack_list <- function(x) {
  if (is.null(names(x)) == T) {
    names(x) <- as.character(1:length(x))
  }
  stack(x)
}

args = commandArgs(trailingOnly=TRUE)
# Debug is enabled
if (args == "true") {
    coverage_level = 0
} else {
    coverage_level = 20
}


coverage_20 <- (readr::read_tsv("SM_coverage.tsv") %>%
                  dplyr::filter(coverage > coverage_level))$strain

WI <- readr::read_tsv("https://docs.google.com/spreadsheets/d/1V6YHzblaDph01sFDI8YK_fP0H7sVebHQTXypGdiQIjI/pub?output=tsv") 

WI %>% readr::write_tsv("WI_metadata.tsv")

existing_WI <- WI %>%
  dplyr::select(strain, isotype, latitude, longitude) %>%
  dplyr::filter(strain %in% coverage_20) %>%
  dplyr::filter(isotype != "NA")

f <- file("isotype_count.txt")

cutoff <- 0.999

gtcheck <- readr::read_tsv("gtcheck.tsv") %>%
  dplyr::filter((i %in% coverage_20) & (j %in% coverage_20)) %>%
  dplyr::mutate(concordance = 1-(discordance/sites)) %>%
  dplyr::mutate(isotype = concordance > cutoff) %>%
  dplyr::filter(!(i %in% c("LSJ1", "JU2250")) & !(j %in% c("LSJ1", "JU2250")))

# Generate strains that do not group with any other strains (single strains)
single_strains <- gtcheck %>%
  dplyr::select(i, j, isotype) %>%
  tidyr::gather(col, strain, -isotype) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(single_strain = (sum(isotype, na.rm = T) == 0)) %>%
  dplyr::distinct(.keep_all = T) %>%
  dplyr::select(strain, single_strain) %>% 
  dplyr::filter(single_strain)

# Add LSJ1
single_strains <- dplyr::bind_rows(list(strain = "LSJ1", single_strain = T),
                                   list(strain = "JU2250", single_strain = T),
                                   single_strains)

single_strains <- as.data.frame(list(i = single_strains$strain,
                                     j = single_strains$strain,
                                     discordance = 0,
                                     concordance = 1,
                                     isotype = TRUE))

gtcheck <- dplyr::bind_rows(gtcheck, single_strains)

# Filter for only grouped strains
iso_gtcheck <- gtcheck %>% dplyr::filter(isotype == T)

# Generate complete strain list
strain_list <- sort(unique(c(gtcheck$i, gtcheck$j)))

# Generate isotype groups
isotype_groups <- stack_list(unique(lapply(strain_list, function(x) {
  grouped_strains <- dplyr::filter(iso_gtcheck, (i == x | j == x)) %>%
    dplyr::select(i, j)
  sort(unique(unlist(grouped_strains)))
}))) %>%
  dplyr::rename(strain = values, group = ind) %>%
  dplyr::mutate(group = ifelse(strain == "LSJ1", 0, group)) %>%
  dplyr::mutate(group = ifelse(strain == "JU2250", -1, group)) %>%
  dplyr::distinct() %>% 
  dplyr::group_by(group) %>%
  tidyr::nest(strain) %>%
  dplyr::distinct(data, .keep_all = T) %>%
  tidyr::unnest()

SM_coverage <- readr::read_tsv("SM_coverage.tsv")

isotype_groups <- dplyr::left_join(isotype_groups, existing_WI, by = c("strain")) %>%
  dplyr::left_join(SM_coverage) %>%
  dplyr::mutate(group = as.integer(group)) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(avg_lat = abs(as.numeric(latitude) - mean(as.numeric(latitude), na.rm = T)), 
                avg_lon = abs(as.numeric(longitude) - mean(as.numeric(longitude), na.rm = T))) %>%
  dplyr::mutate(unique_isotypes_per_group = length(unique(purrr::discard(isotype, is.na)))) %>%
  dplyr::group_by(isotype) %>%
  dplyr::mutate(unique_groups_per_isotype = length(unique(group))) %>%
  dplyr::group_by(strain) %>%
  dplyr::mutate(strain_in_multiple_isotypes = length(strain) > 1) %>%
  dplyr::mutate(location_issue = (avg_lat > 5 | avg_lon > 5)) %>%
  dplyr::select(-avg_lat,-avg_lon) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(strain_conflict = any(unique_isotypes_per_group > 1,
                                      unique_groups_per_isotype > 1,
                                      location_issue,
                                      na.rm = T))

ggplot(gtcheck) +
  geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.00025) +
  scale_fill_manual(values = c("#808080", "#0080FF"))

ggsave("concordance.pdf", width = 5, height = 5)
ggsave("concordance.png", width = 5, height = 5)

ggplot(gtcheck) +
  geom_histogram(aes(x=concordance, fill = isotype), binwidth = 0.000025) +
  scale_fill_manual(values = c("#808080", "#0080FF")) +
  scale_x_continuous(limits = c(0.99, 1.0)) +
  labs(x = "Concordance", y = "Number of Comparisons") +
  geom_vline(aes(xintercept = cutoff), color = "red") +
  theme(axis.title = ggplot2::element_text(size=14, face="bold", color="black", vjust=5))

ggsave("xconcordance.pdf", width = 5, height = 5)
ggsave("xconcordance.png", width = 5, height = 5)

# Save text files
readr::write_tsv(gtcheck, paste0("gtcheck.tsv"))
readr::write_tsv(isotype_groups, paste0("isotype_groups.tsv"))


# Output problem strains
pr_strains <- unique((isotype_groups %>% 
                        dplyr::filter(strain_conflict))$strain)

pr_strain_comparison <- gtcheck %>% dplyr::filter(i %in% pr_strains, j %in% pr_strains) %>% 
  dplyr::filter(discordance < 10000) %>% 
  dplyr::left_join(existing_WI %>% dplyr::rename(iso_i = isotype) %>% dplyr::select(-latitude, -longitude), by = c("i"= "strain")) %>%
  dplyr::left_join(existing_WI %>% dplyr::rename(iso_j = isotype) %>% dplyr::select(-latitude, -longitude), by = c("j"= "strain")) %>%
  dplyr::filter(!is.na(sites))

readr::write_tsv(pr_strain_comparison, paste0("problem_strains.tsv"))

isotypes <- length(unique(isotype_groups$group))

writeLines(paste0(isotypes, sep = "\n"), f)
