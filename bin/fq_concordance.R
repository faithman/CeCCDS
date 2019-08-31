library(tidyverse)

gt_dict = list("0/0" = 0, "1/1" = 1)

df <- readr::read_tsv("rg_gt.tsv", col_names = c("CHROM", "POS", "gt", "SM", "fq")) %>%
      tidyr::unite("CHROM_POS", CHROM, POS, sep = "_") %>%
      dplyr::filter(gt %in% c("0/0", "1/1")) %>%
      { 
        if (nrow(.) > 0)
            dplyr::rowwise(.) %>%
            dplyr::mutate(gt = gt_dict[gt][[1]]) %>% 
            dplyr::ungroup()
        else 
            .
      }

# No concordant sites
if (nrow(df) == 0) {
    q()
}

SM <- df$SM[[1]]

df <- df %>%
       tidyr::spread(fq, gt) %>%
       dplyr::select(-CHROM_POS, -SM)

# Single run isotype
if(ncol(df) == 1) {
  q()
} else {

    # Calculate the concordance between read groups (individual FASTQ pairs)

    # Alt only
    alt_only <- apply(df, 2, function(x) colSums(x==df, na.rm = T) )
    total <- apply(df, 2, function(x) colSums(!is.na(df) & !is.na(x)) )

    alt_only <- alt_only %>% dplyr::tbl_df() %>% 
      dplyr::mutate(b = colnames(alt_only)) %>%
      dplyr::select(b, dplyr::everything()) %>%
      tidyr::gather(fq, concordant_sites, -b) %>% 
      dplyr::rename(a = fq)

    total <- total %>% dplyr::tbl_df() %>% 
      dplyr::mutate(b = colnames(total)) %>%
      dplyr::select(b, dplyr::everything()) %>%
      tidyr::gather(fq, total_sites, -b) %>% 
      dplyr::rename(a = fq)
      
    df <- dplyr::left_join(alt_only, total) %>%
          dplyr::select(a, b, dplyr::everything()) %>%
          dplyr::mutate(concordance = concordant_sites/total_sites, SM = SM) %>%
          dplyr::filter(a != b) %>%
          dplyr::mutate(concordant_sites = as.integer(concordant_sites),
                        total_sites = as.integer(total_sites))

    readr::write_tsv(df, "out.tsv", col_names = F)

}