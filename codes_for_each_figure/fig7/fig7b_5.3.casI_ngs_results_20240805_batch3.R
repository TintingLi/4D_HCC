
#' ---
#' title: "5.ngs_results of casI (20240805)"
#' author: "LiTingting(ting67@126.com)"
#' output:
#'   html_document:
#'     number_sections: false
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
#suppressPackageStartupMessages(library(Biostrings))
#suppressPackageStartupMessages(library(GenomicAlignments))
#suppressPackageStartupMessages(library(Gviz))


setwd('~/cluster_home/project/liverCancer/')
knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/')

options(width = 140)



#' ## 1. data
#' ### 1.1 sgRNA info
#' selected sgRNA-gene info: 41 sgRNAs, 11 genes
casi.sgRNA.info <- read_csv('./scripts/chrX_activation/casA_I_validate_Library/data/casI_sgRNA_info.csv') %>% 
  suppressMessages()
casi.sgRNA.info
unique(casi.sgRNA.info$genename)


#
nc_sgRNAs <- c('GGATCTAGCTACCTCAAAAG', 'GAGAAGGATGGAAATTAGAA', 'CGCGCACCACGGGCGCGCAC', 
               'TTTACGATCTAGCGGCGTAG', 'TTCGTAGGAACTAAACTGTA')



casI_counts_ex1 <- read.table('./results/crispr/20240805/SG-Lib-I-EX1_R1_sgRNA_count.txt', 
                          col.names  = c('casI_count', 'target_sequence')) %>% as_tibble
casI_counts_ex1

casI_counts_ex2 <- read.table('./results/crispr/20240805/SG-Lib-I-EX2_R1_sgRNA_count.txt', 
                              col.names  = c('casI_count', 'target_sequence')) %>% as_tibble
casI_counts_ex2

casI_nc_counts <- read.table('./results/crispr/20240805/SG-Lib-I-NC_R1_sgRNA_count.txt', 
                             col.names  = c('casI_nc_count', 'target_sequence')) %>% as_tibble
casI_nc_counts


nc_sgRNAs %in% casI_counts_ex1$target_sequence # all true
nc_sgRNAs %in% casI_counts_ex2$target_sequence # all true
nc_sgRNAs %in% casI_nc_counts$target_sequence  # all true


#' **EX1**  
casI_counts_ex1 <- left_join(casI_counts_ex1, casI_nc_counts) %>% 
              filter(target_sequence %in% c(casi.sgRNA.info$target_sequence, nc_sgRNAs)) %>%
              left_join(casi.sgRNA.info) %>% print(n =25)

casI_counts_ex1 <- casI_counts_ex1[1:23, c(2, 4, 1, 3, 5)]
print(casI_counts_ex1, n =25)



casI_counts_ex1$casI_nc_frac <- casI_counts_ex1$casI_nc_count / sum(casI_counts_ex1$casI_nc_count)

casI_counts_ex1$casI_frac <- casI_counts_ex1$casI_count / sum(casI_counts_ex1$casI_count)


casI_counts_ex1$up_ratio <- casI_counts_ex1$casI_frac / casI_counts_ex1$casI_nc_frac


casI_counts_ex1 %<>% arrange(desc(up_ratio))


casI_counts_ex1 %>% print(n =30)




#' **EX2**  
casI_counts_ex2 <- left_join(casI_counts_ex2, casI_nc_counts) %>% 
  filter(target_sequence %in% c(casi.sgRNA.info$target_sequence, nc_sgRNAs)) %>%
  left_join(casi.sgRNA.info) %>% print(n =25)

casI_counts_ex2 <- casI_counts_ex2[1:23, c(2, 4, 1, 3, 5)]
print(casI_counts_ex2, n =25)



casI_counts_ex2$casI_nc_frac <- casI_counts_ex2$casI_nc_count / sum(casI_counts_ex2$casI_nc_count)

casI_counts_ex2$casI_frac <- casI_counts_ex2$casI_count / sum(casI_counts_ex2$casI_count)


casI_counts_ex2$up_ratio <- casI_counts_ex2$casI_frac / casI_counts_ex2$casI_nc_frac


casI_counts_ex2 %<>% arrange(desc(up_ratio))


casI_counts_ex2 %>% print(n =30)







casI_counts_ex2$genename[c(1, 7, 10)] <- paste0('nc_sgRNA', 1:3)
casI_counts_ex2$seqname[c(1, 7, 10)] <- 'nc'




#+ fig.width = 6, fig.height = 5
layout(matrix(1))
par(mar = c(4,4,2,1), oma = c(4,4,2,1))
plot(casI_counts_ex2$casI_count[casI_counts_ex2$seqname=='chrX'],
     casI_counts_ex2$up_ratio[casI_counts_ex2$seqname=='chrX'] %>% log2(), 
     xlim = c(12e3, 75e4), ylim = c(-2, 2),
     type ='p', col ='grey', pch = 16, cex =2,
     xlab = 'read count', ylab = 'Fold change of sgRNA abundance (log2)')



points(casI_counts_ex2$casI_count[casI_counts_ex2$seqname=='nc'],
       casI_counts_ex2$up_ratio[casI_counts_ex2$seqname=='nc'] %>% log2(), 
       col ='black', pch = 16, cex =2)

abline(h = 1, col = 'grey', lty = 2 )


write_csv(casI_counts_ex2, '~/cluster_home/project/liverCancer/scripts/figures_update/fig7/data/fig7b.csv')



























