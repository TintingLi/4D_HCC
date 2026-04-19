
#' ---
#' title: "5.ngs_results of casA (20240805)"
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
#' selected sgRNA-gene info: 42 sgRNAs, 14 genes
casA.sgRNA.info <- read_csv('./scripts/chrX_activation/casA_I_validate_Library/data/casA_sgRNA_info.csv') %>% 
  suppressMessages()
casA.sgRNA.info
unique(casA.sgRNA.info$genename)


#
nc_sgRNAs <- c('GGATCTAGCTACCTCAAAAG', 'GAGAAGGATGGAAATTAGAA', 'CGCGCACCACGGGCGCGCAC', 
               'TTTACGATCTAGCGGCGTAG', 'TTCGTAGGAACTAAACTGTA')



casA_counts_ex1 <- read.table('./results/crispr/20240805/SG-Lib-A-EX1_R1_sgRNA_count.txt', 
                          col.names  = c('casA_count', 'target_sequence')) %>% as_tibble
casA_counts_ex1

casA_counts_ex2 <- read.table('./results/crispr/20240805/SG-Lib-A-EX2_R1_sgRNA_count.txt', 
                              col.names  = c('casA_count', 'target_sequence')) %>% as_tibble
casA_counts_ex2


casA_counts_ex3 <- read.table('./results/crispr/20240805/SG-Lib-A-EX3_R1_sgRNA_count.txt', 
                              col.names  = c('casA_count', 'target_sequence')) %>% as_tibble
casA_counts_ex3

casA_nc_counts <- read.table('./results/crispr/20240805/SG-Lib-A-NC_R1_sgRNA_count.txt', 
                             col.names  = c('casA_nc_count', 'target_sequence')) %>% as_tibble
casA_nc_counts


nc_sgRNAs %in% casA_counts_ex1$target_sequence # all true
nc_sgRNAs %in% casA_counts_ex2$target_sequence # all true
nc_sgRNAs %in% casA_counts_ex3$target_sequence # all true
nc_sgRNAs %in% casA_nc_counts$target_sequence  # all true





#' **EX1**  
casA_counts_ex1 <- left_join(casA_counts_ex1, casA_nc_counts) %>% 
              filter(target_sequence %in% c(casA.sgRNA.info$target_sequence, nc_sgRNAs)) %>%
              left_join(casA.sgRNA.info) %>% print(n =33)

casA_counts_ex1 <- casA_counts_ex1[1:29, c(2, 4, 1, 3, 5)]




casA_counts_ex1$casA_nc_frac <- casA_counts_ex1$casA_nc_count / sum(casA_counts_ex1$casA_nc_count)

casA_counts_ex1$casA_frac <- casA_counts_ex1$casA_count / sum(casA_counts_ex1$casA_count)


casA_counts_ex1$up_ratio <- casA_counts_ex1$casA_frac / casA_counts_ex1$casA_nc_frac


casA_counts_ex1 %<>% arrange(desc(up_ratio))


casA_counts_ex1 %>% print(n =30)



#' **EX2**  
casA_counts_ex2 <- left_join(casA_counts_ex2, casA_nc_counts) %>% 
  filter(target_sequence %in% c(casA.sgRNA.info$target_sequence, nc_sgRNAs)) %>%
  left_join(casA.sgRNA.info) %>% print(n =33)

casA_counts_ex2 <- casA_counts_ex2[1:29, c(2, 4, 1, 3, 5)]




casA_counts_ex2$casA_nc_frac <- casA_counts_ex2$casA_nc_count / sum(casA_counts_ex2$casA_nc_count)

casA_counts_ex2$casA_frac <- casA_counts_ex2$casA_count / sum(casA_counts_ex2$casA_count)


casA_counts_ex2$up_ratio <- casA_counts_ex2$casA_frac / casA_counts_ex2$casA_nc_frac


casA_counts_ex2 %<>% arrange(desc(up_ratio))


casA_counts_ex2 %>% print(n =30)






#' **EX3**  
casA_counts_ex3 <- left_join(casA_counts_ex3, casA_nc_counts) %>% 
  filter(target_sequence %in% c(casA.sgRNA.info$target_sequence, nc_sgRNAs)) %>%
  left_join(casA.sgRNA.info) %>% print(n =33)

casA_counts_ex3 <- casA_counts_ex3[1:29, c(2, 4, 1, 3, 5)]




casA_counts_ex3$casA_nc_frac <- casA_counts_ex3$casA_nc_count / sum(casA_counts_ex3$casA_nc_count)

casA_counts_ex3$casA_frac <- casA_counts_ex3$casA_count / sum(casA_counts_ex3$casA_count)


casA_counts_ex3$up_ratio <- casA_counts_ex3$casA_frac / casA_counts_ex3$casA_nc_frac


casA_counts_ex3 %<>% arrange(desc(up_ratio))


casA_counts_ex3 %>% print(n =30)





casA_counts_ex3$genename[c(9, 13, 26:28)] <- paste0('nc_sgRNA', 1:5)
casA_counts_ex3$seqname[c(9, 13, 26:28)] <- 'nc'



#'
#casA_counts$genename[25:29] <- paste0('nc_sgRNA', 1:5)
#casA_counts$seqname[25:29] <- 'nc'
#
#
#+ fig.width = 6, fig.height = 5
layout(matrix(1))
par(mar = c(4,4,2,1), oma = c(4,4,2,1))
plot(casA_counts_ex3$casA_count[casA_counts_ex3$seqname=='chrX'],
     casA_counts_ex3$up_ratio[casA_counts_ex3$seqname=='chrX'] %>% log2(), 
     xlim = c(12e3, 65e4), ylim = c(-2, 3),
     type ='p', col ='grey', pch = 16, cex =2,
     xlab = 'read count', ylab = 'Fold change of sgRNA abundance (log2)')


points(casA_counts_ex3$casA_count[casA_counts_ex3$up_ratio > 2],
       casA_counts_ex3$up_ratio[casA_counts_ex3$up_ratio > 2] %>% log2(),  
       col ='#A72126', pch = 16, cex =2)



points(casA_counts_ex3$casA_count[casA_counts_ex3$seqname=='nc'],
       casA_counts_ex3$up_ratio[casA_counts_ex3$seqname=='nc'] %>% log2(), 
      col ='black', pch = 16, cex =2)

abline(h = 1, col = 'grey', lty = 2 )



write_csv(casA_counts_ex3, '~/cluster_home/project/liverCancer/scripts/figures_update/fig7/data/fig7a.csv')










