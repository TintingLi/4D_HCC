
#' ---
#' title: "3.visualization_of_sex_dependent_DEGs_of_LIHC"
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
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library(FirebrowseR))
suppressPackageStartupMessages(library(stringr))

setwd('~/cluster_home/project/liverCancer/scripts/chrX_activation/tcga/')
knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/scripts/chrX_activation/tcga/')

options(width = 150)


#-----------------------------------------------------#
# 1. stats of TCGA rna-seq data, including sample size, gender of 33 cancer types.  
# 2. 


#-----------------------------------------------------#


stat.cancer_types.female <- read_rds('./data/stat.cancer_types.female.rds')
stat.cancer_types.male <- read_rds('./data/stat.cancer_types.male.rds')



#' ##  1. Female LIHC
#' ### 1.1 UP_DEG ratios per chromosome in female LIHC 
gene_per_chr <- stat.cancer_types.female$LIHC %>% group_by(seqname) %>% summarise(gene_num= n())
gene_per_chr
up_degs_in_LIHC_per_chr <- filter(stat.cancer_types.female$LIHC,  foldchage >=2) %>% 
                    group_by(seqname) %>% summarise(up_deg_num= n())

lihc_stat <- left_join(up_degs_in_LIHC_per_chr, gene_per_chr)
lihc_stat$ratio <- lihc_stat$up_deg_num / lihc_stat$gene_num
lihc_stat <- arrange(lihc_stat, desc(ratio))
lihc_stat
#+ fig.width=8, fig.height=4
barplot(lihc_stat$ratio, ylim = c(0, 0.08), names.arg = lihc_stat$seqname, las =2,
        col = c('#E5736E', rep('#3E5885', 22)),
        main = 'UP_DEGs / total_gene_num per chromosome in female LIHC')

#' ### 1.2 Down_DEG ratios per chromosome in female LIHC 
gene_per_chr <- stat.cancer_types.female$LIHC %>% group_by(seqname) %>% summarise(gene_num= n())
down_degs_in_LIHC_per_chr <- filter(stat.cancer_types.female$LIHC,  foldchage <= 0.5) %>% 
  group_by(seqname) %>% summarise(up_deg_num= n())

lihc_stat <- left_join(down_degs_in_LIHC_per_chr, gene_per_chr)
lihc_stat$ratio <- lihc_stat$up_deg_num / lihc_stat$gene_num

lihc_stat <- arrange(lihc_stat, desc(ratio))

barplot(lihc_stat$ratio, ylim = c(0, 0.08), names.arg = lihc_stat$seqname, las =2,
        main = 'DOWN_DEGs / total_gene_num per chromosome in female LIHC')

female_data <- lihc_stat
colnames(female_data)[4] <- 'female_ratio'


#' ##  2. Male LIHC
#' ### 1.1 UP_DEG ratios per chromosome in male LIHC 
gene_per_chr <- stat.cancer_types.male$LIHC %>% group_by(seqname) %>% summarise(gene_num= n())
up_degs_in_LIHC_per_chr <- filter(stat.cancer_types.male$LIHC,  foldchage >=2) %>% 
  group_by(seqname) %>% summarise(up_deg_num= n())

lihc_stat <- left_join(up_degs_in_LIHC_per_chr, gene_per_chr)
lihc_stat$ratio <- lihc_stat$up_deg_num / lihc_stat$gene_num

lihc_stat <- arrange(lihc_stat, desc(ratio))

#+ fig.width=8, fig.height=4
barplot(lihc_stat$ratio, ylim = c(0, 0.06), names.arg = lihc_stat$seqname, las =2,
        col = c(rep('#3E5885', 2), '#E5736E', rep('#3E5885', 20)),
        main = 'UP_DEGs / total_gene_num per chromosome in male LIHC')

male_data <- lihc_stat

colnames(male_data)[4] <- 'male_ratio'


merged_data <- left_join(female_data[, c(1,4)], male_data[, c(1, 4)])

write_csv(merged_data, '~/cluster_home/project/liverCancer/scripts/figures_update/ext_fig8/data/ext_fig8e.csv')



#' ### 1.2 Down_DEG ratios per chromosome in male LIHC 
gene_per_chr <- stat.cancer_types.male$LIHC %>% group_by(seqname) %>% summarise(gene_num= n())
down_degs_in_LIHC_per_chr <- filter(stat.cancer_types.male$LIHC,  foldchage <= 0.5) %>% 
  group_by(seqname) %>% summarise(up_deg_num= n())

lihc_stat <- left_join(down_degs_in_LIHC_per_chr, gene_per_chr)
lihc_stat$ratio <- lihc_stat$up_deg_num / lihc_stat$gene_num

lihc_stat <- arrange(lihc_stat, desc(ratio))

barplot(lihc_stat$ratio, ylim = c(0, 0.08), names.arg = lihc_stat$seqname, las =2,
        main = 'DOWN_DEGs / total_gene_num per chromosome in male LIHC')











