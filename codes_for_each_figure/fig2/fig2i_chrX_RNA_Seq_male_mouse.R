
#' ---
#' title: "1.1 chrX activation in liver Cancer mouse model (male)"
#' author: "LiTingting(ting67@126.com)"
#' output:
#'   html_document:
#'     number_sections: false
#'     highlight: pygments
#'     theme: cerulean
#'     toc: yes
#' ---

options(width = 120)


suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HiTC))
suppressPackageStartupMessages(library(GENOVA))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(txdbmaker))


setwd('~/cluster_home/project/liverCancer/results/rnaseq/star_deseq2/')
knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/results/rnaseq/star_deseq2/')

options(width =120)

RNAseq.ccl4_18w_vs_ctrl_18w <- read_csv('./male_mouse/output/deseq2/ccl4_18w_vs_ctrl_18w.tpm.csv') %>% suppressMessages()
dim(RNAseq.ccl4_18w_vs_ctrl_18w)


deg.ccl4_18w_vs_ctrl_18w <- filter(RNAseq.ccl4_18w_vs_ctrl_18w,
								   abs(log2FC_ccl4_18w_vs_ctrl_18w) > 1 & padj_ccl4_18w_vs_ctrl_18w < 0.05)
dim(deg.ccl4_18w_vs_ctrl_18w)
names(deg.ccl4_18w_vs_ctrl_18w)
deg.ccl4_18w_vs_ctrl_18w[1:3, 1:7]



gtffile <- makeTxDbFromGFF('~/cluster_home/ref/GENCODE/m38_vM25/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf') %>% suppressMessages()
mm10.gtf <- genes(gtffile) %>% as_tibble()
mm10.gtf

table(deg.ccl4_18w_vs_ctrl_18w$ensembl_id %in% mm10.gtf$gene_id)

deg.ccl4_18w_vs_ctrl_18w.chrInfo <- left_join(deg.ccl4_18w_vs_ctrl_18w[, 1:5 ], mm10.gtf, by = c('ensembl_id' = 'gene_id'))
deg.ccl4_18w_vs_ctrl_18w.chrInfo

#' total DEG numbers per chromosome.   
table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames)[1:22]

#' ### UP & DOWN DEG numbers per chromosome.  
#' 
#' up DEGs    
#' 
table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w > 1])[1:22]
#' down DEGs 
table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w < -1])[1:22]

#' up/down DEGs
table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w > 1])[1:21] /
table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w < -1])[1:21]

#' plot of ratio of up/down DEGs per chromosome  
#' **chrX showed the highest up/down DEG gene number ratio across the whole genome**  

(table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w > 1])[1:21] /
		table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w < -1])[1:21]) %>% 
	barplot(las =2, ylab = 'up/down DEG gene number', ylim = c(0, 5), col = c(rep('#528fad', 19), '#f7aa58','#528fad' ))


tmp <- (table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w > 1])[1:21] /
    table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w < -1])[1:21]) %>% as_vector()

data_upload <- tibble(seqname = names(tmp), ratio = tmp %>% as.numeric() )

write_csv(data_upload, '~/cluster_home/project/liverCancer/scripts/figures_update/fig2/data/fig2i.csv')




layout(matrix(1))
par(mar = c(4,5,2,1), oma = c(3,4,2,1))

#' DEG ratios compared to total gene numbers per chromosome  
#' 
#+ fig.width =8, fig.height =5  

(table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames)[1:21] / table(mm10.gtf$seqnames)[1:21]) %>% 
  barplot(las =2, ylab = 'total_DEGs / total_gene_numbers',  col = c(rep('#528fad', 19), '#f7aa58','#528fad' ))


(table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w > 1])[1:21] / 
    table(mm10.gtf$seqnames)[1:21]) %>% 
  barplot(las =2, ylab = 'up_DEGs / total_gene_numbers',  col = c(rep('#528fad', 19), '#f7aa58','#528fad' ))


(table(deg.ccl4_18w_vs_ctrl_18w.chrInfo$seqnames[deg.ccl4_18w_vs_ctrl_18w.chrInfo$log2FC_ccl4_18w_vs_ctrl_18w < -1])[1:21] / 
    table(mm10.gtf$seqnames)[1:21]) %>% 
  barplot(las =2, ylab = 'down_DEGs / total_gene_numbers',  col = c(rep('#528fad', 19), '#f7aa58','#528fad' ))




