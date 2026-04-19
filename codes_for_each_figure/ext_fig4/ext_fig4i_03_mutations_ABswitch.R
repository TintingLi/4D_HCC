

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse)) %>% suppressWarnings()
suppressPackageStartupMessages(library(bedtoolsr)) %>% suppressWarnings()
options(bedtools.path = "/home/littt/miniforge-pypy3/bin/")
knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/results/mutants/ngs/SNV/')
setwd('~/cluster_home/project/liverCancer/results/mutants/ngs/SNV/')



snp.hc <- read_tsv('./annovar/results/ccl4_18w.snp.Somatic.hc.vcf.mm10_multianno.txt') %>% suppressMessages()
indel.hc <- read_tsv('./annovar/results/ccl4_18w.indel.Somatic.hc.vcf.mm10_multianno.txt') %>% suppressMessages()

comp.files <- list.files('~/cluster_home/project/liverCancer/results/hic/hicexplorer/male_mouse/calAB/200kb/cleanResults/', pattern = '*.bed', full.names = T)


comp <- lapply(comp.files, read.table, header =T)

names(comp) <- basename(comp.files) %>% sub('pca1_', '', x=.) %>% sub('_200000.bg.compartment.bed', '', x=.)

comp <- lapply(comp, function(x){
  return( drop_na( x[(x$seqname != 'chrY'), ] ))
})

sapply(comp, dim)
comp.switch <- list()

ccl4_18w.comp <- comp$ccl4_18w
con_18w.comp <- comp$con_18w


genome.switch.18w <- bt.intersect(comp$ccl4_18w, comp$con_18w, loj =T) %>% dplyr::select(c(1,2,3,5, 6, 12)) %>% 
  set_colnames(c('seqname', 'start', 'end', 'geneNumber', 'comp_CCL4_18w', 'comp_CTRL_18w')) 

genome.switch.18w

genome.switch.18w.SNP <- bt.intersect(genome.switch.18w, snp.hc, c =T) %>% 
  set_colnames(c('seqname', 'start', 'end', 'geneNumber', 'comp_CCL4_18w', 'comp_CTRL_18w', 'numSNP'))


genome.switch.18w.SNP$locus <- seq(1, dim(genome.switch.18w.SNP)[1])


genome.switch.18w.SNP <- as_tibble(genome.switch.18w.SNP)
genome.switch.18w.SNP

group_by(genome.switch.18w.SNP, comp_CCL4_18w,comp_CTRL_18w ) %>% summarise(n =n())


genome.switch.18w.SNP.noswitch.A <- filter(genome.switch.18w.SNP, comp_CCL4_18w == 'A' & comp_CTRL_18w == 'A')
genome.switch.18w.SNP.noswitch.B <- filter(genome.switch.18w.SNP, comp_CCL4_18w == 'B' & comp_CTRL_18w == 'B')


genome.switch.18w.SNP.BtoA <- filter(genome.switch.18w.SNP, comp_CCL4_18w == 'A' & comp_CTRL_18w == 'B')
genome.switch.18w.SNP.AtoB <- filter(genome.switch.18w.SNP, comp_CCL4_18w == 'B' & comp_CTRL_18w == 'A')


plot(genome.switch.18w.SNP.noswitch.B$locus, genome.switch.18w.SNP.noswitch.B$numSNP, pch =19,
     col = '#c2cae3', ylim = c(0, 50), xlab = 'Genome (200 Kb bin)', ylab = 'No. of mutations / 500 Kb')

points(genome.switch.18w.SNP.noswitch.A$locus, genome.switch.18w.SNP.noswitch.A$numSNP, pch =19,
     col = '#7d87b2')



points(genome.switch.18w.SNP.BtoA$locus, genome.switch.18w.SNP.BtoA$numSNP, pch =19,
       col ='#9f6e71')


points(genome.switch.18w.SNP.AtoB$locus, genome.switch.18w.SNP.AtoB$numSNP, pch =19,
       col ='#749e89')


legend("topright", legend = c('B', "A", "BtoA", "AtoB"), pch =19,
        col = c('#c2cae3','#7d87b2', '#9f6e71','#749e89'  ) ,bty = 'n')

write_csv(genome.switch.18w.SNP.AtoB, '~/cluster_home/project/liverCancer/scripts/figures_update/ext_fig4/data/ext_fig4l.csv')


getwd()


#' ### SVs
#' 
ccl4_18w_specific_SV <- read_csv( '../../sv/ccl4_18w_specific_SV.csv')
ccl4_18w_specific_SV$seqname <- gsub('Chr', 'chr', ccl4_18w_specific_SV$seqname)

genome.switch.18w.SV <- bt.intersect(genome.switch.18w, ccl4_18w_specific_SV, c =T) %>% 
  set_colnames(c('seqname', 'start', 'end', 'geneNumber', 'comp_CCL4_18w', 'comp_CTRL_18w', 'numSV'))
head(genome.switch.18w.SV)

genome.switch.18w.SV$locus <- seq(1, dim(genome.switch.18w.SV)[1])


genome.switch.18w.SV <- as_tibble(genome.switch.18w.SV)
genome.switch.18w.SV

group_by(genome.switch.18w.SV, comp_CCL4_18w,comp_CTRL_18w ) %>% summarise(n =n())


genome.switch.18w.SV.noswitch.A <- filter(genome.switch.18w.SV, comp_CCL4_18w == 'A' & comp_CTRL_18w == 'A')
genome.switch.18w.SV.noswitch.B <- filter(genome.switch.18w.SV, comp_CCL4_18w == 'B' & comp_CTRL_18w == 'B')


genome.switch.18w.SV.BtoA <- filter(genome.switch.18w.SV, comp_CCL4_18w == 'A' & comp_CTRL_18w == 'B')
genome.switch.18w.SV.AtoB <- filter(genome.switch.18w.SV, comp_CCL4_18w == 'B' & comp_CTRL_18w == 'A')


plot(genome.switch.18w.SV.noswitch.B$locus, genome.switch.18w.SV.noswitch.B$numSV, pch =19,
     col = '#c2cae3', ylim = c(0, 8), xlab = 'Genome (200 Kb bin)', ylab = 'No. of SVs / 200 Kb')

points(genome.switch.18w.SV.noswitch.A$locus, genome.switch.18w.SV.noswitch.A$numSV, pch =19,
       col = '#7d87b2')



points(genome.switch.18w.SV.BtoA$locus, genome.switch.18w.SV.BtoA$numSV, pch =19,
       col ='#9f6e71')


points(genome.switch.18w.SV.AtoB$locus, genome.switch.18w.SV.AtoB$numSV, pch =19,
       col ='#749e89')


legend("topright", legend = c('B', "A", "BtoA", "AtoB"), pch =19,
       col = c('#c2cae3','#7d87b2', '#9f6e71','#749e89'  ) ,bty = 'n')
