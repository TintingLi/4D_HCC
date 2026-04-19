suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HiTC))
library(bedtoolsr)
knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/results/')

setwd('~/cluster_home/project/liverCancer/results/')
options(bedtools.path = "~/miniforge-pypy3/bin/")

comp.files <- list.files('./hic/hicexplorer/male_mouse/calAB/200kb/cleanResults/', pattern = '*.bed', full.names = T)

comp <- lapply(comp.files, read.table, header =T)

names(comp) <- basename(comp.files) %>% sub('pca1_', '', x=.) %>% sub('_200000.bg.compartment.bed', '', x=.)

comp <- lapply(comp, function(x){
  return( drop_na( x[(x$seqname != 'chrY'), ] ))
})

sapply(comp, dim)


comp.switch <- list()

ccl4_6w.comp <- comp$ccl4_6w
ccl4_10w.comp <- comp$ccl4_10w
ccl4_18w.comp <- comp$ccl4_18w
con_6w.comp <- comp$con_6w
con_10w.comp <- comp$con_10w
con_18w.comp <- comp$con_18w


# 9 groups with a->b & b-> a switch:
#
# ctrl_10w_vs_ctrl_6w
# ctrl_18w_vs_ctrl_10w
# ctrl_18w_vs_ctrl_6w
# 
comp.switch[["ctrl_10w_vs_ctrl_6w_a2b"]] <- bt.intersect(a = con_10w.comp[con_10w.comp$compartment=='A',], 
                                                         b = con_6w.comp[con_6w.comp$compartment == 'B',] )
comp.switch[["ctrl_10w_vs_ctrl_6w_b2a"]] <- bt.intersect(a = con_10w.comp[con_10w.comp$compartment=='B',], 
                                                         b = con_6w.comp[con_6w.comp$compartment == 'A',] )

comp.switch[["ctrl_18w_vs_ctrl_10w_a2b"]] <- bt.intersect(a = con_18w.comp[con_18w.comp$compartment=='A',], 
                                                          b = con_10w.comp[con_10w.comp$compartment == 'B',] )
comp.switch[["ctrl_18w_vs_ctrl_10w_b2a"]] <- bt.intersect(a = con_18w.comp[con_18w.comp$compartment=='B',], 
                                                          b = con_10w.comp[con_10w.comp$compartment == 'A',] )

comp.switch[["ctrl_18w_vs_ctrl_6w_a2b"]] <- bt.intersect(a = con_18w.comp[con_18w.comp$compartment=='A',], 
                                                         b = con_6w.comp[con_6w.comp$compartment == 'B',] )
comp.switch[["ctrl_18w_vs_ctrl_6w_b2a"]] <- bt.intersect(a = con_18w.comp[con_18w.comp$compartment=='B',], 
                                                         b = con_6w.comp[con_6w.comp$compartment == 'A',] )
# ccl4_6w_vs_ctrl_6w
# ccl4_10w_vs_ctrl_10w
# ccl4_18w_vs_ctrl_18w
# 
comp.switch[["ccl4_6w_vs_ctrl_6w_a2b"]] <- bt.intersect(a = ccl4_6w.comp[ccl4_6w.comp$compartment=='A',], 
                                                        b = con_6w.comp[con_6w.comp$compartment == 'B',] )
comp.switch[["ccl4_6w_vs_ctrl_6w_b2a"]] <- bt.intersect(a = ccl4_6w.comp[ccl4_6w.comp$compartment=='B',], 
                                                        b = con_6w.comp[con_6w.comp$compartment == 'A',] )

comp.switch[["ccl4_10w_vs_ctrl_10w_a2b"]] <- bt.intersect(a = ccl4_10w.comp[ccl4_10w.comp$compartment=='A',], 
                                                          b = con_10w.comp[con_10w.comp$compartment == 'B',] )
comp.switch[["ccl4_10w_vs_ctrl_10w_b2a"]] <- bt.intersect(a = ccl4_10w.comp[ccl4_10w.comp$compartment=='B',], 
                                                          b = con_10w.comp[con_10w.comp$compartment == 'A',] )

comp.switch[["ccl4_18w_vs_ctrl_18w_a2b"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='A',], 
                                                          b = con_18w.comp[con_18w.comp$compartment == 'B',] )
comp.switch[["ccl4_18w_vs_ctrl_18w_b2a"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='B',], 
                                                          b = con_18w.comp[con_18w.comp$compartment == 'A',] )

# ccl4_10w_vs_ccl4_6w
# ccl4_18w_vs_ccl4_6w
# ccl4_18w_vs_ccl4_10w
#
comp.switch[["ccl4_10w_vs_ccl4_6w_a2b"]] <- bt.intersect(a = ccl4_10w.comp[ccl4_10w.comp$compartment=='A',], 
                                                         b = ccl4_6w.comp[ccl4_6w.comp$compartment == 'B',] )
comp.switch[["ccl4_10w_vs_ccl4_6w_b2a"]] <- bt.intersect(a = ccl4_10w.comp[ccl4_10w.comp$compartment=='B',], 
                                                         b = ccl4_6w.comp[ccl4_6w.comp$compartment == 'A',] )

comp.switch[["ccl4_18w_vs_ccl4_6w_a2b"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='A',], 
                                                         b = ccl4_6w.comp[ccl4_6w.comp$compartment == 'B',] )
comp.switch[["ccl4_18w_vs_ccl4_6w_b2a"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='B',], 
                                                         b = ccl4_6w.comp[ccl4_6w.comp$compartment == 'A',] )

comp.switch[["ccl4_18w_vs_ccl4_10w_a2b"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='A',], 
                                                          b = ccl4_10w.comp[ccl4_10w.comp$compartment == 'B',] )
comp.switch[["ccl4_18w_vs_ccl4_10w_b2a"]] <- bt.intersect(a = ccl4_18w.comp[ccl4_18w.comp$compartment=='B',], 
                                                          b = ccl4_10w.comp[ccl4_10w.comp$compartment == 'A',] )


sapply(comp.switch, dim)







#RNA-seq results:

rnaseqfiles <- list.files('./rnaseq/star_deseq2/male_mouse/output/deseq2/',
						  pattern = '*tpm.csv', recursive = T, full.names = T)
rnaseqfiles


rnaseq.results <- lapply(rnaseqfiles, read_csv) %>% suppressMessages()
rnaseq.results[[1]]

rnaseq.merged <- rnaseq.results[[1]][, c(1,2)]
for (i in 2:9) {
  rnaseq.merged <- left_join(rnaseq.merged, rnaseq.results[[i]][, c(1,2)])
}
rnaseq.merged


#genecode:
#mm10.gen <- read_tsv('~/Desktop/code/database/GENCODE/fetched_info/mm10_genecode/genes_gencode_clean_mm10.txt')

mm10.gen <- read_tsv('./metadata/mm10_genecode/genes_gencode_clean_mm10.txt')



table(rnaseq.merged$ensembl_id %in% mm10.gen$gene_id) # T 29704

rnaseq.merged <- left_join(rnaseq.merged, mm10.gen, by = c('ensembl_id' = 'gene_id' ))



#
groupnames <- basename(rnaseqfiles) %>% gsub('.tpm.csv', '', x = .)



logFc_geneloci <- list()
for (i in groupnames) {
  logFc_geneloci[[i]] <- rnaseq.merged[, c('seqname','start', 'end', 'ensembl_id', paste0('log2FC_', i))]
}

logFc_geneloci


t1 <- bt.intersect(a = logFc_geneloci$ccl4_10w_vs_ccl4_6w, b = comp.switch$ccl4_10w_vs_ccl4_6w_a2b) 
t2 <- bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_6w, b = comp.switch$ccl4_18w_vs_ccl4_6w_a2b) 
t3 <- bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_10w, b = comp.switch$ccl4_18w_vs_ccl4_10w_a2b) 


median(t1$V5)
median(t2$V5)


boxplot(list(t1$V5, t2$V5, t3$V5))


compSwitch_genelog2FC <- list(
  ctrl_10w_vs_ctrl_6w_a2b = bt.intersect(a = logFc_geneloci$ctrl_10w_vs_ctrl_6w, b = comp.switch$ctrl_10w_vs_ctrl_6w_a2b)$V5,
  ctrl_10w_vs_ctrl_6w_b2a = bt.intersect(a = logFc_geneloci$ctrl_10w_vs_ctrl_6w, b = comp.switch$ctrl_10w_vs_ctrl_6w_b2a)$V5,
  
  ctrl_18w_vs_ctrl_6w_a2b = bt.intersect(a = logFc_geneloci$ctrl_18w_vs_ctrl_6w, b = comp.switch$ctrl_18w_vs_ctrl_6w_a2b)$V5,
  ctrl_18w_vs_ctrl_6w_b2a = bt.intersect(a = logFc_geneloci$ctrl_18w_vs_ctrl_6w, b = comp.switch$ctrl_18w_vs_ctrl_6w_b2a)$V5,
  
  ctrl_18w_vs_ctrl_10w_a2b = bt.intersect(a = logFc_geneloci$ctrl_18w_vs_ctrl_10w, b = comp.switch$ctrl_18w_vs_ctrl_10w_a2b)$V5,
  ctrl_18w_vs_ctrl_10w_b2a = bt.intersect(a = logFc_geneloci$ctrl_18w_vs_ctrl_10w, b = comp.switch$ctrl_18w_vs_ctrl_10w_b2a)$V5,
  
  
  ccl4_10w_vs_ccl4_6w_a2b = bt.intersect(a = logFc_geneloci$ccl4_10w_vs_ccl4_6w, b = comp.switch$ccl4_10w_vs_ccl4_6w_a2b)$V5,
  ccl4_10w_vs_ccl4_6w_b2a = bt.intersect(a = logFc_geneloci$ccl4_10w_vs_ccl4_6w, b = comp.switch$ccl4_10w_vs_ccl4_6w_b2a)$V5,
  
  ccl4_18w_vs_ccl4_6w_a2b = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_6w, b = comp.switch$ccl4_18w_vs_ccl4_6w_a2b)$V5,
  ccl4_18w_vs_ccl4_6w_b2a = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_6w, b = comp.switch$ccl4_18w_vs_ccl4_6w_b2a)$V5,
  
  ccl4_18w_vs_ccl4_10w_a2b = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_10w, b = comp.switch$ccl4_18w_vs_ccl4_10w_a2b)$V5,
  ccl4_18w_vs_ccl4_10w_b2a = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ccl4_10w, b = comp.switch$ccl4_18w_vs_ccl4_10w_b2a)$V5,
  
  
  
  
  ccl4_6w_vs_ctrl_6w_a2b = bt.intersect(a = logFc_geneloci$ccl4_6w_vs_ctrl_6w, b = comp.switch$ccl4_6w_vs_ctrl_6w_a2b)$V5,
  ccl4_6w_vs_ctrl_6w_b2a = bt.intersect(a = logFc_geneloci$ccl4_6w_vs_ctrl_6w, b = comp.switch$ccl4_6w_vs_ctrl_6w_b2a)$V5,
  
  ccl4_10w_vs_ctrl_10w_a2b = bt.intersect(a = logFc_geneloci$ccl4_10w_vs_ctrl_10w, b = comp.switch$ccl4_10w_vs_ctrl_10w_a2b)$V5,
  ccl4_10w_vs_ctrl_10w_b2a = bt.intersect(a = logFc_geneloci$ccl4_10w_vs_ctrl_10w, b = comp.switch$ccl4_10w_vs_ctrl_10w_b2a)$V5,
  
  ccl4_18w_vs_ctrl_18w_a2b = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ctrl_18w, b = comp.switch$ccl4_18w_vs_ctrl_18w_a2b)$V5,
  ccl4_18w_vs_ctrl_18w_b2a = bt.intersect(a = logFc_geneloci$ccl4_18w_vs_ctrl_18w, b = comp.switch$ccl4_18w_vs_ctrl_18w_b2a)$V5
  
  
)


par(oma = c(8,2,2,1), oma =c(8,2,2,1))
boxplot(compSwitch_genelog2FC, las =2, cex =0.7, col = rep(c('#C25160', '#2A6E3F'), times = 9))



# ctrl_10w_vs_ctrl_6w
# ctrl_18w_vs_ctrl_10w
# ctrl_18w_vs_ctrl_6w

# ccl4_10w_vs_ccl4_6w
# ccl4_18w_vs_ccl4_6w
# ccl4_18w_vs_ccl4_10w

# ccl4_6w_vs_ctrl_6w
# ccl4_10w_vs_ctrl_10w
# ccl4_18w_vs_ctrl_18w
