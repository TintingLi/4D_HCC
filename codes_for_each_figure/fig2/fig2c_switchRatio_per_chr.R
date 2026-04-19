

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(HiTC))
options(bedtools.path = "~/miniforge-pypy3/bin/")
suppressPackageStartupMessages(library(bedtoolsr))
#setwd('~/project/liverCancer/results/hic/hicexplorer/male_mouse/calAB/')

setwd('~/cluster_home/project/liverCancer/results/hic/hicexplorer/male_mouse/calAB/')

knitr::opts_knit$set(root.dir = '~/cluster_home/project/liverCancer/results/hic/hicexplorer/male_mouse/calAB/')



#' 

#' ### 1. switched A/B compartment results with 200 kb bin.
#' 
#' compartments of six groups.
comp.files <- list.files('./200kb/cleanResults/', pattern = '*.bed', full.names = T)
comp.files


comp.list <- lapply(comp.files, read.table, header =T)

names(comp.list) <- basename(comp.files) %>% sub('pca1_', '', x=.) %>% 
  sub('_200000.bg.compartment.bed', '', x=.)

comp.list <- lapply(comp.list, function(x){ return(drop_na(x)) } )

sapply(comp.list, dim)


ccl4_18w_a_ccl4_6w_b <- bt.intersect(a = comp.list$ccl4_18w[comp.list$ccl4_18w$compartment =='A',],
                                     b = comp.list$ccl4_6w[comp.list$ccl4_6w$compartment =='B', ])
dim(ccl4_18w_a_ccl4_6w_b) # 953 (200 kb window)

ccl4_6w_a_ctrl_6w_b <- bt.intersect(a = comp.list$ccl4_6w[comp.list$ccl4_6w$compartment =='A',],
                                    b = comp.list$con_6w[comp.list$con_6w$compartment =='B', ])
dim(ccl4_6w_a_ctrl_6w_b) # 771 (200 kb window)

ctrl_18w_a_ctrl_6w_b <- bt.intersect(a = comp.list$con_18w[comp.list$con_18w$compartment =='A',],
                                     b = comp.list$con_6w[comp.list$con_6w$compartment =='B', ])
dim(ctrl_18w_a_ctrl_6w_b) # 525

ccl4_18w_a_ctrl_18w_b <- bt.intersect(a = comp.list$ccl4_18w[comp.list$ccl4_18w$compartment =='A',],
                                      b = comp.list$con_18w[comp.list$con_18w$compartment =='B', ])
dim(ccl4_18w_a_ctrl_18w_b) # 1129


num.ctrl_18w_a_ctrl_6w_b <- tibble(ctrl_18w_a_ctrl_6w_b) %>% group_by(V1) %>% tally() %>% dplyr::rename( "seqname" = 'V1',  "ctrl_18w_a_ctrl_6w_b" ='n')
num.ccl4_6w_a_ctrl_6w_b <- tibble(ccl4_6w_a_ctrl_6w_b) %>% group_by(V1) %>% tally() %>% dplyr::rename("seqname" = 'V1', "ccl4_6w_a_ctrl_6w_b"='n')
num.ccl4_18w_a_ctrl_18w_b <- tibble(ccl4_18w_a_ctrl_18w_b) %>% group_by(V1) %>% tally() %>% dplyr::rename("seqname" = 'V1', "ccl4_18w_a_ctrl_18w_b"='n')
num.ccl4_18w_a_ccl4_6w_b <- tibble(ccl4_18w_a_ccl4_6w_b) %>% group_by(V1) %>% tally() %>% dplyr::rename("seqname" = 'V1', "ccl4_18w_a_ccl4_6w_b"='n')



merged <- 
left_join(num.ctrl_18w_a_ctrl_6w_b, num.ccl4_6w_a_ctrl_6w_b) %>%
  left_join(num.ccl4_18w_a_ctrl_18w_b) %>%
  left_join(num.ccl4_18w_a_ccl4_6w_b)



layout(matrix(c(1:4), nrow = 4))

par(mar = c(1,8,0,2), oma = c(6,0,6,0))



barplot(merged$ctrl_18w_a_ctrl_6w_b, ylim = c(0, 200), las =2)
title(main ='ctrl_18w_a_ctrl_6w_b', line = -3 )

barplot(merged$ccl4_6w_a_ctrl_6w_b, ylim = c(0, 200), las =2)
title(main ='ccl4_6w_a_ctrl_6w_b', line = -3 )


barplot(merged$ccl4_18w_a_ctrl_18w_b, ylim = c(0, 200), las =2)
title(main ='ccl4_18w_a_ctrl_18w_b', line = -3 )


barplot(merged$ccl4_18w_a_ccl4_6w_b, ylim = c(0, 200),  names.arg = merged$seqname, las =2)
title(main ='ccl4_18w_a_ccl4_6w_b', line = -3 )




# normalize by chr length
mm10.chrLength <- read_tsv('/public/reference/mm10/mm10.chrom.sizes', col_names = c('seqname', 'length'))


mm10.chrLength <- filter(mm10.chrLength, !grepl('random', seqname)) %>% filter(!grepl('chrUn', seqname))


merged <- left_join(merged, mm10.chrLength)

merged$ctrl_18w_a_ctrl_6w_b.norm <- merged$ctrl_18w_a_ctrl_6w_b/merged$length
merged$ccl4_6w_a_ctrl_6w_b.norrm <- merged$ccl4_6w_a_ctrl_6w_b/merged$length

merged$ccl4_18w_a_ctrl_18w_b.norm <- merged$ccl4_18w_a_ctrl_18w_b/merged$length
merged$ccl4_18w_a_ccl4_6w_b.norm <- merged$ccl4_18w_a_ccl4_6w_b/merged$length


barplot(merged$ctrl_18w_a_ctrl_6w_b.norm, las =2)
title(main ='ctrl_18w_a_ctrl_6w_b', line = -3 )


barplot(merged$ccl4_6w_a_ctrl_6w_b.norrm,  las =2, col = '#e76254')
title(main ='ccl4_6w_a_ctrl_6w_b', line = -3 )


barplot(merged$ccl4_18w_a_ctrl_18w_b, ylim = c(0, 200), las =2)
title(main ='ccl4_18w_a_ctrl_18w_b', line = -3 )



barplot(merged$ccl4_18w_a_ccl4_6w_b, ylim = c(0, 200),  names.arg = merged$seqname, las =2)
title(main ='ccl4_18w_a_ccl4_6w_b', line = -3 )


layout(matrix(1))

barplot(merged$ccl4_6w_a_ctrl_6w_b.norrm,  las = 2, col = c(rep('#376795', 19), '#e76254','#376795' ))


write_csv(merged, '~/cluster_home/project/liverCancer/scripts/figures_update/fig2/data/fig2c.csv')


