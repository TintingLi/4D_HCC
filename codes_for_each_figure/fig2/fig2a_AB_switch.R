
#' **calculate the A/B switches between 6 groups and visualization**
#' 

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(bedtoolsr))
options(bedtools.path = "/opt/bedtools2/bin/")



knitr::opts_knit$set(root.dir = '~/ltt/projects/liverCancer/results/hicexplorer/calAB/200kb/')
setwd('~/ltt/projects/liverCancer/results/hicexplorer/calAB/200kb/')



#' ### 1. clean A/B compartment results with 200 kb bin.
#' 

comp.files <- list.files('./cleanResults', pattern = '*.bed', full.names = T)
comp.files


comp.list <- lapply(comp.files, read.table, header =T)

names(comp.list) <- basename(comp.files) %>% sub('pca1_', '', x=.) %>% 
  sub('_200000.bg.compartment.bed', '', x=.)

comp <- lapply(comp.list, function(x){
  return( drop_na( x[(x$seqname != 'chrY'), ] ))
})
names(comp) <- names(comp.list)

sapply(comp, dim)


#' **6w example**. 
ctrl_6w.a <- comp$con_6w[comp$con_6w$compartment =='A', ]
ctrl_6w.b <- comp$con_6w[comp$con_6w$compartment =='B',]

ccl4_6w.a <- comp$ccl4_6w[comp$ccl4_6w$compartment =='A',]
ccl4_6w.b <- comp$ccl4_6w[comp$ccl4_6w$compartment =='B',]

bt.intersect(a = ccl4_6w.b, b = ctrl_6w.a) %>% dim    # a->b 258
bt.intersect(a = ccl4_6w.a, b = ctrl_6w.a) %>% dim    # a->a 4952
bt.intersect(a = ccl4_6w.b, b = ctrl_6w.b) %>% dim    # b->b 6875
bt.intersect(a = ccl4_6w.a, b = ctrl_6w.b) %>% dim    # b->a 740


#' ### 2. get the switch matrix among the six groups
#' 
switch.mat <- matrix(data = NA, nrow = 6, ncol = 6)
rownames(switch.mat) <- c('con_6w', 'con_10w', 'con_18w', 'ccl4_6w', 'ccl4_10w', 'ccl4_18w')
colnames(switch.mat) <- rownames(switch.mat)
switch.mat

#rows: a --> b
#cols: b --> a

for (i in rownames(switch.mat)) {
  for (j in colnames(switch.mat)) {
   # if(i ==j){ next}
    #a->b
    row.a <- comp[[i]][comp[[i]]$compartment =='A',]
    row.b <- comp[[j]][comp[[j]]$compartment =='B',]
    row.aTob <-  dim(bt.intersect(a = row.a, b = row.b))[1]
    switch.mat[i,j] <- row.aTob
    
    #b->a
#    col.b <- comp[[i]][comp[[i]]$compartment =='B',]
#    col.a <- comp[[j]][comp[[j]]$compartment =='A',]
#    col.bToa <-  dim(bt.intersect(a = col.b, b = col.a))[1]
#    switch.mat[j,i] <- col.bToa
  
  }
  
}



rownames(switch.mat) <- paste0( c('con_6w', 'con_10w', 'con_18w', 'ccl4_6w', 'ccl4_10w', 'ccl4_18w'), '_A_comp')
colnames(switch.mat) <- paste0( c('con_6w', 'con_10w', 'con_18w', 'ccl4_6w', 'ccl4_10w', 'ccl4_18w'), '_B_comp')



write.csv(switch.mat, './cleanResults/AB_switch.200kb.matrix.csv')


#' **visualization**
#+ fig.width=6, fig.height=6, fig.align = 'center'
switch.mat <- read.csv('/home/littt/project/liverCancer/results/hic/hicexplorer/calAB/200kb/cleanResults/AB_switch.200kb.matrix.csv',row.names = 1)
pheatmap(switch.mat, cluster_rows = F, cluster_cols = F,
         display_numbers = T, number_format= '%.0f')

pheatmap(switch.mat, cluster_rows = F, cluster_cols = F,
		 display_numbers = T, number_format= '%.0f', cex =2, show_rownames =F,show_colnames=F )


