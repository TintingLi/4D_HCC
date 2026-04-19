

suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tidyverse))
#suppressPackageStartupMessages(library(HiTC))
#suppressPackageStartupMessages(library(GENOVA))
#suppressPackageStartupMessages(library(GenomicFeatures))


setwd('~/cluster_home/project/liverCancer/results/specific/AA_BB_interaction_intensity/00.bin_bedpe/male_mouse/')
knitr::opts_knit$set(root.dir = '~/cluster_home/project/project/liverCancer/results/specific/AA_BB_interaction_intensity/00.bin_bedpe/male_mouse/')


#' ### bedpe of 200kb bins   
#' 

bedpe.files <- list.files('./compAB_annotated/male_mouse/', pattern = "*.bedpe", full.names = T)
bedpe.files
compAB.bedpe.200kb.bins <- lapply(bedpe.files, read_tsv) %>% suppressMessages()
names(compAB.bedpe.200kb.bins)	 <- basename(bedpe.files)
compAB.bedpe.200kb.bins$ccl4_10w_norm_corrected_200kb.intra.bedpe.compAB.bedpe


par(oma =c(4,4,2,1), mar = c(4,4,2,1))

#' ### B-B interaction frequency
test1 <- compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe %>%
	mutate(distance = abs(start2 -start1)) %>% 
	filter( seq1 == 'chrX' & seq2 == 'chrX') %>% 
	filter(bin1.comp == 'B' & bin2.comp == 'B') %>%
	mutate(ranges = cut(distance, 100)) %>%
	group_by(ranges) %>%
	summarise(sums = sum(read_count), n =n() ) %>%
	mutate(average_count = sums/(sum(sums)) )

test1 %>% tail	
#plot(test1$average_count, las =2, ylim = c(0, 0.005), type = 'n')	



#+ fig.width = 8, fig.height = 6
plot(test1$average_count, ylim = c(0, 0.01), type ='n', col = 'red', lwd =2.5,
	 main = 'B-B interaction frequencys & distance (con_6w)',
	 xlab = 'chromosomes are split into 100 bins')

for (i in unique(compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe$seq1)[1:19]) {
	test <- compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe %>%
		mutate(distance = abs(start2 -start1)) %>% 
		filter( seq1 == i & seq2 == i) %>% 
		filter(bin1.comp == 'B' & bin2.comp == 'B') %>%
		mutate(ranges = cut(distance, 100)) %>%
		group_by(ranges) %>%
		summarise(sums = sum(read_count), n =n() ) %>%
		mutate(average_count = sums/(sum(sums)) )
	
	points(test$average_count %>% jitter, las =2, cex = 0.8, pch =21, col = 'blue', bg = 'grey')	
	
}

points(test1$average_count, ylim = c(0, 0.01), type ='l', col = 'red', lwd =2.5)
points(test1$average_count, ylim = c(0, 0.01), type ='p', pch =21, col = 'red', cex =1.1, bg ='orange')



#' ### **B-B  interactions**   
#' 
#+ fig.width = 8, fig.height = 6
for (i in names(compAB.bedpe.200kb.bins)[c(6,3,4,1,5,2)]) {
		test1 <- compAB.bedpe.200kb.bins[[i]] %>%
			mutate(distance = abs(start2 -start1)) %>% 
			filter( seq1 == 'chrX' & seq2 == 'chrX') %>% 
			filter(bin1.comp == 'B' & bin2.comp == 'B') %>%
			mutate(ranges = cut(distance, 100)) %>%
			group_by(ranges) %>%
			summarise(sums = sum(read_count), n =n() ) %>%
			mutate(average_count = sums/(sum(sums)) )
		
				
#		print(comp)
		group_info <- gsub(pattern = '_norm_corrected_200kb.intra.bedpe.compAB.bedpe',
						   replacement = '', x = i)
		main_text <- sprintf('B_B interaction frequencys & distance (%s)', group_info)
		
		
		plot(test1$average_count, ylim = c(0, 0.01), type ='n', col = 'red', lwd =2.5,
			 main = main_text,
			 xlab = 'chromosomes are split into 100 bins')
		
		for (chr in unique(compAB.bedpe.200kb.bins[[i]]$seq1)[1:19]) {
			test <- compAB.bedpe.200kb.bins[[i]] %>%
				mutate(distance = abs(start2 -start1)) %>% 
				filter( seq1 == chr & seq2 == chr) %>% 
				filter(bin1.comp == 'B' & bin2.comp == 'B') %>%
				mutate(ranges = cut(distance, 100)) %>%
				group_by(ranges) %>%
				summarise(sums = sum(read_count), n =n() ) %>%
				mutate(average_count = sums/(sum(sums)) )
			
			points(test$average_count %>% jitter, las =2, cex = 0.8, pch =21, col = 'blue', bg = 'grey')	
			
		}
		
		points(test1$average_count, ylim = c(0, 0.01), type ='l', col = 'red', lwd =2.5)
		points(test1$average_count, ylim = c(0, 0.01), type ='p', pch =21, col = 'red', cex =1.1, bg ='orange')
		
	}











#' ### **A-A interaction frequency**  

test2 <- compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe %>%
	mutate(distance = abs(start2 -start1)) %>% 
	filter( seq1 == 'chrX' & seq2 == 'chrX') %>% 
	filter(bin1.comp == 'A' & bin2.comp == 'A') %>%
	mutate(ranges = cut(distance, 100)) %>%
	group_by(ranges) %>%
	summarise(sums = sum(read_count), n =n() ) %>%
	mutate(average_count = sums/(sum(sums)) )

test2 %>% tail	

write_csv(test2, '~/cluster_home/project/liverCancer/scripts/figures_update/ext_fig4/data/ext_fig4h.csv')


test2 <- compAB.bedpe.200kb.bins$ccl4_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe %>%
  mutate(distance = abs(start2 -start1)) %>% 
  filter( seq1 == 'chrX' & seq2 == 'chrX') %>% 
  filter(bin1.comp == 'A' & bin2.comp == 'A') %>%
  mutate(ranges = cut(distance, 100)) %>%
  group_by(ranges) %>%
  summarise(sums = sum(read_count), n =n() ) %>%
  mutate(average_count = sums/(sum(sums)) )

test2 %>% tail	

write_csv(test2, '~/cluster_home/project/liverCancer/scripts/figures_update/ext_fig4/data/ext_fig4i.csv')






#+ fig.width = 8, fig.height = 6
plot(test2$average_count, ylim = c(0, 0.01), type ='n', col = 'red', lwd =2.5,
	 main = 'A-A interaction frequencys & distance (con_6w)',
	 xlab = 'chromosomes are split into 100 bins')

for (i in unique(compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe$seq1)[1:19]) {
	test <- compAB.bedpe.200kb.bins$con_6w_norm_corrected_200kb.intra.bedpe.compAB.bedpe %>%
		mutate(distance = abs(start2 -start1)) %>% 
		filter( seq1 == i & seq2 == i) %>% 
		filter(bin1.comp == 'A' & bin2.comp == 'A') %>%
		mutate(ranges = cut(distance, 100)) %>%
		group_by(ranges) %>%
		summarise(sums = sum(read_count), n =n() ) %>%
		mutate(average_count = sums/(sum(sums)) )
	
	points(test$average_count %>% jitter, las =2, cex = 0.8, pch =21, col = 'blue', bg = 'grey')	
	
}

points(test2$average_count, ylim = c(0, 0.01), type ='l', col = 'red', lwd =2.5)
points(test2$average_count, ylim = c(0, 0.01), type ='p', pch =21, col = 'red', cex =1.1, bg ='orange')


  
#' 
#+ fig.width = 8, fig.height = 6
for (i in names(compAB.bedpe.200kb.bins)[c(6,3,4,1,5,2)]) {
	test2 <- compAB.bedpe.200kb.bins[[i]] %>%
		mutate(distance = abs(start2 -start1)) %>% 
		filter( seq1 == 'chrX' & seq2 == 'chrX') %>% 
		filter(bin1.comp == 'A' & bin2.comp == 'A') %>%
		mutate(ranges = cut(distance, 100)) %>%
		group_by(ranges) %>%
		summarise(sums = sum(read_count), n =n() ) %>%
		mutate(average_count = sums/(sum(sums)) )
	
	
	#		print(comp)
	group_info <- gsub(pattern = '_norm_corrected_200kb.intra.bedpe.compAB.bedpe',
					   replacement = '', x = i)
	main_text <- sprintf('A_A interaction frequencys & distance (%s)', group_info)
	
	
	plot(test1$average_count, ylim = c(0, 0.01), type ='n', col = 'red', lwd =2.5,
		 main = main_text,
		 xlab = 'chromosomes are split into 100 bins')
	
	for (chr in unique(compAB.bedpe.200kb.bins[[i]]$seq1)[1:19]) {
		test <- compAB.bedpe.200kb.bins[[i]] %>%
			mutate(distance = abs(start2 -start1)) %>% 
			filter( seq1 == chr & seq2 == chr) %>% 
			filter(bin1.comp == 'A' & bin2.comp == 'A') %>%
			mutate(ranges = cut(distance, 100)) %>%
			group_by(ranges) %>%
			summarise(sums = sum(read_count), n =n() ) %>%
			mutate(average_count = sums/(sum(sums)) )
		
		points(test$average_count %>% jitter, las =2, cex = 0.8, pch =21, col = 'blue', bg = 'grey')	
		
	}
	
	points(test1$average_count, ylim = c(0, 0.01), type ='l', col = 'red', lwd =2.5)
	points(test1$average_count, ylim = c(0, 0.01), type ='p', pch =21, col = 'red', cex =1.1, bg ='orange')
	
}















