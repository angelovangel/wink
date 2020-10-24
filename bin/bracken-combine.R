#!/usr/bin/env Rscript

# arg1 is a folder containing bracken output files (tsv) with results from one 4k ONT fastq file
# arg2 is the name of the returned tsv file
# returns a tibble with data summarized from all files
require(dplyr)
require(data.table)

arg <- commandArgs(trailingOnly = TRUE)

files <- list.files(arg[1], pattern = "*.tsv", full.names = TRUE)
dfs <- lapply(files, fread)

bind_rows(dfs) %>% 
	group_by(name, taxonomy_id) %>% 
	summarise_at(c("kraken_assigned_reads", "new_est_reads"), sum, na.rm = TRUE) %>% 
	ungroup() %>% mutate(freq = new_est_reads/sum(new_est_reads, na.rm = TRUE)) %>% 
	arrange(desc(new_est_reads)) %>%
	fwrite(file = arg[2], sep = "\t")
