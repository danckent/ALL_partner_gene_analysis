#washU ready files
library(dplyr)
#had to edit this for each  'select("bait", "otherend", "tB")' section. Can't figure out a way to select the final colum (like you would with tail(1)).

washU_formatn <- function(PCHiC_file){
  y <- read.table(file = PCHiC_file, header = T)
  y$bait <- paste(y$baitChr, y$baitStart, y$baitEnd, sep = ",")
  y$otherend <- paste(y$oeChr, y$oeStart, y$oeEnd, sep = ",")
  y <- y %>%
    select("bait", "otherend", "nB")
  print(y[grep(pattern = "MT",x = y$otherend),])
  y <- y[-grep(pattern = "MT",x = y$otherend),]
  readr::write_tsv(x = y, file = paste(basename(PCHiC_file), ".bed"), col_names = F)
}

washU_formatt <- function(PCHiC_file){
  y <- read.table(file = PCHiC_file, header = T)
  y$bait <- paste(y$baitChr, y$baitStart, y$baitEnd, sep = ",")
  y$otherend <- paste(y$oeChr, y$oeStart, y$oeEnd, sep = ",")
  y <- y %>%
    select("bait", "otherend", "tB")
  print(y[grep(pattern = "MT",x = y$otherend),])
  y <- y[-grep(pattern = "MT",x = y$otherend),]
  readr::write_tsv(x = y, file = paste(basename(PCHiC_file), ".bed"), col_names = F)
}

nB <- washU_formatn(PCHiC_file =  "~/Data/ALL_partner_gene/PCHi-C/nB_interactions.txt")
tB <- washU_formatt(PCHiC_file =  "~/Data/ALL_partner_gene/PCHi-C/tB_interactions.txt")

setwd(dir = "~/Data/ALL_partner_gene/PCHi-C/hematopoietic_precursors/")

washU_format_heam_precursors <- function(ibed_file){
  y <- read.table(file = ibed_file, header = T)
  y$bait <- paste(y$bait_chr, y$bait_start, y$bait_end, sep = ",")
  y$otherend <- paste(y$otherEnd_chr, y$otherEnd_start, y$otherEnd_end, sep = ",")
  y <- y %>%
    select("bait", "otherend", "score")
  y$bait <- paste0( 'chr',y$bait)
  y$otherend <- paste0( 'chr',y$otherend)
  readr::write_tsv(x = y, file = paste(basename(ibed_file),".bed"), col_names = F)
}

washU_format_heam_precursors(ibed_file = "~/Data/ALL_partner_gene/PCHi-C/hematopoietic_precursors/ProB_WT_merge_cutoff_5.ibed")
washU_format_heam_precursors(ibed_file = "~/Data/ALL_partner_gene/PCHi-C/hematopoietic_precursors/PreB_WT_merge_cutoff_5.ibed")
washU_format_heam_precursors(ibed_file = "~/Data/ALL_partner_gene/PCHi-C/hematopoietic_precursors/CMP_WT_merge_cutoff_5.ibed")
washU_format_heam_precursors(ibed_file = "~/Data/ALL_partner_gene/PCHi-C/hematopoietic_precursors/HSC_WT_merge_cutoff_5.ibed")

washU_format_loops <- function(Chrom_loops_file){
  y <- read.table(file = Chrom_loops_file, header = F)
  y$A <- paste(y$V1, y$V2, y$V3, sep = ",")
  y$B <- paste(y$V4, y$V5, y$V6, sep = ",")
  y <- y %>%
    select("A", "B", "V7")
  readr::write_tsv(x = y, file = paste(basename(Chrom_loops_file), ".bed"), col_names = F)
}

washU_format_loops(Chrom_loops_file = "Rao_2014.GM12878.hg19.peakachu-merged.loops")
