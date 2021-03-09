library( GenomicRanges)
library( tibble)
library( plyranges)
library( dplyr)

Genes       <- read.table( "/home/dankent/Data/HG_Esembl98_Genes.tab", header = T, sep = "\t")
Want                              <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X") # no TADs on the Y chromosome in the data 
Genes                             <- subset(Genes, Genes$Chromosome.scaffold.name %in% Want)
Genes$Chromosome.scaffold.name <- paste0( 'chr',Genes$Chromosome.scaffold.name)
Genes$width <- Genes$Gene.end..bp.- Genes$Gene.start..bp.
Genes$mid_point <- (Genes$width/2)+Genes$Gene.start..bp.
HG <- Genes$HGNC.symbol

B_cell <- as.factor(c("MYC", "CEBPG", "CEBPB", "EPOR", "DDX6", "BCL2", "ETV6", "CRLF2", "CEBPA",
                                     "CEBPD","IL3", "CEBPE", "ID4", "LHX4", "MIR125B1")) 

#removed BCL9 
Lisa_genes <- as.factor(c("MYC", "BCL6", "NFKB2", "BCL3", "BCL11A", "CCND3", "BCL10", "CHST11", "LAPTM5",
                          "CEBPG","CEBPB", "DDX6","MALT1", "RHOH", "PAFAH1B2", "PCSK7", "EPOR", "CCND1", "FOXP1",
                          "CDK6", "FCGR2B", "PAX5", "IRF4", "MUC1", "BCL2", "CCND2", "ETV6", "MAF", "CRLF2", 
                          "CEBPA", "CEBPD", "IL3", "CEBPE", "ID4", "NBEAP1", "LHX4", "FGFR3", "MAFB"))

Lisa_genes <- Lisa_genes[!Lisa_genes %in% B_cell]

NSD2 <-  Genes[Genes$HGNC.symbol=="WHSC1",]
IGL <- Genes[Genes$Gene.start..bp.==22380474,]

Lisa_genes_match <- na.omit(match(Lisa_genes, HG))
Lisa_genes_match <- Genes[Lisa_genes_match,]
Lisa_genes_match <- rbind(Lisa_genes_match, NSD2, IGL)
Lisa_genes_match <- makeGRangesFromDataFrame(Lisa_genes_match, keep.extra.columns = TRUE, ignore.strand = TRUE, start.field = "Gene.start..bp.", end.field = "Gene.end..bp.", seqnames.field = "Chromosome.scaffold.name")

TADs_BI_GM12878 <- read.table("/home/dankent/Data/GM12878/TAD_domains/MARCO_GM_TADs.bed", header = TRUE)
TADs_BI_GM12878 <- TADs_BI_GM12878[ ,c(1,2,3,8)]
colnames(TADs_BI_GM12878) <- c("chromosome", "start", "end", "corner_score")
TADs_BI_GM12878$chromosome<- paste0( 'chr',TADs_BI_GM12878$chromosome)

TADs_process_for_overlap <- function( input_variable){
  
  input_variable              <- tibble::rowid_to_column( input_variable, "ID")
  input_variable$widthTAD     <- ( input_variable$end - input_variable$start)
  input_variable$startTAD     <- input_variable$start
  input_variable$endTAD       <- input_variable$end
  makeGRangesFromDataFrame( input_variable, keep.extra.columns = TRUE, ignore.strand = TRUE)
}

TADs_GM12878 <- TADs_process_for_overlap(TADs_BI_GM12878)
TADs_GM12878            <- data.frame(reduce_ranges(TADs_GM12878))
TADs_GM12878$widthTAD   <- ( TADs_GM12878$end - TADs_GM12878$start)
TADs_GM12878$startTAD   <- TADs_GM12878$start
TADs_GM12878$endTAD     <- TADs_GM12878$end
TADs_GM12878            <- tibble::rowid_to_column( TADs_GM12878, "ID")
TADs_GM12878            <- makeGRangesFromDataFrame( TADs_GM12878, keep.extra.columns = T )
##### get random genes and split into 500 chunks of 41 genes (same length as original) #####

real_random_ranges <- sample_n(tbl = Genes, size = 20500)
real_random_ranges <- split(real_random_ranges, sample(rep(500)))

#since all data is now help in nested lists, I build fucntions with basic functions within, to apply across the whole list
#create Granges 
real_random_ranges_to_GRanges_func <- function(list){
  makeGRangesFromDataFrame(list, keep.extra.columns = TRUE, ignore.strand = TRUE, start.field = "Gene.start..bp.", end.field = "Gene.end..bp.", seqnames.field = "Chromosome.scaffold.name")
}
GR_real_random_ranges <- lapply(real_random_ranges,real_random_ranges_to_GRanges_func)

#find nearest TAD to each gene
nearest_for_list <- function(list){GenomicRanges::nearest(list, TADs_GM12878)}
near <- lapply(GR_real_random_ranges, nearest_for_list)

#identify the nearest TAD using the above
Tad_Extract <- function(list){
  TADs_GM12878[list]
}
TAD_match_random <- lapply(near, Tad_Extract)

#create data frames
TAD_match_random <- lapply(TAD_match_random, data.frame)
real_random_ranges <- lapply(real_random_ranges, data.frame)

#extract important information
midpoint_all  <- unlist(lapply(real_random_ranges, `[`, 8))
TADstart      <- unlist(lapply(TAD_match_random, `[`, 8))
TADend        <- unlist(lapply(TAD_match_random, `[`, 9))

#calculate the absolute distance and find the minimum distance (either TAD start of TAD end)
Distance_Random_mid_to_TAD <- data.frame(cbind(abs(midpoint_all - TADstart),abs(midpoint_all - TADend)))
Distance_Random_mid_to_TAD <- apply(Distance_Random_mid_to_TAD,1,min)

#mean average 
median(Distance_Random_mid_to_TAD)

Distance_Random_mid_to_TAD_df <- data.frame(Distance_Random_mid_to_TAD)

############################################
### calculate distance for actual genes ####
############################################

TAD_match <- TADs_GM12878[nearest(Lisa_genes_match, TADs_GM12878)]

Distance_mid_to_TAD <- data.frame(cbind(abs(Lisa_genes_match$mid_point - TAD_match$startTAD),abs(Lisa_genes_match$mid_point - TAD_match$endTAD)))
Distance_mid_to_TAD <- apply(Distance_mid_to_TAD,1,min)

median(Distance_mid_to_TAD)

Distance_mid_to_TAD_df <- data.frame(Distance_mid_to_TAD)

Clean_up <- cbind(data.frame(Lisa_genes_match), Distance_mid_to_TAD_df, TAD_match)
Clean_up <- Clean_up[c(1:3,7,9:10,12,13)]
Clean_up %>% arrange(Distance_mid_to_TAD)
#### Visualise

boxplot(Distance_mid_to_TAD_df$Distance_mid_to_TAD, Distance_Random_mid_to_TAD_df$Distance_Random_mid_to_TAD, outline=FALSE, 
        ylab = "Distance from gene mid-point to TAD start/end", xlab = "Observed                 Random", main = "Distance from Gene to TAD boundary", col = c("#E69F00", "#56B4E9"))


