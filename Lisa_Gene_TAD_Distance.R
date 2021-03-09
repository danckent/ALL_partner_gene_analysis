library( GenomicRanges)
library( tibble)
library( plyranges)

Genes       <- read.table( "/home/dankent/Data/HG_Esembl98_Genes.tab", header = T, sep = "\t")
Genes$Chromosome.scaffold.name <- paste0( 'chr',Genes$Chromosome.scaffold.name)
Genes$width <- Genes$Gene.end..bp.- Genes$Gene.start..bp.
Genes$mid_point <- (Genes$width/2)+Genes$Gene.start..bp.
HG <- Genes$HGNC.symbol
Lisa_genes <- as.factor(c("MYC", "BCL6", "NFKB2", "BCL3", "BCL11A", "CCND3", "BCL10", "CHST11", "LAPTM5",
                            "CEBPG", "BCL9","CEBPB", "DDX6","MALT1", "RHOH", "PAFAH1B2", "PCSK7", "EPOR", "CCND1", "FOXP1",
                           "CDK6", "FCGR2B", "PAX5", "IRF4", "MUC1", "BCL2", "CCND2", "ETV6", "MAF", "CRLF2", 
                           "CEBPA", "CEBPD", "IL3", "CEBPE", "ID4", "NBEAP1", "LHX4", "FGFR3", "MAFB"))

NSD2 <-  Genes[Genes$HGNC.symbol=="WHSC1",]
IGL <- Genes[Genes$Gene.start..bp.==22380474,]

Lisa_genes_match <- na.omit(match(Lisa_genes, HG))
Lisa_genes_match <- Genes[Lisa_genes_match,]
Lisa_genes_match <- rbind(Lisa_genes_match, NSD2, IGL)
Lisa_genes_match_DF <- Lisa_genes_match 
Lisa_genes_match <- makeGRangesFromDataFrame(Lisa_genes_match, keep.extra.columns = TRUE, ignore.strand = TRUE, start.field = "Gene.start..bp.", end.field = "Gene.end..bp.", seqnames.field = "Chromosome.scaffold.name")

TADs_BI_GM12878 <- read.table("/home/dankent/Data/GM12878/TAD_domains/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt", header = TRUE)
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

##### RANDOM RANGES ##### 
#Global background model:
#SEs & BDs are relocated in the genome This script generates a matrix with the random starts and random chromosomes for all SEs

library( "regioneR" )
library( "BSgenome.Hsapiens.UCSC.hg19.masked" )
library( "gtools" )
library( "seqbias")


setwd( "/home/dankent/Data/ALL_partner_genes/")
if ( !file.exists("partner_gene_global_random_intervals")) dir.create( file.path( ".", "partner_gene_global_random_intervals"))
if ( !file.exists("partner_gene_global_random_intervals/output"))  dir.create( file.path( "partner_gene_global_random_intervals", "output"))

input_dir     <- "partner_gene_global_random_intervals/"
output_dir    <- "partner_gene_global_random_intervals/output/"
total_permuts <- 1000

#Datasets
both_sets <- c( "Lisa_genes_match")

## Retrieving genome and mask
# Get genome and mask
human.genome    <- getGenomeAndMask( genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$genome

# Filter genome and mask
human.mask <- getGenomeAndMask( genome = "BSgenome.Hsapiens.UCSC.hg19.masked")$mask

#  dataset preparation
# Filtering dataset
for( set in both_sets) {
  set_ranges <- get( set)
  # Pass all ranges to UCSC seqname style
  seqlevelsStyle( set_ranges)      <- "UCSC"
}

# CREATE GRanges WITH AVAILABLE REGIONS FOR RANDOMIZATION (what is not low-mappability)
# unmasked genome
masked_coverage <- coverage( append( human.genome, reduce( human.mask)))

coverages_RANGES <- GRanges( )
for( i in paste0( "chr", 1:22)) {
  coverage_Rle <- masked_coverage[[i]]
  chr_RANGES <- GRanges( seqnames = as.character( i), ranges= ranges( coverage_Rle), coverage = coverage_Rle@values)
  coverages_RANGES <- append( coverages_RANGES, chr_RANGES)
}

# Coverage 2 when there is mask 
# Coverage 1 when no mask : available regions for  randomization
unmasked_genome <- coverages_RANGES[ coverages_RANGES$coverage == 1]
unmasked_genome_for_randomization <- unmasked_genome
unmasked_genome_for_randomization$original_chr <- seqnames( unmasked_genome)

# # Transform chromosomes from chr---to number
seqlevelsStyle( unmasked_genome_for_randomization) <- "NCBI"
chrs_in_unmasked_genome <- factor( as.character( seqnames( unmasked_genome_for_randomization) ),                               
                                   levels = 1:length( unmasked_genome_for_randomization) )

# One seqlevel per range 
seqlevels( unmasked_genome_for_randomization) <- as.character( 1:length( unmasked_genome_for_randomization) )
seqnames( unmasked_genome_for_randomization)  <- 1:length( unmasked_genome_for_randomization)
names( unmasked_genome_for_randomization)     <- as.character( seqnames( unmasked_genome_for_randomization) )

#  low mappability needs to remain 

## Randomization function
global_randomization <- function( set_map, seed) {
  
  set.seed( seed)
  random_int0 <- random.intervals( unmasked_genome_for_randomization,
                                   n=length( set_map),
                                   ms=width( set_map)-1)
  
  # This function has crated a randomized set giving the coordinates within the range. For example,
  # a CNV starting in position 2500 from  a range that starts with 2000 will be given as start  position 500.
  # Needs to be fixed. And also recover original seqnames
  random_int <- GRanges( seqnames = unmasked_genome_for_randomization [ as.character( seqnames( random_int0) ) ]$original_chr,
                         ranges = IRanges( start = start( random_int0) + start( unmasked_genome_for_randomization[ seqnames( random_int0) ] ),
                                           end   = end( random_int0) + start( unmasked_genome_for_randomization[ seqnames( random_int0) ] ) ) )
  
  random_int
}


for( map_num in 1) {
  set <- both_sets[ map_num]
  set_ranges <- get( set)
  
  all_random_starts <- c( )
  all_random_chrs   <- c( )
  
  for( i in 1:total_permuts) {
    # Randomize
    random_ranges <- global_randomization( set_map = set_ranges, seed = map_num*1000000 + i)
    # We will save all random starts in a matrix
    all_random_starts <- cbind( all_random_starts, start( random_ranges))
    # We will save all random chromosomes in a matrix
    all_random_chrs   <- cbind( all_random_chrs,   as.character( seqnames( random_ranges) ) )
    if( i%%500 == 0) print( i)
  }
  
  write.table( all_random_starts, quote = F, row.names = F, col.names = F,
               file = paste0( output_dir, "partner_gene_global_random_intervals_random_starts_", set, "_",
                              total_permuts, ".txt"))
  rm( all_random_starts)  
  
  write.table( all_random_chrs, quote = F, row.names = F, col.names = F,
               file = paste0( output_dir, "partner_gene_global_random_intervals_random_chrs_", set, "_",
                              total_permuts, ".txt"))
  rm( all_random_chrs)
  
  save( set_ranges, 
        file = paste0( output_dir, "partner_gene_global_random_intervals_original_", set, "_GRanges.RData"))
  
  
  print( set)
}

save( both_sets, total_permuts, file = "/home/dankent/Data/ALL_partner_genes/partner_genes.RData")

#### RANDOM PARTNER GENE RANGES IN GENES #######
load("/home/dankent/Data/ALL_partner_genes/partner_gene_global_random_intervals/output/Lisa_genes_match_1000permuts.RData")

## Loading RANDOM STARTS and CHROMOSOMES
chromosomes <- read.table("~/Data/ALL_partner_genes/partner_gene_global_random_intervals/output/partner_gene_global_random_intervals_random_chrs_Lisa_genes_match_1000.txt", stringsAsFactors = F, header = F)
print("chromosomes.loaded")

starts <- read.table("~/Data/ALL_partner_genes/partner_gene_global_random_intervals/output/partner_gene_global_random_intervals_random_starts_Lisa_genes_match_1000.txt", stringsAsFactors = F, header = F)
print("starts.loaded")

# original_ranges
load("partner_gene_global_random_intervals/output/partner_gene_global_random_intervals_original_Lisa_genes_match_GRanges.RData")

print("original_ranges.loaded")

# Random values
ran  <- c()

while(ncol(chromosomes) > 0) {
  random <- GRanges(seqnames = chromosomes[, 1],
                    ranges   = IRanges(start = starts[, 1],
                                       width = width( set_ranges)))

  ran <- c(ran, nearest(random, TADs_GM12878))
  
  chromosomes <- chromosomes[, -1]
  starts      <- starts[, -1]
  
  if(class(starts) == "integer") {
    chromosomes <- matrix(chromosomes, ncol = 1)
    starts      <- matrix(starts, ncol = 1)
  }
  if(ncol(chromosomes)%%100 == 0) print(ncol(chromosomes))
}

save(ran,
      file = paste0( "partner_gene_global_random_intervals/output/", 
                     set, "_",
                     total_permuts, "permuts.RData"))


### calculate distance for random ####
Nearest_Random_TAD <- TADs_GM12878[ran]
Nearest_Random_TAD_df <- data.frame(Nearest_Random_TAD)

#creat mid-point for all random genes 
#re-load starts
starts <- read.table("~/Data/ALL_partner_genes/partner_gene_global_random_intervals/output/partner_gene_global_random_intervals_random_starts_Lisa_genes_match_1000.txt", stringsAsFactors = F, header = F)
midpoint_all <- starts+(Lisa_genes_match_DF$width/2)
#all in one row
midpoint_all <- data.frame(unlist(midpoint_all))

Distance_Random_mid_to_TAD <- data.frame(cbind(abs(midpoint_all - Nearest_Random_TAD_df$startTAD),abs(midpoint_all - Nearest_Random_TAD_df$endTAD)))
Distance_Random_mid_to_TAD <- apply(Distance_Random_mid_to_TAD,1,min)

#mean average 
median(Distance_Random_mid_to_TAD)

Distance_Random_mid_to_TAD_df <- data.frame(Distance_Random_mid_to_TAD)

### calculate distance for actual genes ####
obs <- nearest(set_ranges, TADs_GM12878)
obs <- TADs_GM12878[obs] 

Lisa_genes_match_DF <- data.frame(Lisa_genes_match)

Distance_mid_to_TAD <- data.frame(cbind(abs(Lisa_genes_match_DF$mid_point - obs$startTAD),abs(Lisa_genes_match_DF$mid_point- obs$endTAD)))
Distance_mid_to_TAD <- apply(Distance_mid_to_TAD,1,min)

#mean 
median(Distance_mid_to_TAD)

Distance_mid_to_TAD_df <- data.frame(Distance_mid_to_TAD)

#### Visualise

boxplot(Distance_mid_to_TAD_df$Distance_mid_to_TAD, Distance_Random_mid_to_TAD_df$Distance_Random_mid_to_TAD, outline=FALSE, ylab = "Distance between the mid point of the gene and the start/end of the TAD", xlab = "Observed                                                                                      Random")

save(Distance_mid_to_TAD_df, Distance_Random_mid_to_TAD_df, file = "/home/dankent/Data/GM12878/Partner Gene Proximity TAD Boundary/for_Marco.RData")
