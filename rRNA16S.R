library(seqinr)
library(optparse)
path <- list.files(path= "~", full.names=TRUE,
                   recursive=TRUE,pattern="(Func.R)")
pathTrue <- path[grep(pattern = "Function/Func.R", path, fixed = T)]
source(pathTrue) #path to Func.R file
#func--------------------------------------------------------------
modTable <- function(taxa_16S){
  for (i in 1:ncol(taxa_16S)) {
    taxa_16S[which(taxa_16S[,i]==""),i] <- NA
  }
  name <- strsplit(taxa_16S$V1, "_", fixed = TRUE)
  name2 <- unlist(lapply(name, "[", 1))
  taxa_16S$V1 <- name2
  taxa_16S <- taxa_16S[,-c(2,3)]
  taxaA <- taxa_16S[which(taxa_16S$V4=="rootrank"),]
  taxaA <- taxaA[,-c(2,4,6,8,10,12,14)]
  colnames(taxaA) <- c("subject",	"superkingdom",	"phylum",	"class",	"order",	"family",	"genus")
  taxaB <- taxa_16S[which(taxa_16S$V4=="Lineage=Root"),]
  taxaB <- taxaB[,-c(2,3,5,7,9,11,13)]
  colnames(taxaB) <- c("subject",	"superkingdom",	"phylum",	"class",	"order",	"family", "genus")
  taxaC <- taxa_16S[which(taxa_16S$V4=="Bacteria"),]
  taxaC <- taxaC[,-c(3,5,7,9,11,13,14)]
  colnames(taxaC) <- c("subject",	"superkingdom",	"phylum",	"class",	"order",	"family", "genus")
  taxa16s <- rbind(taxaA, taxaB, taxaC)
  return(taxa16s)
}


#Options------------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-r", "--blastRDP"), action = "store", type = "character", default = "16sContigvsRDP.txt",
              help = "output contigs with 16S vs RDP in path/filename.ext - [default %default]"),
  make_option(c("-t", "--taxaRDP"), action = "store", type = "character", default = "RDP16s_taxa_mod.txt",
              help = "RDP subject taxonomy path/filename.ext - [default %default]"),
  make_option(c("-l", "--minContigLen"), action = "store", type = "integer", default = 1000, 
              help = "minimum contig length considered - [default %default]"),
  make_option(c("-s", "--min16SLen"), action = "store", type = "integer", default = 500, 
              help = "minimum 16S gene length considered - [default %default]")
  )
opt <- parse_args(OptionParser(option_list = option_list))

#Tables----------------------------------------------------------------------------------------------------------
rRNA16S <- read.table(opt$blastRDP,
                      col.names = c("Contig",	"subject",	"X.id",	"ali_Length",	"X1",	"X2",	
                                    "STARTQ",	"ENDQ",	"STARTS",
                                    "ENDS",	"EVALUE",	"BITSCORE",	"QLEN",	"SLEN",	"QCOVS",	
                                    "GAPS",	"QCOVHSP"),
                      header= FALSE, stringsAsFactors = FALSE)

taxa_rRNA16S <- read.table(opt$taxaRDP,
                           header= FALSE, stringsAsFactors = FALSE,
                           sep = "\t", fill = T, comment.char = "")
taxa16s <- modTable(taxa_rRNA16S)
#16S-------------------------------------------------------------------------------------------------------------
rRNA16sTaxonomy <- merge(rRNA16S, taxa16s, by = "subject")
length <- unlist(lapply(strsplit(rRNA16sTaxonomy$Contig,split = "_"), "[[", 4))
rRNA16sTaxonomy$length <- length
rRNA16sTaxonomy <- rRNA16sTaxonomy[which(rRNA16sTaxonomy$length>=opt$minContigLen),]

rRNA16sTaxonomy <- rRNA16sTaxonomy[which(rRNA16sTaxonomy$ali_Length>opt$min16SLen),]
rRNA16sTaxonomy2 <- BestTaxa(rRNA16sTaxonomy)

#rRNA16sTaxonomy3 <- rRNA16sTaxonomy2[which(rRNA16sTaxonomy2$ali_Length>600),]
write.table(rRNA16sTaxonomy2, "rRNA16sTaxonomy2.txt", sep = "\t", quote = F, row.names = F)





