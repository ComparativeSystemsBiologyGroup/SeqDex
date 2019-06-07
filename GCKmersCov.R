path <- list.files(path= "~", full.names=TRUE,
                   recursive=TRUE,pattern="(Func.R)")
pathTrue <- path[grep(pattern = "Function/Func.R", path, fixed = T)]
source(pathTrue) #path to Func.R file

library(optparse)
library(seqinr)
#Options--------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-c", "--contigs"), action = "store", type = "character", default = NA,
              help = "path/filename.ext of contigs file - [default detect fasta file in folder]"),
  make_option(c("-K", "--Kmers"), action = "store", type = "integer", default = 3,
              help = "Kmers length to be used - [default %default]"),
  make_option(c("-s", "--covSingle"), action = "store", type = "character", default = "Coverage/onlymapping_sorted_SINGLE.bed",
              help = "path/filename.ext coverage file of single paired reads - [default %default]"),
  make_option(c("-p", "--covPaired"), action = "store", type = "character", default = "Coverage/onlymapping_sorted_PAIRED.bed",
              help = "path/filename.ext coverage file of paired reads - [default %default]"),
  make_option(c("-T", "--threads"), action = "store", type = "integer", default = 8,
              help = "number of threads to be used - [default %default]")
  )
opt <- parse_args(OptionParser(option_list = option_list))


#Tables-----------------------------------------------------------------------------------------------------
if(is.na(opt$contigs)){
  stop("--contigs is mandatory")
  } else {
    contigs <- read.fasta(file = opt$contigs, seqtype = "DNA")
  }

covfileSingle<- read.table(opt$covSingle, header=FALSE, stringsAsFactors = FALSE)

covfilePaired<- read.table( opt$covPaired, header=FALSE, stringsAsFactors = FALSE)
#USER defined parameters --------------------------------------------------------------------
KLEN=opt$Kmers
#Kmers count----------------------------------------------------------------------------------
gck <- KmersCountNorm(KLEN, contigs, cluster = opt$threads)
gckMult <- cbind(gck[,1:2], gck[,3:ncol(gck)]*1000)
#write.table(gckMult, "gcx.txt", sep = "\t")


#Coverage-------------------------------------------------------------------------------
covTableSum <- data.frame(Contigs=covfileSingle[,1],
                          TotalCoverage=(covfileSingle$V4+(covfilePaired$V4/2)),
                          length = covfileSingle$V3,
                          covLen = (covfileSingle$V4+(covfilePaired$V4/2))/covfileSingle$V3)
gckCov <- gckMult[gckMult$Contig %in% covTableSum$Contig,] #only contig with coverage
covOrd <- covTableSum[match(gckCov$Contig, covTableSum$Contig),] #ordering the contig

gckCovTable <- data.frame(gckCov, covLen = covOrd[,4], length = covOrd$length, stringsAsFactors = F)

write.table(gckCovTable, "gckCovTable.txt", sep = "\t", quote = F, row.names = F)
