path <- list.files(path= "~", full.names=TRUE,
                   recursive=TRUE,pattern="(Func.R)")
pathTrue <- path[grep(pattern = "Function/Func.R", path, fixed = T)]
source(pathTrue) #path to Func.R file
suppressPackageStartupMessages(library(optparse, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))
suppressPackageStartupMessages(library(dbscan, quietly = T))
suppressPackageStartupMessages(library(uwot, quietly = TRUE))
suppressPackageStartupMessages(library(igraph, quietly = T))

#Options-----------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-m", "--modelOutput"), action = "store", type = "character", default = "../SVMoutput/1outputSVM.txt", 
              help = "output of SVM/RF prediction in path/file.ext format or \n\t\t\t comma separated list of taxonomy level/s used in ML modeling- no default set "),
  make_option(c("-y", "--modelType"), action = "store", type = "character", default = NA, 
              help = "which machine learning algorithm output are being clustered. \n\t\t\t Only if modelOutput is in taxonomy levels format. Only allowed: SVM or RF (not both) - no default set "),
  make_option(c("-G", "--gcCovKmersTable"), action = "store", type = "character", default = "../Coverage/gckCovTable.txt",
              help = "GCKmersCovDef.R output in path/filename.ext - [default %default]"),
  make_option(c("-r", "--rRNA16S"), action = "store", type = "character", default = "../Taxonomy/rRNA16sTaxonomy2.txt",
              help = "rRNA16S.R output table in path/file.ext - [default %default]"),
  make_option(c("-t", "--taxonomy"), action = "store", type = "character", default = "../Taxonomy/1taxonomyIteration.txt",
              help = "Taxonomy.R output/s in path/file.ext. If multiple iteration performed,\n\t\t\t then use the higher order iteration file [default %default]"),
  make_option(c("-R", "--network"), action = "store", type = "character", default = "../Coverage/contigNet.txt",
              help = "reads based graph network path/filename.ext - [default %default]"),
  make_option(c("-i", "--input"), action = "store", type = "character", default = "Kmers",
              help = "type of input to be used in clusters analysis. \n\t\t\t Possible input type are: \n\t\t\t\t Kmers \n\t\t\t\t Kmers, 'GC' and/or 'cov'. \n\t\t [default %default]"),
  make_option(c("-l", "--minContigLen"), action = "store", type = "integer", default = 1000, 
              help = "minimum contig length considered - [default %default]"),
  make_option(c("-n", "--targetName"), action = "store", type = "character", default = "Betaproteobacteria",
              help = "target to be isolated using umap-dbscan; \n\t\t taxonomic category at class level of the target to be searched on clustering [default %default]"),
  make_option(c("-T", "--threads"), action = "store", type = "integer", default = 8,
              help = "number of threads to be used - [default %default]"),
  make_option(c("-C", "--ncomponents"), action = "store", type = "integer", default = 2,
              help = "number of final dimensions to be obtained with umap dimensions reduction - [default %default]"),
  make_option(c("-E", "--Edges"), action = "store", type = "integer", default = 10,
              help = "minimum edges value to be considered to construct network to extend the SVM predicted taxonomy - [default %default]"),
  make_option(c("-V", "--VerticesDegree"), action = "store", type = "integer", default = 5,
              help = "maximum vertex degree of the network to be considered to extract component and extend the SVM predicted taxonomy - [default %default]"),
  make_option(c("-S", "--componentSize"), action = "store", type = "integer", default = 2,
              help = "minimum size of the components to be considered to extend the SVM predicted taxonomy - [default %default]"),
  make_option(c("-M", "--mixedComponents"), action = "store", type = "double", default = 0.2,
              help = "threshold of frequency of alternative SVM predicted taxonomy over the total of the component. \n\t\t\t Under this frequency value the alternative will be corrected and substituded by majority - [default %default]")
)
#options(error=stop)
opt <- parse_args(OptionParser(option_list = option_list))

#USER defined parameters-----------------------------------------------------------------

if(opt$modelOutput=="NA") {
  stop("ERROR in CLustering.R = model output table need to be provided")
}

if(grepl(opt$modelOutput, pattern = "output")){
  if (grepl(opt$modelOutput, pattern = ",")) {
    stop("Only one prediction/classification file (SVM por RF) is allowed ")
  } else {
    out.file <- as.vector(opt$modelOutput)
  }
  
} else {
  flagEND <- paste(paste("Output", opt$modelType, sep = ""), ".txt", sep = "")
  flagBEG <- paste(paste("../", opt$modelType, sep = ""), "output/", sep = "")
  if (grepl(opt$modelOutput, pattern = ",")) {
    out.file <- tolower(as.character(unlist(strsplit(opt$modelOutput, split = ",", fixed = T))))
    out.file <- out.file[length(out.file)]
    out.file <- paste(flagBEG, paste(out.file, flagEND, sep = ""), sep = "")
  } else {
    out.file <- tolower(as.character(opt$modelOutput)) 
    out.file <- paste(flagBEG, paste(out.file, flagEND, sep = ""), sep = "")
  } 
}
#outputTable <- read.table("../SVMoutput/1outputSVM.txt", sep = "\t", stringsAsFactors = F, header = T)
outputTable <- read.table(out.file, sep = "\t", stringsAsFactors = F, header = T)
taxa.name <- colnames(outputTable)[(ncol(outputTable)-1)]

if(grepl(opt$taxonomy, pattern = "TaxonomyIteration.txt")){
  if (grepl(opt$taxonomy, pattern = ",")) {
    stop("Only one taxonomy file is allowed ")
  } else {
    taxa.file <- as.vector(opt$taxonomy)
  }
} else {
  if (grepl(opt$taxonomy, pattern = ",")) {
    taxa.file <- tolower(as.character(unlist(strsplit(opt$taxonomy, split = ",", fixed = T))))
    taxa.file <- taxa.name[which(taxa.file == taxa.name)]
    taxa.file <- paste("../Taxonomy/", paste(taxa.file, "TaxonomyIteration.txt", sep = ""), sep = "")
  } else {
    taxa.file <- tolower(as.character(opt$taxonomy)) 
    taxa.file <- paste("../Taxonomy/", paste(taxa.file, "TaxonomyIteration.txt", sep = ""), sep = "")
  } 
}



#Data------------------------------------------------------------------------
#for coverage info of the target----------------------------------------------------------
mergedTaxonomy <-read.table(opt$taxonomy, sep = "\t", stringsAsFactors = F, header = T)
gckCovTable <- read.table(opt$gcCovKmersTable, sep = "\t", stringsAsFactors = F, header = T)
gckCovTable <- gckCovTable[which(gckCovTable$length>opt$minContigLen),]

rRNA16sTaxonomy2 <- read.table(opt$rRNA16S, sep = "\t", stringsAsFactors = F, header = T)
rRNA16sTaxonomy2 <- rRNA16sTaxonomy2[which(rRNA16sTaxonomy2$ali_Length>100 ),]
rRNA16sTaxonomy3 <- data.frame(Contig = rRNA16sTaxonomy2[,2],
                               TaxonDensity = 1,
                               taxa.name = rRNA16sTaxonomy2[,which(tolower(colnames(rRNA16sTaxonomy2))==taxa.name)], stringsAsFactors = F)
colnames(rRNA16sTaxonomy3)[3] <- taxa.name
mergedTaxonomy <- rbind.data.frame(mergedTaxonomy,
                                   rRNA16sTaxonomy3, stringsAsFactors = F)
mergedTaxonomy <- unique(mergedTaxonomy)

freq <- as.data.frame(table(mergedTaxonomy[,which(tolower(colnames(mergedTaxonomy))==taxa.name)]))
x <- freq[which(freq$Freq>5),]
mergedTaxonomy <- mergedTaxonomy[which(mergedTaxonomy[,which(tolower(colnames(mergedTaxonomy))==taxa.name)] %in% x$Var1),]

completeTable <- merge(x=gckCovTable, y=mergedTaxonomy, #then merge the two tables
                       by= "Contig", all.x =TRUE) #merge put NA's where there is not the value
completeTable[is.na(completeTable)] <-"NoBlastHit" #Substitute NA's with another variabe name
coverage16s <- merge(rRNA16sTaxonomy2[,c(2:4, which(colnames(rRNA16sTaxonomy2)=="superkingdom"), which(colnames(rRNA16sTaxonomy2)=="class"))], completeTable[,c(1,2,(ncol(completeTable)-3), (ncol(completeTable)-2))], by = "Contig", all.x = T)

#UMAP on Selected target table--------------------------------------------------------------------------------------
#TAXA = which(colnames(completeTableNum)==taxa.name)
COLUMN = which(colnames(outputTable) %in% colnames(outputTable)[-which(colnames(outputTable)=="Contig" |
                                                                           colnames(outputTable)=="GCcont" |
                                                                           colnames(outputTable)=="covLen" |
                                                                           colnames(outputTable)=="length" |
                                                                           colnames(outputTable)=="TaxonDensity" |
                                                                           tolower(colnames(outputTable))==taxa.name |
                                                                           colnames(outputTable)=="percent")])



umapT <- umap(outputTable[,COLUMN], n_threads = opt$threads, n_components = opt$ncomponents)
plot(umapT[,1], umapT[,2])
umapT.layout <- as.data.frame(umapT, stringsAsFactors = F)

if (grepl(opt$input, pattern = ",")) {
  input.type <- unlist(strsplit(opt$input, split = ",", fixed = T))
  for (g in 1:length(input.type)) {
    if (tolower(input.type[g])==tolower("GC")) {
      gcCont <- which(colnames(outputTable)=="GCcont")
    } else if (tolower(input.type[g])==tolower("cov")) {
      coverage <- which(colnames(outputTable)=="covLen")
    } else {
      NULL
    }
  }
  if(get0("gcCont", ifnotfound = F)) {
    umpaT.umapT.layout <- cbind(umapT.layout, outputTable[,gcCont])
    }
  if(get0("coverage", ifnotfound = F)) {
    umpaT.umapT.layout <- cbind(umapT.layout, outputTable[,coverage])
    } else {
        NULL
    }
}

  


#kNNdistplot(umapT.layout, k=round(log(nrow(outputTable))))
y <- kNNdist(umapT.layout, k=round(log(nrow(outputTable))))

eps <- as.data.frame(table(round(y, digits = 1)), stringsAsFactors = F)
eps$delta <- rep(NA, nrow(eps))

for (i in 2:(nrow(eps))) {
  eps[i,3] <-   eps[(i-1),2] - eps[(i),2] 
}

#print(eps[which(eps$delta==max(eps$delta, na.rm = T)),1])
cl.dbscan <- dbscan(umapT.layout, 
                    minPts = round(log(nrow(outputTable))), 
                    eps = as.numeric(eps[which(eps$delta==max(eps$delta, na.rm = T)),1])
                    )

umapT.layout$group <- cl.dbscan$cluster
umapT.layout$Contig <- outputTable$Contig
umapT.layout$percent <- outputTable$percent

ggplot(umapT.layout, aes(x=umapT.layout[,1], y = umapT.layout[,2], 
                         color = as.character(umapT.layout$group)
                         )) + geom_point() +
 scale_colour_discrete(name="Clusters")



Ngroup <- umapT.layout[which(umapT.layout$Contig==coverage16s[which(coverage16s$covLen==max(coverage16s[which(coverage16s$class==opt$targetName),]$covLen, na.rm = T)),]$Contig),]$group
print(cl.dbscan)
print(paste(opt$targetName, "organism group", sep = " "))
print(Ngroup)

Target <- umapT.layout[which(umapT.layout$group==Ngroup),]
gckOutputTable <- outputTable[which(outputTable$Contig %in% Target$Contig),]
gckOutputTable$cluster <- Target$group

#---------------------------------------
#building graph
g <- read.table(opt$network, stringsAsFactors = F)
g <- g[,c(2,3,1)]
#graph
G <- graph.data.frame(g, directed = F)
if (!is_simple(G)){
  G <- simplify(G, edge.attr.comb = "first")
}
#deleting edges of values lower than N (here N=10)
E<-as.matrix(edge_attr(G)[[1]])
G.copy <- delete.edges(G,which(E<opt$Edges))
#deleting vertices with degree higher than 3 (so the remaining vertices will have degree comprised between 1 to 3)
G.copy <- delete.vertices(G.copy,which(igraph::degree(G.copy)>opt$VerticesDegree))

G.copy <- delete.vertices(G.copy,which(igraph::degree(G.copy)==0))
#summary(G.copy)

C<-components(G.copy)
u<-as.data.frame(unlist(C[["membership"]]))
colnames(u)<-"membership"
unames<-rownames(u)
lcomp<-which(C$csize>=opt$componentSize)
if (order(unique(umapT.layout$group==0), decreasing = T)[1] %in% 2) {
  umapT.layout <- umapT.layout[-which(umapT.layout$group==0),]
}



compTax <- as.data.frame(matrix(data = NA, nrow = length(lcomp), 
                                ncol = (2 + length(unique(umapT.layout$group)))))
colnames(compTax) <- c("comp", "n", unique(umapT.layout$group))


for (l in 1:length(lcomp)) {
  compTax[l,1] <- lcomp[l]
  compTax[l,2] <- length(which(u$membership==l))
  tmp <-which(u$membership==l)
  for (m in 3:ncol(compTax)) {
    compTax[l,m] <- length(which(umapT.layout$Contig %in% rownames(u)[tmp] & 
                                   umapT.layout$group==colnames(compTax)[m]))
  }
  
}

predCCcomp <- data.frame(Contig = rownames(u), 
                         CC = u$membership, stringsAsFactors = F)
predCCcomp <- merge(predCCcomp, compTax, by.x="CC", by.y="comp", all.x = TRUE)
predCCcomp <- merge(predCCcomp, umapT.layout, by = "Contig")

modTaxPredCC <- NULL

for(z in lcomp){
  contigs_in_comp<-predCCcomp[which(predCCcomp$CC==z),]
  l<-length(unique(contigs_in_comp$group))
  
  if (l>1){
    tmp <- as.data.frame(table(contigs_in_comp$group),
                         stringsAsFactors = F)
    tmp$perc <- tmp[,2]/sum(tmp[,2])
    if (nrow(contigs_in_comp)==length(unique(umapT.layout$group))) {
      contigs_in_comp$group <- "misclassified"
      } else if (l==2 && min(tmp$perc)<opt$mixedComponents){
        contigs_in_comp$group <- tmp$Var1[which(tmp$perc==max(tmp$perc))]
        } else if (l>3) {
          contigs_in_comp$group <- "misclassified"
        }
    }
  modTaxPredCC <- rbind(modTaxPredCC, contigs_in_comp)
}


predTable <- merge(umapT.layout[,-which(colnames(umapT.layout)=="group")], modTaxPredCC[,c(1:2, (ncol(modTaxPredCC)-3):(ncol(modTaxPredCC)-1))], by = c("Contig", "V1", "V2"),
                   all.y = TRUE)
umapMod <- umapT.layout[-which(umapT.layout$Contig %in% predTable$Contig),]
umapMod$CC <- rep(-1, nrow(umapMod))

extUmap <- rbind(predTable, umapMod)


homoCC <- 0
compTaxMod <- compTax[-which(rowSums(compTax[,3:ncol(compTax)])==0),]
for (k in 3:ncol(compTaxMod)) {
  for (n in 1:nrow(compTaxMod)) {
    if(compTaxMod[n,k]==sum(compTaxMod[n, 3:ncol(compTaxMod)], na.rm = T) ) {
      homoCC <- homoCC+1
    }
  }
}
h.index <- homoCC/(nrow(compTaxMod))

print("Homogenety index of the Clustering extended with graph based network is:")
print(h.index)



outputContig <- unique(c(extUmap[which(extUmap$group==Ngroup),]$Contig, 
                  unames[which(u$membership %in% unique(extUmap[which(extUmap$group==Ngroup),]$CC))]))
write.table(extUmap, "extendedClusteringCC.txt", sep = "\t", quote = F, row.names = F)
write.table(outputContig, "OutputClustering.txt", sep = "\t", quote = F, row.names = F, col.names = F)

cl <- cl.dbscan[["cluster"]]
clT <- table(cl)

cat("Clustering summary", 
    file = "output_statsClustering.txt", append = T)
cat("\n\n", file =  "output_statsClustering.txt", append = T)
cat("BDscan clustering for objects", file = "output_statsClustering.txt", append = T)
suppressWarnings(write(nrow(outputTable), file = "output_statsClustering.txt", append = T))
cat("\n\n", file = "output_statsClustering.txt", append = T)
cat("the clustering (row 1) contian (row2):\n", file = "output_statsClustering.txt", append = T)
suppressWarnings(write.table(t(clT), "output_statsClustering.txt", 
                             append = T, sep = "\t", quote = F, row.names = F))
cat("\n\n", file = "output_statsClustering.txt", append = T)
cat(paste(opt$targetName, "organism group is ", sep = " "), file = "output_statsClustering.txt", append = T)
suppressWarnings(write(Ngroup, file = "output_statsClustering.txt", append = T))
cat("\n\n", file = "output_statsClustering.txt", append = T)
cat("Homogenety index of the Clustering extended with graph based network is:",  
    file = "output_statsClustering.txt", append = T)
cat("\n\n", file ="output_statsClustering.txt", append = T)
suppressWarnings(write(h.index, "output_statsClustering.txt", 
                             append = T, sep = "\t"))












