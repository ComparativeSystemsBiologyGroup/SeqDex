suppressPackageStartupMessages(library(randomForest, quietly = T))
path <- list.files(path= "~", full.names=TRUE,
                   recursive=TRUE,pattern="(Func.R)")
pathTrue <- path[grep(pattern = "Function/Func.R", path, fixed = T)]
source(pathTrue) #path to Func.R file
suppressPackageStartupMessages(library(optparse, quietly = T))
suppressPackageStartupMessages(library(igraph, quietly = T))
#Options-----------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-G", "--gcCovKmersTable"), action = "store", type = "character", default = "../Coverage/gckCovTable.txt",
              help = "GCKmersCovDef.R output in path/filename.ext - [default %default]"),
  make_option(c("-r", "--rRNA16S"), action = "store", type = "character", default = "../Taxonomy/rRNA16sTaxonomy2.txt",
              help = "rRNA16S.R output table in path/file.ext - [default %default]"),
  make_option(c("-t", "--taxonomy"), action = "store", type = "character", default = NA,
              help = "Taxonomy.R output/s in path/file.ext. If multiple iteration performed,\n\t\t\t then a ordered comma separated list of path/file.ext must be provided \n\t\t [default %default]"),
  make_option(c("-R", "--network"), action = "store", type = "character", default = "NetTax/contigNet.txt",
              help = "reads based graph network path/filename.ext - [default %default]"),
  make_option(c("-l", "--minContigLen"), action = "store", type = "integer", default = 1000, 
              help = "minimum contig length considered - [default %default]"),
  make_option(c("-n", "--targetName"), action = "store", type = "character", default = "Bacteria",
              help = "taxonomy target to be isolated using RF;\n\t\t\t a ordered list of comma separated value for each itariation \n\t\t [default %default]"),
  make_option(c("-N", "--TaxaName"), action = "store", type = "character", default = "superkingdom",
              help = "Taxonomic level on which perfomr RF; \n\t\t\t a ordered list of comma separated levels for each iteration \n\t\t [default %default]"),
  make_option(c("-m", "--nmodel"), action = "store", type = "integer", default = 100,
              help = "Number of SVM model iteratively constructed- [default %default]"),
  make_option(c("-p", "--replace"), action = "store", type = "logical", default = T, 
              help = "Random Forest sampling of cases should be done with or without replacement? [default %default]"),
  make_option(c("-e", "--ntree"), action = "store", type = "integer", default = 500, 
              help = "Numbe of tree of the forest to be grown by Random Forest [default %default]"),
  make_option(c("-T", "--threads"), action = "store", type = "integer", default = 8,
              help = "number of threads to be used - [default %default]"),
  make_option(c("-E", "--Edges"), action = "store", type = "integer", default = 10,
              help = "minimum edges value to be considered to construct network to extend the RF predicted taxonomy - [default %default]"),
  make_option(c("-V", "--VerticesDegree"), action = "store", type = "integer", default = 5,
              help = "maximum vertex degree of the network to be considered to extract component and extend the RF predicted taxonomy - [default %default]"),
  make_option(c("-S", "--componentSize"), action = "store", type = "integer", default = 2,
              help = "minimum size of the components to be considered to extend the RF predicted taxonomy - [default %default]"),
  make_option(c("-M", "--mixedComponents"), action = "store", type = "double", default = 0.2,
              help = "threshold of frequency of alternative RF predicted taxonomy over the total of the component. \n\t\t\t Under this frequency value the alternative will be corrected and substituded by majority - [default %default]"),
  make_option(c("-v", "--verbose"), action = "store", type = "logical", default = T,
              help = "print RF output stats? - [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (grepl(opt$targetName, pattern = ",")) {
  target.name <- as.character(unlist(strsplit(opt$targetName, split = ",", fixed = T)))
} else {
  target.name <- as.character(opt$targetName)
}

if (grepl(opt$TaxaName, pattern = ",")) {
  taxa.name <- tolower(as.character(unlist(strsplit(opt$TaxaName, split = ",", fixed = T))))
} else {
  taxa.name <- tolower(as.character(opt$TaxaName))
}

if(grepl(opt$taxonomy, pattern = "TaxonomyIteration.txt")){
  if (grepl(opt$taxonomy, pattern = ",")) {
    taxa.file <- as.vector(unlist(strsplit(opt$taxonomy, split = ",", fixed = T)))
  } else {
    taxa.file <- as.vector(opt$taxonomy)
  }
} else {
  if (grepl(opt$taxonomy, pattern = ",")) {
    taxa.file <- tolower(as.character(unlist(strsplit(opt$taxonomy, split = ",", fixed = T))))
    taxa.file <- paste("../Taxonomy/", paste(taxa.file, "TaxonomyIteration.txt", sep = ""), sep = "")
  } else {
    taxa.file <- tolower(as.character(opt$taxonomy)) 
    taxa.file <- paste("../Taxonomy/", paste(taxa.file, "TaxonomyIteration.txt", sep = ""), sep = "")
  } 
}
#graph---------------------------------------------------
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

tmpT <- read.table(taxa.file[1], sep = "\t", header = T, stringsAsFactors = F)


#Data------------------------------------------------------------------------
for (j in 1:length(taxa.file)) {
  mergedTaxonomy <- read.table(taxa.file[j], sep = "\t", stringsAsFactors = F, header = T)
  mergedTaxonomy <- mergedTaxonomy[which(mergedTaxonomy$Contig %in% tmpT$Contig),]
  gckCovTable <- read.table(opt$gcCovKmersTable, sep = "\t", stringsAsFactors = F, header = T)
  gckCovTable <- gckCovTable[which(gckCovTable$length>opt$minContigLen),]
  
  rRNA16sTaxonomy2 <- read.table(opt$rRNA16S, sep = "\t", stringsAsFactors = F, header = T)
  rRNA16sTaxonomy2 <- rRNA16sTaxonomy2[which(rRNA16sTaxonomy2$ali_Length>100 ),]
  rRNA16sTaxonomy3 <- data.frame(Contig = rRNA16sTaxonomy2[,2],
                                 TaxonDensity = 1,
                                 taxa.name = rRNA16sTaxonomy2[,which(tolower(colnames(rRNA16sTaxonomy2))==taxa.name[j])], stringsAsFactors = F)
  colnames(rRNA16sTaxonomy3)[3] <- taxa.name[j]
  mergedTaxonomy <- rbind.data.frame(mergedTaxonomy,
                                     rRNA16sTaxonomy3, stringsAsFactors = F)
  mergedTaxonomy <- unique(mergedTaxonomy)
  
  freq <- as.data.frame(table(mergedTaxonomy[,which(tolower(colnames(mergedTaxonomy))==taxa.name[j])]))
  x <- freq[which(freq$Freq>5),]
  mergedTaxonomy <- mergedTaxonomy[which(mergedTaxonomy[,which(tolower(colnames(mergedTaxonomy))==taxa.name[j])] %in% x$Var1),]
  
  completeTable <- merge(x=gckCovTable, y=mergedTaxonomy, #then merge the two tables
                         by= "Contig", all.x =TRUE) #merge put NA's where there is not the value
  completeTable[is.na(completeTable)] <-"NoBlastHit" #Substitute NA's with another variabe name
  completeTableF <- completeTable[which(completeTable$length>=opt$minContigLen),]
  
  completeTableNum <- taxonNumber(completeTableF, taxa.name[j], "NoBlastHit", 1)
  taxaColumn = which(tolower(colnames(completeTableF))==taxa.name[j])
  name <- unique(completeTable[-which(completeTable[,taxaColumn]=="NoBlastHit"),taxaColumn]) 
  completeTableNB <- completeTable[which(completeTable[,taxaColumn]=="NoBlastHit"),]
  BacNum <- which(tolower(unique(completeTable[-which(completeTable[,taxaColumn]=="NoBlastHit"), taxaColumn]))==tolower(target.name[j]))
  print(paste(target.name[j], "name is converted to number:", sep = " "))
  print(BacNum)
  
  TAXA = which(tolower(colnames(completeTableNum))==taxa.name[j])
  
  
  COLUMN = which(colnames(completeTableNum) %in% colnames(completeTableNum)[-which(colnames(completeTableNum)=="Contig" |
                                                                                     colnames(completeTableNum)=="GCcont" |
                                                                                     colnames(completeTableNum)=="covLen" |
                                                                                     colnames(completeTableNum)=="length" |
                                                                                     colnames(completeTableNum)=="TaxonDensity" |
                                                                                     tolower(colnames(completeTableNum))==taxa.name[j])])
  
  
  RF.model <- RFModelParIter(table = completeTableNum, column = COLUMN, 
                             taxa = TAXA, cluster = opt$threads, replace = opt$replace, ntree = opt$ntree,
                            nmodel = opt$nmodel,  verbose = opt$verbose) 
  
  outptStats_N <- data.frame(matrix(data = NA, nrow = 4, ncol = (length(RF.model)+1)), 
                             stringsAsFactors = F)
  outptStats_L <- data.frame(matrix(data = NA, nrow = 4, ncol = (length(RF.model)+1)), 
                             stringsAsFactors = F)
  
  
  for (i in 1:(length(RF.model)-1)) {
    pred.stats <- data.frame(true = as.vector(RF.model[[i]][2]),
                             pred = as.vector(RF.model[[i]][3]),
                             length = as.vector(RF.model[[i]][4]), stringsAsFactors = F)
    tmpN <- outputStats(pred.stats, real_column = 1, real_name = BacNum, 
                        predicted_colum = 2, predicted_name = BacNum)[c(5,9,10,11),]
    tmpL <- outputStatsLen(pred.stats, real_column = 1, real_name = BacNum, 
                           predicted_colum = 2, predicted_name = BacNum, length_column = 3)[c(5,9,10,11),]
    if (i==1) {
      outptStats_N[,1:3] <- tmpN
      outptStats_L[,1:3] <- tmpL
    } else {
      outptStats_N[,(i+2)] <- tmpN[,3]
      outptStats_L[,(i+2)] <- tmpL[,3]
    }
    
  } 
  outptStats_N$X2 <- rep(name[as.numeric(unique(outptStats_N$X2))], nrow(outptStats_N))
  outptStats_L$X2 <- rep(name[as.numeric(unique(outptStats_L$X2))], nrow(outptStats_L))
  
  statisticsN <- as.data.frame(matrix(NA, nrow = 4, ncol = 6), stringsAsFactors = F)
  colnames(statisticsN) <- c("stats", "test", "N","mean", "sd", "se")
  statisticsN$stats <- outptStats_N$X1
  statisticsN$test <- outptStats_N$X2
  statisticsL <- statisticsN
  
  
  for (n in 1:nrow(outptStats_N)) {
    x <- outptStats_N[n,-c(1, 2)]
    
    statisticsN[n,3] <- length(t(x))
    statisticsN[n,4] <- mean(t(x))
    statisticsN[n,5] <- sd(t(x))
    statisticsN[n,6] <- statisticsN[n,5]/length(t(x))
    
    
    y <- outptStats_L[n,-c(1, 2)]
    
    statisticsL[n,3] <- length(t(y))
    statisticsL[n,4] <- mean(t(y))
    statisticsL[n,5] <- sd(t(y))
    statisticsL[n,6] <- statisticsN[n,5]/length(t(y))
    
  }
  
  error.predN <- as.data.frame(RF.model$error[1], stringsAsFactors = F)
  col <- gsub("X", "", colnames(error.predN))
  col <- col[-length(col)]
  col <- as.numeric(col)
  row.names(error.predN) <- name
  colnames(error.predN)[-ncol(error.predN)] <- name[col]
  
  error.predL <- as.data.frame(RF.model$error[2], stringsAsFactors = F)
  row.names(error.predL) <- name
  colnames(error.predL)[-ncol(error.predL)] <- name[col]
  
  cat("percentage of error committed in reclassification, considering number of sequences", 
      file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  suppressWarnings(write.table(error.predN, paste(taxa.name[j],"Output_statsRF.txt", sep = ""), 
                               append = T, sep = "\t", quote = F, row.names = T ))
  cat("\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("percentage of error committed in reclassification, considering length of the sequences",  
      file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  suppressWarnings(write.table(error.predL, paste(taxa.name[j],"Output_statsRF.txt", sep = ""), 
                               append = T, sep = "\t", quote = F, row.names = T ))
  cat("\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("performance statistics of RF model calculated on number of target sequences correctly reclassified",
      file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  suppressWarnings(write.table(statisticsN, paste(taxa.name[j],"Output_statsRF.txt", sep = ""), 
                               append = T, sep = "\t", quote = F, row.names = F ))
  cat("\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("performance statistics of RF model calculated on cumulative amount of target bp correctly reclassified", 
      file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
  suppressWarnings(write.table(statisticsL, paste(taxa.name[j],"Output_statsRF.txt", sep = ""), 
                               append = T, sep = "\t", quote = F, row.names = F ))
  
  pred <- predictRFIter(RF.model, completeTableNB, COLUMN)
  
  pred.filt <- as.data.frame(matrix(data=NA, nrow = 0, ncol = 3), stringsAsFactors = F)
  for (k in 1:length(unique(pred$V1))) {
    tmp <- pred[which(pred$V1==unique(pred$V1)[k]),]
    tmp <- tmp[which(tmp$V3 == max(tmp$V3)),]
    
    pred.filt <- rbind.data.frame(pred.filt, tmp, stringsAsFactors = F)
  }
  colnames(pred.filt) <- c("Contig", taxa.name[j], "percent")
  
  
  pTAXA = which(tolower(colnames(pred.filt))==taxa.name[j])
  compTax <- as.data.frame(matrix(data = NA, nrow = length(lcomp), 
                                  ncol = (2 + length(unique(pred.filt[,pTAXA])))))
  colnames(compTax) <- c("comp", "n", unique(pred.filt[,pTAXA]))
  
  
  for (l in 1:length(lcomp)) {
    compTax[l,1] <- lcomp[l]
    compTax[l,2] <- length(which(u$membership==l))
    tmp <-which(u$membership==l)
    for (m in 3:ncol(compTax)) {
      compTax[l,m] <- length(which(pred.filt$Contig%in%rownames(u)[tmp] & 
                                     pred.filt[,pTAXA]==colnames(compTax)[m]))
    }
    
  }
  
  predCCcomp <- data.frame(Contig = rownames(u), 
                           CC = u$membership, stringsAsFactors = F)
  predCCcomp <- merge(predCCcomp, compTax, by.x="CC", by.y="comp", all.x = TRUE)
  predCCcomp <- merge(predCCcomp, pred.filt, by = "Contig")
  
  modTaxPredCC <- NULL
  
  for(z in lcomp){
    contigs_in_comp<-predCCcomp[which(predCCcomp$CC==z),]
    l<-length(unique(contigs_in_comp[,colnames(contigs_in_comp)==taxa.name[j]]))
    
    if (l>1){
      tmp <- as.data.frame(table(contigs_in_comp[,colnames(contigs_in_comp)==taxa.name[j]]),
                           stringsAsFactors = F)
      tmp$perc <- tmp[,2]/sum(tmp[,2])
      if (nrow(contigs_in_comp)==length(unique(completeTableNum[,TAXA]))) {
        contigs_in_comp[,colnames(contigs_in_comp)==taxa.name[j]] <- "misclassified"
      } else if (l==2 && min(tmp$perc) <opt$mixedComponents){
        contigs_in_comp[,colnames(contigs_in_comp)==taxa.name[j]] <- tmp$Var1[which(tmp$perc==max(tmp$perc))]
      } else if (l>3) {
        contigs_in_comp[,colnames(contigs_in_comp)==taxa.name[j]] <- "misclassified"
      }
      
      
    }
    modTaxPredCC <- rbind(modTaxPredCC, contigs_in_comp)
  }
  
  predTable <- merge(completeTableNB[,-ncol(completeTableNB)], pred.filt, by = "Contig")
  tableNB <- merge(predTable[,1:(ncol(predTable)-2)], modTaxPredCC[,c(1,(ncol(modTaxPredCC)-1):ncol(modTaxPredCC))],
                   by="Contig", all.y = TRUE)
  tableNB <- rbind(tableNB, predTable[-which(predTable$Contig %in% tableNB$Contig),])
  for (k in 1:length(unique(tableNB$Contig))) {
    tmp <- tableNB[which(tableNB$Contig==unique(tableNB$Contig)[k]),]
    if (nrow(tmp)==2) {
      tableNB <- tableNB[-which(tableNB$Contig==unique(tableNB$Contig)[k])[2],]
      tableNB[k,c(TAXA, (TAXA+1))] <- c("misclassified", "misclassified")
    }
  }
  
  completeTableNum$percent <- rep(-1, nrow(completeTableNum))
  output <- rbind.data.frame(tableNB[which(tableNB[,taxaColumn]==BacNum),],
                             tableNB[which(tableNB[,taxaColumn]=="misclassified"),],
                             completeTableNum[which(completeTableNum[,taxaColumn]==BacNum),],
                             stringsAsFactors = F)
  
  write.table(output, paste(taxa.name[j],"OutputRF.txt", sep = ""), sep = "\t",quote = F, row.names = F)
  
  homoCC <- 0
  compTaxMod <- compTax[-which(compTax[,3:ncol(compTax)]==0),]
  if(nrow(compTaxMod)>0){
    for (k in 3:ncol(compTaxMod)) {
      for (n in 1:nrow(compTaxMod)) {
        if(compTaxMod[n,k]==sum(compTaxMod[n, 3:ncol(compTaxMod)])) {
          homoCC <- homoCC+1
        }
      }
    }
    h.index <- homoCC/(nrow(compTaxMod))
    print("Homogenety index of the RF prediction extended with graph based network is:")
    print(h.index)
    
    cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
    cat("Homogenety index for RF machine learning predictions. \n\nThis value have been calculated by dividing number of non zero homogeneous components to the total number of components", 
        file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T )
    cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
    suppressWarnings(write(h.index, file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""),
                           append = T))
  } else {
    h.index <- NA
    print("Homogenety index of the RF prediction extended with graph based network is:")
    print("no component found")
    
    cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
    cat("Homogenety index for RF machine learning predictions. \n\nThis value have been calculated by dividing number of non zero homogeneous components to the total number of components", 
        file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T )
    cat("\n\n", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), append = T)
    cat("No component found", file = paste(taxa.name[j],"Output_statsRF.txt", sep = ""), 
        append = T)
  }
  
  save.image(paste(taxa.name[j],"RFModel.RData", sep = ""))
  
  tmpT <- output[,-ncol(output)]
}






