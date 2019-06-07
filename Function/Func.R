#Best taxa choice-------------------------------------------------------------------------------

BestTaxa <- function(table) {
  table[is.na(table)] <- "TaxaNoDef"
  query.uniq <- unique(table$Contig)
  table2 <- as.data.frame(matrix(data = 0, nrow = 0, ncol = ncol(table)))
  colnames(table2) <- colnames(table)
  for (n in 1:length(query.uniq)) {
    Idx <- unique(table[which(table$Contig==query.uniq[n]),])
    if (nrow(Idx)>1) {
      if (length(unique(Idx$order))>1){
        #print(Idx$Contig)
        a <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(table)))
        colnames(a) <- colnames(table)
          a <- Idx[which(Idx$X.id==max(Idx$X.id)),][1,]
        #for (k in 1:length(unique(Idx$order))) {
         # Idx2 <- Idx[which(Idx$order==unique(Idx$order)[k]),]
          #}
        table2 <- rbind(table2, a)
        } else {
          Idx2 <- Idx[which(Idx$order==unique(Idx$order)),]
          a <- as.data.frame(matrix(data = NA, nrow = 1, ncol = ncol(table)))
          colnames(a) <- colnames(table)
          a <- Idx2[which(Idx2$X.id==max(Idx2$X.id)),][1,]
          table2 <- rbind(table2, a)
          }
      } else {
        table2 <- rbind(table2, Idx)
      }
    }
  return(table2)
}

#all taxonomic affiliation parallelized---------------------------------------------------------
BlastTaxaAll <- function(table, blast_input, taxa){
  
  if (blast_input=="nt"){
    subject <- strsplit(table$subject, "|", fixed = TRUE)
    db <- unlist(lapply(subject, "[[", 4)) 
    table$subject <- db
  } else if (blast_input=="diamond") {
    name <- strsplit(table$query, "|", fixed = TRUE)
    name2 <- unlist(lapply(name, "[", 2))
    name3 <- strsplit(name2, ":", fixed = TRUE)
    name4 <- unlist(lapply(name3, "[", 1))
    name5 <- strsplit(name4, "_", fixed = TRUE)
    name6 <- paste(unlist(lapply(name5, "[", 2)), 
                   unlist(lapply(name5, "[", 3)),
                   unlist(lapply(name5, "[", 4)),
                   unlist(lapply(name5, "[", 5)),
                   unlist(lapply(name5, "[", 6)),
                   unlist(lapply(name5, "[", 7)),
                   sep = "_"
    )
    table$query <- name6
  } else if (blast_input=="chromosome") {
    NULL
  } else {
    print("error1-blast_input bad names called")
  }
  
  Blast <- subset.data.frame(table,
                             select = c("query", "subject", "X.id","ali_Length", "QLEN", "QCOVHSP"))
  #TX<-matrix(data=0,nrow=nrow(Blast),ncol=length(taxa)) #empty matrix for the loop
  U <- unique(Blast$subject) #vector with unique subject in Blast
  Blast2 <- as.data.frame(matrix(data=0, 
                                 nrow = nrow(Blast), 
                                 ncol = (ncol(Blast) + length(taxa))
  ))
  Blast2[,1:ncol(Blast)] <- Blast
  path <- list.files(path= "~", full.names=TRUE,
                     recursive=TRUE,pattern="(accessionTaxa.sql)")
  taxaId<-taxonomizr::accessionToTaxa(Blast2$V2, path) #retiving accession number
  X=taxonomizr::getTaxonomy(taxaId, path) #get taxonomic affiliation
  X <- X[,taxa]
  Blast2[,(ncol(Blast)+1):(ncol(Blast2))] <- X
  colnames(Blast2) <- c(colnames(Blast), colnames(X))
  return(Blast2)
  
}

#obtaining taxonomic density for Blast Table with accession number already converted------------
TaxonomyDensity <- function(table, taxaColumn) {
  #table[is.na(table)] <- "TaxaNoDef"
  query.uniq <- unique(table$Contig)
  table2 <- as.data.frame(matrix(data = 0, nrow = 0, ncol = 3))
  colnames(table2) <- c("Contig", "TaxonDensity", colnames(table)[taxaColumn] )
  for (n in 1:length(query.uniq)) {
    Idx <- unique(table[which(table$Contig==query.uniq[n]),])
    nTaxa <- as.data.frame(table(Idx[,taxaColumn]), stringsAsFactors = FALSE)
    nTaxa$Percent <- nTaxa[,2]/sum(nTaxa[,2])
    table3 <- as.data.frame(matrix(data = 0, nrow = nrow(nTaxa), ncol = 3),
                            stringsAsFactors = FALSE)
    colnames(table3) <- c("Contig", "TaxonDensity", colnames(table)[taxaColumn] )
    for (i in 1:nrow(nTaxa)) {
      table3[i, 1] <- query.uniq[n]
      table3[i, 2] <- nTaxa[i,3]
      table3[i, 3] <- nTaxa[i, 1]
    }
    if (max(table3$TaxonDensity)>0.3) {
      table4 <- table3[which(table3$TaxonDensity>=0.001),]
      table2 <- rbind(table2, table4)
    } else {
      NULL
    }
    
  }
  return(table2)
}



##Kmers by length-------------------------------------------------------------------------------
KmersCountNorm <- function(KMER_SIZE, contigs, cluster ){
  gcx<-as.data.frame(matrix(data=0,nrow=length(contigs),ncol=2)) #table for the loop
  colnames(gcx) <- c("Contig", "GCcont")
  Kmers <- as.data.frame(matrix(data = 0, nrow = length(contigs), ncol = 4^KMER_SIZE))
  cl = parallel::makeCluster(cluster)
  parallel::clusterExport(cl, "contigs")
  Kmers <- t(parallel::parSapply(cl, contigs, seqinr::count, KMER_SIZE))
  #kmers tabel have to be trasposed due to exit order
  #retriving polymer order to be used as colnames in kmers table
  #setting colnames of gcx
  gcx[ ,1] <- rownames(Kmers)
  gcx[, 2] <- parallel::parSapply(cl, contigs, seqinr::GC)
  parallel::stopCluster(cl)
  #counting reverse complement to sum. We avoid counting diretly the reverse complement to avoid
  ##time consuming script. So each polynucleotide have to be summed 
  ###to its complement and reverse complement.
  ####Palindrome have to be summed only to itself (they are already counted twice).
  #####So the for loop is divided in tosum= not palindrome, and sametosum=palindrome
  tosum <- matrix(data = NA, nrow = 2, ncol = 4^KMER_SIZE) #with NA to be able to remove efficently
  sametosum <- matrix(data = NA, nrow = 2, ncol = 4^KMER_SIZE)
  x <- c()
  for(i in 1:4^KMER_SIZE){
    x[i] <-which(paste(unlist(rev(seqinr::comp(strsplit(colnames(Kmers), "")[[i]]))), collapse ="") == colnames(Kmers))
    if(x[i]<i){
      i -> tosum[1, i]
      x[i]->tosum[2, i]
    } else if(x[i]==i){
      i -> sametosum[1, i]
      x[i]->sametosum[2, i]
    } else {
      NULL
    }
    
  }
  tosum<-tosum[,!colSums(!is.finite(tosum))] #removing NA
  sametosum<-sametosum[,!colSums(!is.finite(sametosum))]
  if (length(sametosum)==0) {
    tosumKmers <- as.data.frame(matrix(data=0,nrow=length(contigs),ncol=(ncol(tosum))))
    colnames(tosumKmers) <-colnames(Kmers[,tosum[2,]]) #creating table with plynucleotide to be summed
    for (i in 1:(ncol(tosum))) {
      tosumKmers[,i] <- (2*Kmers[,tosum[1,i]]) + (2*Kmers[, tosum[2,i]]) 
    }
    gcKmers <- cbind(gcx, tosumKmers) #putting table all together
    gcKmers$Num <- rowSums(gcKmers[, 3:ncol(gcKmers)])
    
  } else if (length(sametosum)!=0) {
    tosumKmers <- as.data.frame(matrix(data=0,nrow=length(contigs),ncol=(ncol(tosum))))
    colnames(tosumKmers) <-colnames(Kmers[,tosum[2,]]) #creating table with plynucleotide to be summed
    for (i in 1:(ncol(tosum))) {
      tosumKmers[,i] <- (2*Kmers[,tosum[1,i]]) + (2*Kmers[, tosum[2,i]]) 
    }
    sameKmers <- as.data.frame(matrix(data=0,nrow=length(contigs),ncol=(ncol(sametosum))))
    colnames(sameKmers) <-colnames(Kmers[,sametosum[2,]])
    for (i in 1:(ncol(sametosum))) {
      sameKmers[,i] <- Kmers[,sametosum[1,i]] + Kmers[, sametosum[2,i]] 
    }
    
    
    gcKmers <- cbind(gcx, tosumKmers, sameKmers) #putting table all together
    gcKmers$Num <- rowSums(gcKmers[, 3:ncol(gcKmers)])
  }
  
  ###Kmer count have to be standardized by all the count 
  ###(have to be similar to the contig length, but it is not the same)
  
  gcxlen<- data.frame(gcKmers[,1:2], gcKmers[,3:(ncol(gcKmers)-1)]/gcKmers[,ncol(gcKmers)])
  #toc()
  return(gcxlen)
}
#converting taxonomy to number------------------------------------------------------------------
taxonNumber <- function(table, Taxa, noblast, n){
  
  x<- which(colnames(table) == Taxa)
  NoBlastHit <- which(table[,x] == noblast)
  taxonomy <- matrix(data=0, nrow = nrow(table), ncol = 1)
  uorder<-unique(table[,x])
  print(uorder)
  taxonomy[NoBlastHit]<-NA
  if (length(NoBlastHit)!=0){
    uorder<-uorder[-which(uorder==noblast)]
  } else if (length(NoBlastHit)==0){
    NULL
  }
  
  for(i in 1:length(uorder)){
    idx<-which(table[,x]==uorder[i])
    taxonomy[idx]<-i
  }
  table[,x] <- taxonomy
  if (n==1){
    table <- table[-NoBlastHit,]
  } else if (n==2){
    NULL
  }
  return(table)
}

#SVM--------------------------------------------------------------------------------------------
SvmModelParIter <- function(table,column, taxa, cost, gamma, scale, cross, cluster, verbose = T ){
  # tic()
  index <- 1:nrow(table)
  testindex <- sample(index, trunc(length(index)/3))
  testset <- table[testindex,]
  
  trainset <- table[-testindex,]
  
  
  grid <- data.frame(gamma = expand.grid(gamma, cost)[,1],
                     cost = expand.grid(gamma, cost)[,2],
                     stringsAsFactors = F)
  suppressPackageStartupMessages(library(doParallel, quietly = T))
  cl = parallel::makeCluster(cluster)
  doParallel::registerDoParallel(cl)
  model.tune <- foreach::foreach(i=1:nrow(grid)) %dopar% {
    e1071::tune.svm(trainset[,column],
                    as.factor(trainset[,taxa]),
                    cost = grid$cost[i], gamma = grid$gamma[i], scale = scale)
    
  }
  
  param.tune <- data.frame(t(sapply(model.tune, "[[", 7 )), stringsAsFactors = F )
  param.tune2 <- data.frame(gamma = unlist(param.tune$gamma),
                            cost = unlist(param.tune$cost),
                            error = unlist(param.tune$error),
                            dispersion = unlist(param.tune$dispersion),
                            stringsAsFactors = F)
  
  svm.model <- foreach::foreach(n=1:100) %dopar% {
    index <- 1:nrow(table)
    testindex <- sample(index, trunc(length(index)/3))
    testset <- table[testindex,]
    
    trainset <- table[-testindex,]
    x <- e1071::svm( as.factor(trainset[,taxa]) ~., data=trainset[,column],
                     type = "C-classification",
                     cost = param.tune2[which(param.tune2$error==min(param.tune2$error)),which(colnames(param.tune2)=="cost")], 
                     gamma = param.tune2[which(param.tune2$error==min(param.tune2$error)),which(colnames(param.tune2)=="gamma")], 
                     cross=cross, scale = scale)
    
    y <- as.character(predict(x, testset[,column]))
    z <- testset[,taxa]
    l <- testset[, which(colnames(testset)=="length")]
    list(x,true = z, pred =y, length = l)
  }
  #svm.pred <- lapply(svm.model, "[[", c(2,3))
  
  
  pred <- data.frame(matrix(data = 0, nrow = 0, ncol = 0), stringsAsFactors = F)
  
  parallel::stopCluster(cl)
  
  for (k in 1:100) {
    tmp <- data.frame(true = svm.model[[k]][2], pred = svm.model[[k]][3], length = svm.model[[k]][4],
                      stringsAsFactors = F)
    #tmp <- cbind.data.frame(testset[,c(1, taxa)], tmp)
    pred <- rbind.data.frame(pred, tmp, stringsAsFactors = F)
  }
  #print(table(pred=pred$tmp, true = pred[,which(colnames(pred)==colnames(testset)[taxa])]))
  freq.pred <- as.data.frame(table(true = pred[,1], pred=pred[,2]) ,
                             stringsAsFactors = F)
  
  error.pred <- data.frame(matrix(data=0, nrow = length(unique(freq.pred$true)), ncol = length(unique(freq.pred$pred))))
  rownames(error.pred) <- sort(unique(freq.pred$true))
  colnames(error.pred) <- sort(unique(freq.pred$pred))
  
  for (n in 1:ncol(error.pred)) {
    for (l in 1:nrow(error.pred)) {
      
        error.pred[l,n] <- freq.pred[which(freq.pred$true==rownames(error.pred)[l] & 
                                             freq.pred$pred ==colnames(error.pred)[n]),3] 
      
    }
  }
  
  error.pred$percent.error <- rep(NA, nrow(error.pred))
  for (j in 1:nrow(error.pred)){
    error.pred[j,ncol(error.pred)] <- sum(error.pred[j, which(colnames(error.pred)[-which(colnames(error.pred)=="percent.error")]!=row.names(error.pred)[j])]/sum((error.pred[j,]), na.rm = T))*100
  }
  
  error.length <- data.frame(matrix(data = 0, nrow = length(unique(pred$true)), 
                                    ncol = length(unique(pred$pred))), stringsAsFactors = F)
  colnames(error.length) <- sort(unique(pred$pred))
  row.names(error.length) <- sort(unique(pred$true))
  for (r in 1:nrow(error.length)) {
    for (c in 1:ncol(error.length)) {
      error.length[r,c] <- sum(pred[which(pred$true==row.names(error.length)[r] & 
                                            pred$pred==colnames(error.length)[c]),]$length)
    }
  }
  
  error.length$percent.error <- rep(NA, nrow(error.length))
  for (j in 1:nrow(error.length)){
    error.length[j,ncol(error.length)] <- sum(error.length[j, which(colnames(error.length)[-which(colnames(error.length)=="percent.error")]!=row.names(error.length)[j])]/sum(error.length[j,], na.rm = T))*100
  }
  
  if (verbose==T) {
    print("percentage of error committed in reclassification, considering number of sequences")
    print(error.pred)
    print("percentage of error committed in reclassification, considering length of sequences")
    print(error.length)
  }
  svm.model[["error"]] <- list(error.pred, error.length)
  
  return(svm.model)
}

predictSVMIter <- function(model, table, column){
  if (attr(model[[1]][3], which = "name")!="pred") {
    stop("ERROR: SVM model needed is from SvmModelParIter function")
  } else {
    NULL
  }
  
  pred <- as.data.frame(matrix(data = NA, nrow = nrow(table), 
                               ncol = (length(model))), stringsAsFactors = F)
  pred[,1] <- table[,which(colnames(table)=="Contig")]
  
  for (i in 1:(length(model)-1)) {
    pred[,(i+1)] <-as.vector(predict(model[[i]][[1]], table[,column]))
  }
  
  outpt <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  for (l in 1:nrow(pred)) {
    unq <- unique(t(pred[l, 2:ncol(pred)]))
    tmp <- as.data.frame(matrix(data = NA, nrow = length(unq), ncol = 3))
    tmp[,1] <- rep(pred[l,1], length(unq))
    for (j in 1:length(unq)) {
      tmp[j,2] <- unq[j]
      tmp[j,3] <- length(which(pred[l,]==unq[j]))/(ncol(pred)-1)
      
    }
    outpt <- rbind.data.frame(outpt, tmp, stringsAsFactors = F)
  }
  return(outpt)
}

#Random Forset function------------------------------------------------------------------------
RFModelParIter <- function(table,column, taxa, cluster, replace, ntree, verbose = T ){
  suppressPackageStartupMessages(library(doParallel, quietly = T))
  cl = parallel::makeCluster(cluster)
  doParallel::registerDoParallel(cl)
  
  RF.model <- foreach::foreach(n=1:100) %dopar% {
    index <- 1:nrow(table)
    testindex <- sample(index, trunc(length(index)/3))
    testset <- table[testindex,]
    
    trainset <- table[-testindex,]
    x <- randomForest::randomForest(as.factor(trainset[,taxa]) ~., 
                                    data=trainset[,column],
                                    ntree = ntree,
                                    replace= replace
    )
    
    y <- as.character(predict(x, testset[,column]))
    z <- testset[,taxa]
    l <- testset[, which(colnames(testset)=="length")]
    list(x,true = z, pred =y, length = l)
  }
  
  
  pred <- data.frame(matrix(data = 0, nrow = 0, ncol = 0), stringsAsFactors = F)
  
  parallel::stopCluster(cl)
  
  for (k in 1:100) {
    tmp <- data.frame( true = RF.model[[k]][2], pred = RF.model[[k]][3], length=RF.model[[k]][4],stringsAsFactors = F)
    #tmp <- cbind.data.frame(testset[,c(1, taxa)], tmp)
    pred <- rbind.data.frame(pred, tmp, stringsAsFactors = F)
  }
  #print(table(pred=pred$tmp, true = pred[,which(colnames(pred)==colnames(testset)[taxa])]))
  freq.pred <- as.data.frame(table(true = pred[,1], pred=pred[,2]) ,
                             stringsAsFactors = F)
  
  error.pred <- data.frame(matrix(data=0, nrow = length(unique(freq.pred$true)), ncol = length(unique(freq.pred$pred))))
  rownames(error.pred) <- sort(unique(freq.pred$true))
  colnames(error.pred) <- sort(unique(freq.pred$pred))
  
  for (n in 1:ncol(error.pred)) {
    for (l in 1:nrow(error.pred)) {
      error.pred[l,n] <- freq.pred[which(freq.pred$true==rownames(error.pred)[l] & 
                                           freq.pred$pred ==colnames(error.pred)[n]),3]
    }
  }
  
  error.pred$percent.error <- rep(NA, nrow(error.pred))
  for (j in 1:nrow(error.pred)){
    error.pred[j,ncol(error.pred)] <- sum(error.pred[j, which(colnames(error.pred)[-which(colnames(error.pred)=="percent.error")]!=row.names(error.pred)[j])]/sum(error.pred[j,], na.rm = T))*100
  }
  
  error.length <- data.frame(matrix(data = 0, nrow = length(unique(pred$true)), 
                                    ncol = length(unique(pred$pred))), stringsAsFactors = F)
  colnames(error.length) <- sort(unique(pred$pred))
  row.names(error.length) <- sort(unique(pred$true))
  for (r in 1:nrow(error.length)) {
    for (c in 1:ncol(error.length)) {
      error.length[r,c] <- sum(pred[which(pred$true==row.names(error.length)[r] & 
                                            pred$pred==colnames(error.length)[c]),]$length)
    }
  }
  
  error.length$percent.error <- rep(NA, nrow(error.length))
  for (j in 1:nrow(error.length)){
    error.length[j,ncol(error.length)] <- sum(error.length[j, which(colnames(error.length)[-which(colnames(error.length)=="percent.error")]!=row.names(error.length)[j])]/sum(error.length[j,], na.rm = T))*100
  }
  
  
  
  if (verbose==T) {
    print("percentage of error committed in reclassification, considering number of sequences")
    print(error.pred)
    print("percentage of error committed in reclassification, considering length of sequences")
    print(error.length)
  }
  RF.model[["error"]] <- list(error.pred, error.length)
  
  return(RF.model)
}
predictRFIter <- function(model, table, column){
  if (attr(model[[1]][3], which = "name")!="pred") {
    stop("ERROR: RF model needed is from SvmModelParIter function")
  } else {
    NULL
  }
  
  pred <- as.data.frame(matrix(data = NA, nrow = nrow(table), 
                               ncol = (length(model))), stringsAsFactors = F)
  pred[,1] <- table[,which(colnames(table)=="Contig")]
  
  for (i in 1:(length(model)-1)) {
    pred[,(i+1)] <-as.vector(predict(model[[i]][[1]], table[,column]))
  }
  
  outpt <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
  for (l in 1:nrow(pred)) {
    unq <- unique(t(pred[l, 2:ncol(pred)]))
    tmp <- as.data.frame(matrix(data = NA, nrow = length(unq), ncol = 3))
    tmp[,1] <- rep(pred[l,1], length(unq))
    for (j in 1:length(unq)) {
      tmp[j,2] <- unq[j]
      tmp[j,3] <- length(which(pred[l,]==unq[j]))/(ncol(pred)-1)
      
    }
    outpt <- rbind.data.frame(outpt, tmp, stringsAsFactors = F)
  }
  return(outpt)
}

#OUTPUT STATS Ncontig---------------------------------------------------------------------------
outputStats <- function(table, real_column, real_name, predicted_colum, predicted_name){
  output <- as.data.frame(matrix(data = NA, nrow = 14, ncol = 3))
  output[,1] <- c("TP", "FP", "FN", "TN", "Sensitivity_TPR", 
                      "FNR", "FPR", "Specificity_TNR", "precision", 
                      "Accuracy", "F1", "PLR", "NLR", "ODDSR")
  
  colnames(output) <- c("stats", "test", "value")
  output$test <- rep(predicted_name, nrow(output))
  
  output[1,3] <- nrow(table[which(table[,real_column]==real_name & table[,predicted_colum]==predicted_name),])
  output[2,3] <- nrow(table[which(table[,real_column]!=real_name & table[,predicted_colum]==predicted_name),])
  output[3,3] <- nrow(table[which(table[,real_column]==real_name & table[,predicted_colum]!=predicted_name),])
  output[4,3] <- nrow(table[which(table[,real_column]!=real_name & table[,predicted_colum]!=predicted_name),])
  output[5,3] <- output[1,3]/(output[1,3]+output[3,3])
  output[6,3] <- output[3,3]/(output[1,3]+output[3,3])
  output[7,3] <- output[2,3]/(output[2,3]+output[4,3])
  output[8,3] <- output[4,3]/(output[2,3]+output[4,3])
  output[9,3] <- output[1,3]/(output[1,3]+output[2,3])
  output[10,3] <- (output[1,3]+output[4,3])/(output[1,3]+output[2,3]+output[3,3]+output[4,3])
  output[11,3] <- (2*output[1,3])/((2*output[1,3])+output[2,3]+output[3,3])
  output[12,3] <- output[5,3]/output[7,3]
  output[13,3] <- output[6,3]/output[8,3]
  output[14,3] <- output[12,3]/output[13,3]
  
  return(output)
}

#OUTPUT STATS Lcontig---------------------------------------------------------------------------

outputStatsLen <- function(table, real_column, real_name, predicted_colum, predicted_name, length_column){
  output <- as.data.frame(matrix(data = NA, nrow = 14, ncol = 3))
  output[,1] <- c("TP", "FP", "FN", "TN", "Sensitivity_TPR", 
                  "FNR", "FPR", "Specificity_TNR", "precision", 
                  "Accuracy", "F1", "PLR", "NLR", "ODDSR")
  
  colnames(output) <- c("stats", "test", "value")
  output$test <- rep(predicted_name, nrow(output))
  
  output[1,3] <- sum(table[which(table[,real_column]==real_name & table[,predicted_colum]==predicted_name),length_column])
  output[2,3] <- sum(table[which(table[,real_column]!=real_name & table[,predicted_colum]==predicted_name),length_column])
  output[3,3] <- sum(table[which(table[,real_column]==real_name & table[,predicted_colum]!=predicted_name),length_column])
  output[4,3] <- sum(table[which(table[,real_column]!=real_name & table[,predicted_colum]!=predicted_name),length_column])
  output[5,3] <- output[1,3]/(output[1,3]+output[3,3])
  output[6,3] <- output[3,3]/(output[1,3]+output[3,3])
  output[7,3] <- output[2,3]/(output[2,3]+output[4,3])
  output[8,3] <- output[4,3]/(output[2,3]+output[4,3])
  output[9,3] <- output[1,3]/(output[1,3]+output[2,3])
  output[10,3] <- (output[1,3]+output[4,3])/(output[1,3]+output[2,3]+output[3,3]+output[4,3])
  output[11,3] <- (2*output[1,3])/((2*output[1,3])+output[2,3]+output[3,3])
  output[12,3] <- output[5,3]/output[7,3]
  output[13,3] <- output[6,3]/output[8,3]
  output[14,3] <- output[12,3]/output[13,3]
  
  return(output)
}



