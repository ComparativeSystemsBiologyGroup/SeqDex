#!/Desktop/
path <- list.files(path= "~", full.names=TRUE,
                   recursive=TRUE,pattern="(Func.R)")
pathTrue <- path[grep(pattern = "Function/Func.R", path, fixed = T)]
source(pathTrue) #path to Func.R file


suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(igraph, quietly = T))
#Option------------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-b", "--blast"), action = "store", type = "character", default = "ContigsvsNt.txt",
              help = "contigs vs balst nt path/filename.ext - [default %default]"),
  make_option(c("-d", "--diamond"), action = "store", type = "character", default = NA,
              help = "contigs vs diamond nr path/filename.ext"),
  make_option(c("-R", "--network"), action = "store", type = "character", default = "../Coverage/contigNet.txt",
              help = "reads based graph network path/filename.ext - [default %default]"),
  #make_option(c("-t", "--taxa"), action = "store", type = "character", default = "1,2,3,4,5,6,7",
             # help = "taxa levels to be considered when parsing blast/diamond output:\n\t\t\t 1-superkingdom;\n\t\t\t 2-phylum;\n\t\t\t 3-class;\n\t\t\t 4-order;\n\t\t\t 5-family;\n\t\t\t 6-genus;\n\t\t\t 7-species\n\t\t [default %default]"),
  make_option(c("-L", "--taxaLevels"), action = "store", type = "character", default = "Superkingdom",
              help = "taxonomy levels to be considered to perform SeqDex; a comma separated list of values, one for each iteration - [default %default]"),
  make_option(c("-a", "--aliLength"), action = "store", type = "integer", default = 200,
              help = "minimum alignment length (nts) to be considered in blast output parsing - [default %default]"),
  make_option(c("-i", "--Xid"), action = "store", type = "integer", default = 70, 
              help = "minimum percendtage of identity to be considered when parsing blast output - [default %default]"),
  make_option(c("-A", "--aliLengthNr"), action = "store", type = "integer", default = 70,
              help = "minimum alignment length (aa) to be considered in diamond output parsing - [default %default]"),
  make_option(c("-I", "--XidNr"), action = "store", type = "integer", default = 90, 
              help = "minimum percendtage of identity to be considered when parsing diamond output - [default %default]"),
  make_option(c("-D", "--TaxonDensity"), action = "store", type = "double", default = 0.75,
              help = "minimum TaxonDensity level to be considered in the output of Taxonomy.R - [default %default]"),
  make_option(c("-E", "--Edges"), action = "store", type = "integer", default = 10,
              help = "minimum edges value to be considered to construct network to extend the taxonomy - [default %default]"),
  make_option(c("-V", "--VerticesDegree"), action = "store", type = "integer", default = 5,
              help = "maximum vertex degree of the network to be considered to extract component and extend taxonomy - [default %default]"),
  make_option(c("-S", "--componentSize"), action = "store", type = "integer", default = 2,
              help = "minimum size of the components to be considered to extend taxonomy - [default %default]"),
  make_option(c("-M", "--mixedComponents"), action = "store", type = "double", default = 0.2,
              help = "threshold of frequency of alternative taxonomy over the total of the component. \n\t\t\t Under this frequency value the alternative will be corrected and substituded by majority - [default %default]"),
  make_option(c("-x", "--taxTransfer"), action = "store", type = "character", default = "all",
              help = "if All means that taxonomy will be transferred to all contigs in the CC. To limit the transfer write the maximum distance desired - [default %default]")
  
)
opt <- parse_args(OptionParser(option_list = option_list))


#USER difined parameters-------------------------------------------------------------------------

if (grepl(opt$taxaLevels, pattern = ",")) {
  TAXAlevel <- tolower(as.character(unlist(strsplit(opt$taxaLevels, split = ",", fixed = T))))
} else {
  TAXAlevel <- tolower(as.character(opt$taxaLevels))
}

if (opt$TaxonDensity<0.51) {
  opt$TaxonDensity <- 0.51
  print("minimum possible value of --TaxonDensity is 0.51. User defined value changed to this value.")
}

#Blast vs Nt-------------------------------------------------------------------------------------
Blast <-read.table(opt$blast,#"../sacffolds91vsntAll.txt", # , #blast vs Nt output file name (in working dir)
                    header= FALSE,
                    col.names = c("query",	"subject",	"X.id",	"ali_Length",	"X1",	"X2",	"STARTQ",	
                                  "ENDQ",	"STARTS", "ENDS",	"EVALUE",	"BITSCORE",	"QLEN",	"SLEN",
                                  "QCOVS",	"GAPS",	"QCOVHSP"),
                    stringsAsFactors = FALSE)

BlastFilt <- Blast[which(Blast$ali_Length>opt$aliLength),]

suppressWarnings(BlastAll <-BlastTaxaAll(table=BlastFilt,  taxa=TAXAlevel, 
                                         blast_input = "chromosome"
                                         ))

colnames(BlastAll) <-c("Contig", colnames(BlastAll)[2:ncol(BlastAll)]) #colnames
#write.table(BlastAll, "Nt.txt", sep = "\t", quote = F, row.names = F)
taxonomy.df <- data.frame(taxa=c("superkingdom", " phylum", "class", "order", "family", "genus", "species"),
                          ntaxa = c(1:7), stringsAsFactors = F)

taxonDensityLevel = opt$TaxonDensity
#BlastAll <- read.table("Nt.txt", sep = "\t", stringsAsFactors = FALSE) #uncomment to restore work
BlastAllNoNa <- BlastAll[complete.cases(BlastAll),]
BlastNoNaFilterLength <- BlastAllNoNa[which(BlastAllNoNa$ali_Length>opt$aliLength & 
                                              BlastAllNoNa$X.id>=opt$Xid),]

#Diamond vs Nr------------------------------------------------------------------
if (!is.na(opt$diamond)) {
  DiamondTaxa <- read.table(opt$diamond,#"../ProdigalScaffolds91vsnrTaxDiamond.txt", #, 
                            sep = "\t",
                            header = FALSE,
                            stringsAsFactors = FALSE,
                            quote="'\"",
                            col.names = c("query", "QLEN", "subject", "SLEN", "STARTQ", "ENDQ",
                                        "STARTS", "ENDS", "EVALUE", "BITSCORE", "ali_Length",
                                        "X.id", "GAPS", "X", #"STITLE", 
                                        "QCOVHSP")
                            )
  suppressWarnings(DiamondBlastp <-BlastTaxaAll(table=DiamondTaxa, taxa=TAXAlevel, 
                                                blast_input = "chromosome" 
                                                ))
  colnames(DiamondBlastp) <-c("Contig", colnames(DiamondBlastp)[2:ncol(DiamondBlastp)]) #colnames
  name <- strsplit(DiamondBlastp$Contig, "_", fixed = TRUE)
  name2 <- paste(unlist(lapply(name, "[", 1)),
                unlist(lapply(name, "[", 2)), 
                unlist(lapply(name, "[", 3)),
                unlist(lapply(name, "[", 4)),
                unlist(lapply(name, "[", 5)),
                unlist(lapply(name, "[", 6)),
                sep = "_"
                )
  DiamondBlastp$Contig <- name2
  #write.table(DiamondBlastp, "Taxonomy/Nr.txt", sep = "\t", quote = F, row.names = F)
  DiamondBlastNoNa <- DiamondBlastp[complete.cases(DiamondBlastp),]
  DiamondNoNaFilterLength <- DiamondBlastNoNa[which(DiamondBlastNoNa$ali_Length>opt$aliLengthNr & DiamondBlastNoNa$X.id>opt$XidNr),]
}

#taxonomy graph---------------------------------------------------------------------------------------
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



#taxonomy Table ---------------------------------------------------------------------------------------


for (i in 1:length(TAXAlevel)) {
  if (is.na(opt$diamond)) {
    taxa.level= which(colnames(BlastAll)==taxonomy.df[which(taxonomy.df$taxa == TAXAlevel[i]),1])
    BlastTD <- TaxonomyDensity(BlastNoNaFilterLength, taxa.level)
    BTDallFilt75 <- BlastTD[which(BlastTD$TaxonDensity>=opt$TaxonDensity),]
    compTax <- as.data.frame(matrix(data = NA, nrow = length(lcomp), 
                                    ncol = (2 + length(unique(BTDallFilt75[,3])))))
    colnames(compTax) <- c("comp", "n", unique(BTDallFilt75[, 3]))
    
    
    for (j in 1:length(lcomp)) {
      compTax[j,1] <- lcomp[j]
      compTax[j,2] <- length(which(u$membership==j))
      tmp <-which(u$membership==j)
      for (m in 3:ncol(compTax)) {
        compTax[j,m] <- length(which(BTDallFilt75$Contig%in%rownames(u)[tmp] & 
                                       BTDallFilt75[,3]==colnames(compTax)[m]))
      }
      
    }
    
    for (n in 1:nrow(compTax)) {
      if (sum(compTax[n,3:ncol(compTax)])==0) {
        compTax[n,] <- NA
      } else if (sum(compTax[n,3:ncol(compTax)])==compTax[n,2])
        compTax[n,] <- NA
    }
    compTax<- compTax[complete.cases(compTax),]
    compTaxProp <- cbind(compTax[,1:2], compTax[,3:ncol(compTax)]/compTax[,2])
    
    #taxCong<-matrix(0,nrow=nrow(compTax),ncol=1)
    
    newtax<-NULL
    
    for(k in 1:nrow(compTax)){
      
      contigs_in_comp<-unames[which(u$membership==compTax[k,1])]
      tax<-BTDallFilt75[which(BTDallFilt75$Contig %in% contigs_in_comp),]
      l<-length(unique(tax[,3]))
      if (tolower(opt$taxTransfer)=="all"){
        s<-setdiff(contigs_in_comp,tax[,1])
      } else {
        s<-setdiff(contigs_in_comp,tax[,1])  
        ego(G.copy, tax$Contig, mindist = 0, order = opt$taxTransfer)
        egoTOT <- NULL
        for (x in 1:length(egoTax)) {
          egoTaxDF <- as.matrix(unlist(egoTax[[x]]))
          egoTOT <- rbind(egoTOT, egoTaxDF)
        }
        egoTOT <- unique(egoTOT)
        
        s <- s[which(s %in% rownames(egoTOT))]
        }
      
      
      
      if (compTax$n[k]==2 && l==1){
        newtax<-rbind(newtax,c(s,"-1",tax[1,3]))
      } else if (compTax$n[k]>2 && l==1) {
        tmp <- cbind(s, rep("-1", length(s)), rep(tax[1,3], length(s)))
        newtax<-rbind(newtax,tmp)
        } else if (compTax$n[k]>2 && l>1) {
          if(min(compTax[k,3:ncol(compTax)])==1 && compTaxProp[k,which(colnames(compTaxProp)%in% unique(tax[,3]))]>opt$mixedComponents){
          tmp <- cbind(s, rep("-1", length(s)), rep(tax[1,3], length(s)))
          newtax<-rbind(newtax,tmp)
        }
      }
      
    }
    
    colnames(newtax) <- colnames(BTDallFilt75)
    mergedNewTax <- rbind(BTDallFilt75, newtax)
    
    
    
    write.table(mergedNewTax, paste(TAXAlevel[i],"TaxonomyIteration.txt", sep = ""), sep = "\t", quote = F, row.names = F)
  } else {
    taxa.level= which(colnames(BlastAll)==taxonomy.df[which(taxonomy.df$taxa == TAXAlevel[i]),1])
    BlastTD <- TaxonomyDensity(BlastNoNaFilterLength, taxa.level)
    BTDallFilt75 <- BlastTD[which(BlastTD$TaxonDensity>=opt$TaxonDensity),]
    taxa.levelD= which(colnames(DiamondNoNaFilterLength)==taxonomy.df[which(taxonomy.df$taxa == TAXAlevel[i]),1])
    DiamondFiltTD <- TaxonomyDensity (DiamondNoNaFilterLength, taxa.levelD)
    DiamondTDFilt75 <- DiamondFiltTD[which(DiamondFiltTD$TaxonDensity>=opt$TaxonDensity),]
    TaxonomyAll <- merge(BTDallFilt75, DiamondTDFilt75, by = "Contig", all = FALSE, suffixes = c("B", "D"))
    
    BlastTaxaNum = which(colnames(TaxonomyAll)==paste(taxonomy.df[TAXAlevel[i],1], "B", sep = ""))
    DiamondTaxaNum = which(colnames(TaxonomyAll)==paste(taxonomy.df[TAXAlevel[i],1], "D", sep = ""))
    
    Shared <- TaxonomyAll[which(TaxonomyAll[,BlastTaxaNum]==TaxonomyAll[,DiamondTaxaNum]),]
    TaxonMax <- c(rep(NA, nrow(Shared)))
    for (k in 1:nrow(Shared)) {
      TaxonMax[k] <- max(Shared$TaxonDensityB[k], Shared$TaxonDensityR[k], Shared$TaxonDensity[k])
    }
    Shared<- data.frame(Contig = Shared$Contig,
                        TaxonDensity = TaxonMax,
                        taxa.name = Shared[,BlastTaxaNum], stringsAsFactors = FALSE)
    colnames(Shared)[3]<-taxonomy.df[TAXAlevel[i],1] 
    
    
    NotSharedB <- BTDallFilt75[-which(BTDallFilt75$Contig %in% Shared$Contig | 
                                           BTDallFilt75$Contig %in% DiamondTDFilt75$Contig),]
    NotSharedD <- DiamondTDFilt75[-which(DiamondTDFilt75$Contig %in% Shared$Contig |
                                           DiamondTDFilt75$Contig %in% BTDallFilt75$Contig),]
    
    mergedTaxonomy <- rbind.data.frame(Shared, NotSharedB, NotSharedD, 
                                          stringsAsFactors = FALSE)
    
    compTax <- as.data.frame(matrix(data = NA, nrow = length(lcomp), 
                                    ncol = (2 + length(unique(mergedTaxonomy[,3])))))
    colnames(compTax) <- c("comp", "n", unique(mergedTaxonomy[,3]))
    
    
    for (j in 1:length(lcomp)) {
      compTax[j,1] <- lcomp[j]
      compTax[j,2] <- length(which(u$membership==j))
      tmp <-which(u$membership==j)
      for (m in 3:ncol(compTax)) {
        compTax[j,m] <- length(which(mergedTaxonomy$Contig%in%rownames(u)[tmp] & 
                                       mergedTaxonomy[,3]==colnames(compTax)[m]))
      }
      
    }
    
    for (n in 1:nrow(compTax)) {
      if ( sum(compTax[n,3:ncol(compTax)])==0) {
        compTax[n,] <- NA
      } else if (sum(compTax[n,3:ncol(compTax)])==compTax[n,2])
        compTax[n,] <- NA
    }
    compTax<- compTax[complete.cases(compTax),]
    compTaxProp <- cbind(compTax[,1:2], compTax[,3:ncol(compTax)]/compTax[,2])
    
    #taxCong<-matrix(0,nrow=nrow(compTax),ncol=1)
    
    newtax<-NULL
    
    for(k in 1:nrow(compTax)){
      
      contigs_in_comp<-unames[which(u$membership==compTax[k,1])]
      tax<-mergedTaxonomy[which(mergedTaxonomy$Contig %in% contigs_in_comp),]
      l<-length(unique(tax[,3]))
      if (tolower(opt$taxTransfer)=="all"){
        s<-setdiff(contigs_in_comp,tax[,1])
      } else {
        s<-setdiff(contigs_in_comp,tax[,1])  
        ego(G.copy, tax$Contig, mindist = 0, order = opt$taxTransfer)
        egoTOT <- NULL
        for (x in 1:length(egoTax)) {
          egoTaxDF <- as.matrix(unlist(egoTax[[x]]))
          egoTOT <- rbind(egoTOT, egoTaxDF)
        }
        egoTOT <- unique(egoTOT)
        
        s <- s[which(s %in% rownames(egoTOT))]
      }
      
      
      if (compTax$n[k]==2 && l==1){
        newtax<-rbind(newtax,c(s,"-1",tax[1,3]))
      } else if (compTax$n[k]>2 && l==1) {
        tmp <- cbind(s, rep("-1", length(s)), rep(tax[1,3], length(s)))
        newtax<-rbind(newtax,tmp)
        
      } else if (compTax$n[k]>2 && l>1) {
        if(min(compTax[k,3:ncol(compTax)])==1 && compTaxProp[k,which(colnames(compTaxProp) %in% unique(tax[,3]))]>opt$mixedComponents){
          tmp <- cbind(s, rep("-1", length(s)), rep(tax[1,3], length(s)))
          newtax<-rbind(newtax,tmp)
        }
      }
      
    }
    
    colnames(newtax) <- colnames(mergedTaxonomy)
    mergedNewTax <- rbind(mergedTaxonomy, newtax)
    
    
    
    write.table(mergedNewTax, paste(TAXAlevel[i], "TaxonomyIteration.txt", sep = ""), sep = "\t", quote = F, row.names = F)
    
  }
  #print(i)
}


