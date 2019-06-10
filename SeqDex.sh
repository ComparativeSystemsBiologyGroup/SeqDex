#!/bin/bash
#$1 basename alignment file
#$2 basename contigs file
echo $1
echo $2

THREADS=10
#path to folder containing blast database files as downloaded from ftp://ftp.ncbi.nlm.nih.gov/blast/db/
#or any other custom database with sequence titles fulfilling NCBI sequence titles formatting rules
NT=~/database/nt 
#file name of the blast database in $NT
NTI=nt 
#path to folder containing nr database file in diamond format
#or any other custom database with sequence titles fulfilling NCBI sequence titles formatting rules
NR=~/database/
#file name of the nr diamond database in $NR
NRI=nrTaxonomy.dmnd 
#path to blast database built using $RDPF downloaded from https://rdp.cme.msu.edu/misc/resources.jsp
#or any other custom database with sequence titles fulfilling RDP sequence titles formatting rules
RDP=~/database/rdp16s
#fasta file used to build $RDPI
RDPF=current_Bacteria_unaligned.fa
#file name of the blast database in $RDP
RDPI=rdp16S
#home folder of SeqDex
SCRIPT=~/Downloads/
#taxonomy information to be used. Can be NT of NTNR
TAX=NT
#machine learning algorithm to be used. Can be RF, SVM, BOTH
MLALG=BOTH
#target taxonomy name to be used to identify the target cluster
TRG=Alphaproteobacteria
#if it is equal to 'yes', then SeqDex perform the final clustering step
CLUSTER=yes

##coverage calculation
echo "coverage calculation"

if [ -d Coverage ];then 
	echo "Using existing Coverage directory for the following steps. If this is not ok, please rename the directory.";

else 
	echo "The Coverage directory does not exist. Creating..."
	mkdir Coverage
fi


cd Coverage

if [ -f onlymapping_sorted_PAIRED.bed ];then 
	echo "Using existing coverage outputs for the following steps. If this is not ok, please remove the file.";

else 
	echo "The coverage outputs do not exist. Running...This will take a while"

	samtools view --threads "$THREADS" -b ../"${1}".sam -o "${1}".bam

	samtools view -F4 --threads "$THREADS" "${1}".bam -o "${1}"_onlymapping.bam

	samtools view -H ../"${1}".sam > header.sam

	samtools reheader header.sam "${1}"_onlymapping.bam > "${1}"_onlymapping_h.bam

	samtools sort --threads "$THREADS" "${1}"_onlymapping_h.bam > "${1}"_onlymapping_sorted.bam

	grep "@SQ" header.sam | cut -f2,3 -d ":" --output-delimiter "	" | sed -e 's/LN/1/g' > Scaffold_size.bed

	samtools view -F1 --threads "$THREADS" -O BAM "${1}"_onlymapping_sorted.bam | bedtools bamtobed  > "${1}"_SINGLE_end.bed

	samtools view -f1 --threads "$THREADS" -O BAM "${1}"_onlymapping_sorted.bam | bedtools bamtobed > "${1}"_PAIRED_end.bed

	coverageBed -a Scaffold_size.bed -bed -b "${1}"_SINGLE_end.bed > onlymapping_sorted_SINGLE.bed

	coverageBed -a Scaffold_size.bed -bed -b "${1}"_PAIRED_end.bed > onlymapping_sorted_PAIRED.bed

	samtools sort -n ../"${1}".sam -o ReadsRemap.bam -@ "$THREADS"

	samtools view ReadsRemap.bam > ReadsRemap.sam

	cut -f3,7 ReadsRemap.sam | grep -v "=" | grep -vi "*" | sort | uniq -c > contigNet.txt
fi
cd ..


##taxonomic affiliations & 16S
echo "taxonomy affiliation and 16S"

if [ -d Taxonomy ];then 
	echo "Using existing Taxonomy directory for the following steps. If this is not ok, please rename the directory.";

else 
	echo "The Taxonomy directory does not exist. Creating..."
	mkdir Taxonomy; 
fi

 
cd Taxonomy

if [ -f ContigsvsNt.txt ];then 
	echo "Using existing blastn output for the following steps. If this is not ok, please remove the file.";

else 
	echo "The blastn output does not exist. Running blastn...This will take a while"
	#blastn of the contigs file vs the database $NTI
	blastn -query ../"${2}".fasta -out ContigsvsNt.txt -db "$NT"/"${NTI}" -outfmt "6 std qlen slen qcovs gaps qcovhsp" -num_threads "$THREADS"; 
fi

if [ $TAX == "NTNR" ]; then
	if [ -f ContigsvsNr_mod.txt ];then 
		echo "Using existing diamond output for the following steps. If this is not ok, please remove the file.";

	else 
		echo "The diamond output does not exist. Running prodigal and diamond..."
		#blastn of the contigs file vs the database $NTI
		prodigal -i ../"${2}".fasta -a prodigalContigs.faa -f gff -p meta -q -o ProdContigs
diamond blastp --db "$NR"/"${NRI}" --query prodigalContigs.faa --out ContigsvsNr.txt --outfmt 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident gaps staxids stitle qcovhsp -p "$THREADS" --quiet
		cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,16 ContigsvsNr.txt > ContigsvsNr_mod.txt 
	fi
fi

echo "16S rRNA genes in contigs"
if [ -f RDP16s_taxa_mod.txt ];then 
	echo "Using existing 16S genes for the following steps. If this is not ok, please remove the file.";

else 
	echo "The 16S genese file does not exist. Running...This will take a while"

#find rRNAs in contigs
	barrnap -kingdom bac --threads "$THREADS" ../"${2}".fasta --quiet ON > barrnap16s_contigs.gff
#extract contigs with matches to 16S_rRNA
	grep -w 16S_rRNA barrnap16s_contigs.gff | cut -f1 > 16scontigsName.txt ; seqtk subseq ../"${2}".fasta 16scontigsName.txt > 16sContigs.fasta
#blastn contigs containing 16S rRNA against $RDPI
	blastn -query 16sContigs.fasta -out 16sContigvsRDP.txt -db "$RDP"/"${RDPI}" -outfmt "6 std qlen slen qcovs gaps qcovhsp" -num_threads "$THREADS"
#retrieving RDP taxonomy
	cut -f2 16sContigvsRDP.txt > RDP16s.txt; grep "$RDP"/"${RDPF}" -wf RDP16s.txt > RDP16s_taxa.txt
#adapting output
	sed "s/ /_/g" RDP16s_taxa.txt | sed "s/:/_/g"  | sed "s/#/_/g" | sed "s/'/_/g" | sed $"s/;/\t/g" | sed "s/>//g" > RDP16s_taxa_mod.txt
fi
cd ..
##Rscript

echo "Rscripts for kmers frequencies, SVM and final clustering step"
echo "Taxonomy"

cd Taxonomy
#create taxonomy tables
#build Taxonomy table with the following options (for full list run Rscript TaxonomyDef.R -h at the terminal)
#--taxaLevels 1 meaning contigs will be classified at superkingdom level
#--blast output of the blast of the contigs vs $NTI
#--aliLength minimum length of an HSP to be considered for taxonomy in blastn output
#--Xid minimum percentage of identity of an HSP to be considered for taxonomy in blastn output
#--TaxonDensity congruence score for taxonomy from blast
#--diamond output of the diamond of the prodigal predicted protein vs $NRI
#--aliLengthNR minimum length (in aa) of an HSP to be considered for taxonomy in Diamond output
#--XidNr minimum percentage of identity of an HSP to be considered for taxonomy in Diamond output
#--network network file 
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy
if [ $TAX == "NT" ]; then
	Rscript "$SCRIPT"/SeqDex/Taxonomy.R --blast ContigsvsNt.txt --taxaLevels 1 --aliLength 200 --Xid 70 --TaxonDensity 0.75 --network ../Coverage/contigNet.txt --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2
elif [ $TAX == "NTNR" ]; then
	Rscript "$SCRIPT"/SeqDex/Taxonomy.R --blast ContigsvsNt.txt --taxaLevels 1 --aliLength 200 --Xid 70 --TaxonDensity 0.75 --diamond ContigsvsNr_mod.txt --aliLengthNr 70 --XidNr 80 --network ../Coverage/contigNet.txt --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2
fi
 
#create tables with taxonomy assignments for 16S rRNA contigs (for full list run Rscript rRNA16S.R -h at the terminal)
#--blastRDP output of the blast of the 16S contigs vs RDP
#--taxaRDP file with taxonomy from $RDPF
#--minContigLeng minimum length of contigs considered by the model (should be the same of SVM/RF)
Rscript "$SCRIPT"/SeqDex/rRNA16S.R --blastRDP 16sContigvsRDP.txt --taxaRDP RDP16s_taxa_mod.txt --minContigLen 1000
cd ..
echo "gck"
cd Coverage
if [ -f gckCovTable.txt ];then 
	echo "Using existing K-mers frequency table. If this is not ok, please remove the file.";

else 
	echo "The k-mers frequency table does not exists, running Rscript"
#create tables with kmers frequencies, GC content and contigs coverage 
#--contigs contigs fasta file
#--kmers length of the kmer to use for machine learning
#--covSingle --covPaired files created during coverage calculation (see above)
#--threads number of threads to be used 
	Rscript "$SCRIPT"/SeqDex/GCKmersCov.R --contigs ../"${2}".fasta --Kmers 3 --covSingle onlymapping_sorted_SINGLE.bed --covPaired onlymapping_sorted_PAIRED.bed --threads "$THREADS"
fi
cd ..

if [ $MLALG == "SVM" ]; then
	echo "SVM"
	mkdir SVMoutput
	cd SVMoutput
#predicting taxonomic affiliation using SVM based upon output of TaxonomyDef.R (for full list run Rscript SVMDef.R -h at the terminal)
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--cross --cost --gamma --scale parameters of SVM algorithm
#--minContigLeng minimum length of contigs to be considerend by SVM
#--targetName name of the target taxonomy of interest
#--TaxaName taxonomy level of interest
#--threads number of threads to be used
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy
#--verbose set to true to print on the terminal performance statistics and percentage of error committed by SVM model 
	Rscript "$SCRIPT"/SeqDex/SVM.R --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt --network ../Coverage/contigNet.txt --cross 5 --cost 1e-1,1e3 --gamma 1e-5,1e-1 --scale F --minContigLen 1000 --targetName Bacteria --TaxaName superkingdom --threads "$THREADS" --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2 --verbose T
	cd ..

	if [ $CLUSTER == "yes" ]; then

		echo "umap e dbscan"

		mkdir ClusteringOutput
		cd ClusteringOutput
#clustering contigs of the target taxonomy found with SVMDef.R. The  cluster of interest is identified by using 16S rRNA of target taxonomy with higher coverage 
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--input variables to be used in clustering algorithm
#--minContigLeng minimum length of contigs to be considerend by SVM
#--targetName name of the target taxonomy of interest
#--threads number of threads to be used
#--ncomponents number of components to be extract by umap algorithm
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy

		Rscript "$SCRIPT"/SeqDex/Clustering --modelOutput ../SVMoutput/1outputSVM.txt --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt  --network ../Coverage/contigNet.txt --input Kmers --minContigLen 1000 --targetName "$TRG" --threads "$THREADS" --ncomponents 2 --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2 

		seqtk subseq ../"${2}".fasta OutputClustering.txt > OutputClustering.fasta
		cd ..
	else
		cd SVMoutput
		cut -f1 1outputSVM.txt > 1outputSVM_names.txt
		seqtk subseq ../"${2}".fasta 1outputSVM_names.txt > 1outputSVM_contigs.fasta
		cd ..
	fi

elif [ $MLALG == "RF" ]; then

	echo "RF"
	mkdir RFoutput
	cd RFoutput
#predicting taxonomic affiliation using SVM based upon output of TaxonomyDef.R (for full list run Rscript RFDef.R -h at the terminal)
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--minContigLeng minimum length of contigs to be considerend by RF
#--targetName name of the target taxonomy of interest
#--TaxaName taxonomy level of interest
#--replace if sampling cases of RandomForest should be done with or without replacement
#--ntree number of tree to be grown in te forest
#--threads number of threads to be used
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy
#--verbose set to true to print on the terminal performance statistics and percentage of error committed by RF model 
	Rscript "$SCRIPT"/SeqDex/RF.R --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt --network ../Coverage/contigNet.txt --minContigLen 1000 --targetName Bacteria --TaxaName superkingdom --replace T --ntree 500 --threads "$THREADS"  --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2 --verbose T
	cd ..
	if [ $CLUSTER == "yes" ]; then
		echo "umap e dbscan"

		mkdir ClusteringOutput
		cd ClusteringOutput
#clustering contigs of the target taxonomy found with RFDef.R. The  cluster of interest is identified by using 16S rRNA of target taxonomy with higher coverage 
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--input variables to be used in clustering algorithm
#--minContigLeng minimum length of contigs to be considerend by RF
#--targetName name of the target taxonomy of interest
#--threads number of threads to be used
#--ncomponents number of components to be extract by umap algorithm
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy

		Rscript "$SCRIPT"/SeqDex/Clustering.R --modelOutput ../RFoutput/1outputRF.txt --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt  --network ../Coverage/contigNet.txt --input Kmers --minContigLen 1000 --targetName "$TRG" --threads "$THREADS" --ncomponents 2 --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2
		seqtk subseq ../"${2}".fasta OutputClustering.txt > OutputClustering.fasta
		cd ..
	else
		cd RFoutput
		cut -f1 1outputRF.txt > 1outputRF_names.txt
		seqtk subseq ../"${2}".fasta 1outputRF_names.txt > 1outputRF_contigs.fasta
		cd ..
	fi
elif [ $MLALG == "BOTH" ]; then
	mkdir SVMoutput
	cd SVMoutput
	echo "SVM"
#predicting taxonomic affiliation using SVM based upon output of TaxonomyDef.R (for full list run Rscript SVMDef.R -h at the terminal)
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--cross --cost --gamma --scale parameters of SVM algorithm
#--minContigLeng minimum length of contigs to be considerend by SVM
#--targetName name of the target taxonomy of interest
#--TaxaName taxonomy level of interest
#--threads number of threads to be used
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy
#--verbose set to true to print on the terminal performance statistics and percentage of error committed by SVM model 
	Rscript "$SCRIPT"/SeqDex/SVM.R --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt --network ../Coverage/contigNet.txt --cross 5 --cost 1e-1,1e3 --gamma 1e-5,1e-1 --scale F --minContigLen 1000 --targetName Bacteria --TaxaName superkingdom --threads "$THREADS" --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2 --verbose T
	cd ..
	if [ $CLUSTER == "yes" ]; then
		echo "umap e dbscan on SVM"

		mkdir ClusteringOutputSVM
		cd ClusteringOutputSVM
#clustering contigs of the target taxonomy found with SVMDef.R. The  cluster of interest is identified by using 16S rRNA of target taxonomy with higher coverage 
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--input variables to be used in clustering algorithm
#--minContigLeng minimum length of contigs to be considerend by SVM
#--targetName name of the target taxonomy of interest
#--threads number of threads to be used
#--ncomponents number of components to be extract by umap algorithm
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy

		Rscript "$SCRIPT"/SeqDex/Clustering.R --modelOutput ../SVMoutput/1outputSVM.txt --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt  --network ../Coverage/contigNet.txt --input Kmers --minContigLen 1000 --targetName "$TRG" --threads "$THREADS" --ncomponents 2 --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2
		seqtk subseq ../"${2}".fasta OutputClustering.txt > OutputClusteringSVM.fasta
		cd ..
	else
		cd SVMoutput
		cut -f1 1outputSVM.txt > 1outputSVM_names.txt
		seqtk subseq ../"${2}".fasta 1outputSVM_names.txt > 1outputSVM_contigs.fasta
		cd ..
	fi

	echo "RF"
	mkdir RFoutput
	cd RFoutput
#predicting taxonomic affiliation using SVM based upon output of TaxonomyDef.R (for full list run Rscript RFDef.R -h at the terminal)
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--minContigLeng minimum length of contigs to be considerend by RF
#--targetName name of the target taxonomy of interest
#--TaxaName taxonomy level of interest
#--ntree number of tree to be grown in te forest
#--replace if sampling cases of RandomForest should be done with or without replacement
#--threads number of threads to be used
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy
#--verbose set to true to print on the terminal performance statistics and percentage of error committed by RF model 
	Rscript "$SCRIPT"/SeqDex/RF.R --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt --network ../Coverage/contigNet.txt --minContigLen 1000 --targetName Bacteria --TaxaName superkingdom --replace T --ntree 500 --threads "$THREADS" --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2 --verbose T
	cd ..
	if [ $CLUSTER == "yes" ]; then

		echo "umap e dbscan on RF"

		mkdir ClusteringOutputRF
		cd ClusteringOutputRF
#clustering contigs of the target taxonomy found with RFDef.R. The  cluster of interest is identified by using 16S rRNA of target taxonomy with higher coverage 
#--gckCovTable path and filename of the output of GCKmersCovDef.R files
#--rRNA16S path and filename to output of rRNA16S.R
#--taxonomy path and filenameto output of Taxonomy.R
#--network network file 
#--input variables to be used in clustering algorithm
#--minContigLeng minimum length of contigs to be considerend by RF
#--targetName name of the target taxonomy of interest
#--threads number of threads to be used
#--ncomponents number of components to be extract by umap algorithm
#--Edges minimum edge value to be considered 
#--VerticesDegree maximum vertices degree to be considered 
#--componentSize minimum component size 
#--mixedComponents maximum proportion of alternative taxonomy over the total to consider the alternative as erroneous and correct/uniform the CC taxonomy

		Rscript "$SCRIPT"/SeqDex/Clustering.R --modelOutput ../RFoutput/1outputRF.txt --gcCovKmersTable ../Coverage/gckCovTable.txt --rRNA16S ../Taxonomy/rRNA16sTaxonomy2.txt --taxonomy ../Taxonomy/1taxonomyIteration.txt  --network ../Coverage/contigNet.txt --input Kmers --minContigLen 1000 --targetName "$TRG" --threads "$THREADS" --ncomponents 2 --Edges 10 --VerticesDegree 5 --componentSize 2 --mixedComponents 0.2
 
		seqtk subseq ../"${2}".fasta OutputClustering.txt > OutputClusteringRF.fasta
		cd ..
	else
		cd RFoutput
		cut -f1 1outputRF.txt > 1outputRF_names.txt
		seqtk subseq ../"${2}".fasta 1outputRF_names.txt > 1outputRF_contigs.fasta
		cd ..
	fi
fi





