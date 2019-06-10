# SeqDex

A model to deconvolve endosymbionts genomic sequences insides whole genome sequences of both host and symbionts. It is written in bash and is composed by bash package and R scripts. It can be executed by running both the .sh script or script independently in bash. It has been developed in Linux.

## Download
SeqDex does not need to be installed, just download and decompress the directory. To run SeqDex.sh, the first time only run

chmod 755 SeqDex.sh

Inside the SeqDex folder

## Dependencies

### Mandatory:

    • Bash:
        ◦ Samtools
        ◦ Bedtools
        ◦ NCBI-BLAST+
        ◦ Barrnap
        ◦ Seqtk
    • R package:
        ◦ Taxonomizr
            ▪ See vignette at https://cran.r-project.org/web/packages/taxonomizr/vignettes/usage.html . 
            SeqDex needs the accessionTaxa.sql file, so follow the instruction to download it. 
            Note that the command prepareDatabase(‘accessionTaxa.sql’) only downloads 
            nucleotide accession information of NCBI. If you wish to use SeqDex using also 
            proteic taxonomic affiliations, then follow the Manual preparation of database 
            instructions to download a complete accessionTaxa.sql file.
        ◦ Seqinr
        ◦ randomForest
        ◦ e1071
        ◦ Uwot
        ◦ Dbscan
        ◦ Parallel
        ◦ doParallel
        ◦ Foreach
        ◦ Optparse
        ◦ Ggplot2
        ◦ igraph

### Optional:

These dependences need to be installed if you wish to run SeqDex using also using taxonomic affiliation obtained by align 
sequences to protein database
    • Bash
        ◦ Prodigal
        ◦ Diamond

###Databases

SeqDex needs databases to obtain taxonomic affiliations. In detail, two are mandatory: a nucleotidic database, 
either NCBI nt or a custom made database in blast format with NCBI titles; RDP bacterial 16S (unaligned), 
both the fasta file and the blast format database. If you wish to use protein, also a protein database, either nr NCBI or a 
custom database with NCBI titles, in diamond format. Nr NCBI database is available only in blast format, but to improve 
computational time, we choose not to use blastp. 

## Running

SeqDex needs two input files: contigs file resulting from assembly (fasta format)
and the alignment file obtained through alignment of the assembly reads to the contigs (sam format). 
You can choose which program use: we tested SeqDex by assembling with SPADes and mapping using Bowtie, but it needs just an 
alignment file in sam format and assembly in fasta format.

### Quick Run:

To run SeqDex first open the SeqDex.sh file with you text editor to insert needed path to nucleotide, protein, 16S databases and 
homefolder of SeqDex. If you wish, you can save the sh file with another name, just do not forget to rerun chmod command first.
Then copy the sh in the folder where you want to run SeqDex and just run the command below 

    ./SeqDex.sh basename_alignment basename_contig

This will run SeqDex by using:
    • nucleotidic taxonomic affiliations
    • 3-mers frequencies
    • SVM and RF machine learning algorithm only at superkingdom taxonomic level
    • Clustering by searching for Alphaproteobacteria 16S gene with higher coverage

### Custom Run:

You can run SeqDex by using also proteic taxonomic affiliation and changing machine learning algorithm by modifying TAX and MLALG
variables in SeqDex.sh. Also, the R scripts are highly flexible and allow to run iteratively the machine learning predictive step 
first to perform the final clustering. To do so, please first read the manual and/or run in bash terminal 
    
    Rscript name_script.R -h 
    
to see the list of all available option for each R script of SeqDex. 
Then, you can choose to modify the SeqDex.sh file or to run the R scripts independently.

## Output

SeqDex will produce various files. Most of them are used to deconvolving sequences. 
Sequences of putative target endosymbiont are listed in fasta file ClusteringOutput/OutputClustering.fasta

It contains only the sequences of the cluster that contains the 16S rRNA of interest.
