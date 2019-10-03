setwd("/usr/share/sequencing/internships/JonathanBrenton/polyester_simulation/simulate")
library(doMC)
library(polyester)
library(Biostrings)
library(stringr)
library(dplyr)

### read the fasta file and imagine a coverage of 30 FPKM
### i.e. proportional to the transcripts length

writeLines("####### reading transcriptome FASTA")
gencode_fasta_file="/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/gencode.v32.transcripts.fa.gz"
gencode_fasta = readDNAStringSet(gencode_fasta_file)
GRCh38_fasta_file="/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.cdna.all.fa"
GRCh38_fasta = readDNAStringSet(GRCh38_fasta_file)

readspertx = round(30 * width(fasta) / 100)

###get the gene information separated

genes_gencode<- names(gencode_fasta)
genes_GRCh38<-names(GRCh38_fasta)

########Create a function to extract the transcript name, gene name and description from the gencode fasta information

parseGencode<-function(genestring){
  
identifier<-sub("^(ENST.+)\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+$", "\\1", genestring)
geneName<-sub("^ENST.+\\|.+\\|.+\\|.+\\|.+\\|(.+)\\|.+\\|.+$", "\\1", genestring)
description<-sub("^ENST.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|(.+)\\|$", "\\1", genestring)
return(c(identifier, geneName, description))
	  }
######apply the function to the gencode information to extract the terms above

extracted_gencode<-do.call("rbind", lapply(genes,parseGencode))

########Create a function to extract the transcript name, gene name and description from the GRCh36 cdna fasta information

parseGRCh38<- function(genestring){
  identifier<-sub("^(ENST.+)\\s.+\\s.+\\s.+\\s.+\\s.+\\s.+\\sdescription.+\\[.+\\]$", "\\1", genestring)
    geneName<-sub("^ENST.+\\s.+\\s.+\\s.+\\s.+\\s.+\\sgene_symbol:(.+)\\sdescription.+\\[.+\\]$", "\\1", genestring)
      description<-sub("^ENST.+\\s.+\\s.+\\s.+\\s.+\\s.+\\s.+\\sdescription:(.+)\\s\\[.+\\]$", "\\1", genestring)
        return(c(identifier, geneName, description))
	    }
######apply the function to the GRCh38 information to extract the terms above

extracted_GRCh38<-do.call("rbind", lapply(genes2,parseGRCh38))

colnames(extracted_GRCh38)<-c("Transcript_ID", "Gene_Name", "Description")
colnames(extracted_gencode)<-c("Transcript_ID", "Gene_Name", "Description")

extracted_gencode<-as.data.frame(extracted_gencode)

extracted_GRCh38<-as.data.frame(extracted_GRCh38)

merged_transcript_dataframe<-left_join(extracted_GRCh38, extracted_gencode, by=c("Transcript_ID", "Gene_Name"))


merged_transcript_dataframe<-na.omit(merged_transcript_dataframe)

immunoglobulin<-sample(which(grepl("immunoglobulin", merged_transcript_dataframe$Description.x)),200)
potassium<-sample(which(grepl("potassium", merged_transcript_dataframe$Description.x)),150)

#write.table(merged_transcript_dataframe[immunoglobulin,], "Immunoglobulin_genes_upregulated", sep=",",col.names = TRUE)
#write.table(merged_transcript_dataframe[potassium,], "Potassium_genes_upregulated", sep=",",col.names = TRUE)

##protein_coding_genes<-which(grepl("protein_coding", merged_transcript_dataframe$Description.y))
###protein_coding_genes<-merged_transcript_dataframe[protein_coding_genes,]

writeLines("##### preparing fold changes matrix")
## I then create a baseline fold matrix with 1
fold_changes <- matrix(
  rep(1,2*length(genes_gencode)),
    nrow = length(genes_gencode)
    )


## and substitute in the fold matrix the selected genes
## in cases and controls
## with random fold changes between 2 and 4
fold_changes[immunoglobulin,1] <- sample(2:4,length(immunoglobulin), replace = T)
fold_changes[potassium,2] <- sample(2:4,length(potassium), replace = T)

writeLines("#### simulating the reads")
## at this point I can simulate the experiment
simulate_experiment(fasta = gencode_fasta_file,outdir = "reads",num_reps = c(20,20),reads_per_transcript=readspertex, fold_changes=fold_changes,size = 0.2, error_model = "illumina5", gzip=TRUE,cores=16)
