setwd("/usr/share/sequencing/internships/JonathanBrenton/polyester_simulation/simulate/polyester_tests_13.10.19_excit+inhib_genes")
library(doMC)
library(polyester)
library(Biostrings)
library(stringr)
library(dplyr)

### read the fasta file and imagine a coverage of 30 FPKM
### i.e. proportional to the transcripts length

writeLines("####### reading transcriptome FASTA")

gencode_fasta_file="/usr/share/sequencing/references/to_be_deleted/rnaseq/Homo_sapiens/ENSEMBL/gencode.v32.transcripts.fa.gz"
gencode_fasta = readDNAStringSet(gencode_fasta_file)

readspertx = round(30 * width(gencode_fasta) / 100)

###get the gene information separated

genes_gencode<- names(gencode_fasta)

parseGencode<-function(genestring){

  identifier<-sub("^(ENST.+)\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+$", "\\1", genestring)
  geneName<-sub("^ENST.+\\|.+\\|.+\\|.+\\|.+\\|(.+)\\|.+\\|.+$", "\\1", genestring)
  description<-sub("^ENST.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|.+\\|(.+)\\|$", "\\1", genestring)
  return(c(identifier, geneName, description))
}

######apply the function to the gencode information to extract the terms above

extracted_gencode<-do.call("rbind", lapply(genes_gencode,parseGencode))

colnames(extracted_gencode)<-c("Transcript_ID", "Gene_Name", "Description")

extracted_gencode<-as.data.frame(extracted_gencode)

#gencode_txid_to_geneid<-read.delim(file ="genecode_txid_to_geneid.txt", sep=" ", header=FALSE)
#colnames(gencode_txid_to_geneid)<-c("Transcript_ID", "Gene_ID","Gene_Name")


read.csv(file="inhibitory_updated_list.csv", header=TRUE, skip=1)-> inhibitory_list

inhibitory_list<-inhibitory_list %>% filter(!str_detect(Match.type, 'Unmatched'))

read.csv(file="excitatory_updated_list.csv", header=TRUE, skip=1)-> excitatory_list

excitatory_list<-excitatory_list %>% filter(!str_detect(Match.type,"Unmatched"))

### C6ORF165 was duplicated so:
excitatory_list<-excitatory_list[!duplicated(excitatory_list$Approved.symbol),]

####mismatch between ensembl and HUGO-KIF1BP vs KIFBP, because using gencode had to change

levels(inhibitory_list$Approved.symbol) <- c(levels(inhibitory_list$Approved.symbol), 'KIF1BP')
inhibitory_list$Approved.symbol[inhibitory_list$Approved.symbol == 'KIFBP'] <- 'KIF1BP'

## intersect(excitatory_list$Approved.symbol, inhibitory_list$Approved.symbol)
##[1] "GALNT18" "MFSD4B"  "MLIP"    "NECTIN3" "PYHIN1"  "ZBTB18"

##so to remove overlapping :
overlap<-intersect(excitatory_list$Approved.symbol, inhibitory_list$Approved.symbol)
excitatory_list<-excitatory_list[excitatory_list$Approved.symbol %in% setdiff(excitatory_list$Approved.symbol,overlap),]
inhibitory_list<-inhibitory_list[inhibitory_list$Approved.symbol %in% setdiff(inhibitory_list$Approved.symbol,overlap),]
## intersect(inhibitory_list$Approved.symbol,excitatory_list$Approved.symbol)
###character(0)



#encode_txid_to_geneid$Transcript_ID<-sub("\\..+", "", gencode_txid_to_geneid$Transcript_ID)
###to avoid version discrepancies
extracted_gencode$Transcript_ID<-sub("\\..+", "", extracted_gencode$Transcript_ID)

excitatory_names_and_gene_matches<-extracted_gencode %>% filter(Gene_Name %in% excitatory_list$Approved.symbol)

inhibitory_names_and_gene_matches<-extracted_gencode %>% filter(Gene_Name %in% inhibitory_list$Approved.symbol)


write.table(excitatory_names_and_gene_matches, "Excitatory_genes_transcripts_upregulated.csv", col.names =TRUE, sep=" ", row.names=FALSE)
write.table(inhibitory_names_and_gene_matches, "Inhibitory_genes_transcripts_upregulated.csv", col.names =TRUE, sep=" ", row.names=FALSE)

excit_positions<-which(extracted_gencode$Transcript_ID %in% excitatory_names_and_gene_matches$Transcript_ID)
inhib_positions<-which(extracted_gencode$Transcript_ID %in% inhibitory_names_and_gene_matches$Transcript_ID)

###for reference so to show works
#> unique(extracted_gencode[inhib_positions,2])->aa
#> unique(inhibitory_names_and_gene_matches$Gene_Name)->bb
#> identical(aa,bb)
#[1] TRUE
##same with excitatory


writeLines("##### preparing fold changes matrix")
## I then create a baseline fold matrix with 1
fold_changes <- matrix(
  rep(1,2*length(genes_gencode)),
  nrow = length(genes_gencode)
)


fold_changes[excit_positions,1] <- sample(2:4,length(excit_positions), replace = T)
fold_changes[inhib_positions,2]  <- sample(2:4,length(inhib_positions), replace = T)
## identical(which(fold_changes[,2]>1),inhib_positions)
##[1] TRUE


simulate_experiment(fasta = gencode_fasta_file, 
                    outdir = "excit.inhib.polyester.test.size-0.3",
                    num_reps = c(20,20),
                    reads_per_transcript=readspertx,
                    fold_changes=fold_changes,
                    size = 0.3,
                    error_model = "illumina5",
                    gzip=TRUE,
                    ### if doMC is installed here you can add
                    cores=32
                    ### but devtools::install_github('kcha/polyester')
                    ### should be run first
                    )
