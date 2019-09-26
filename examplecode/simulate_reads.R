setwd("/usr/share/sequencing/internships/JonathanBrenton/polyester_simulation/simulate")
library(doMC)
library(polyester)
library(Biostrings)
library(stringr)

### read the fasta file and imagine a coverage of 30 FPKM
### i.e. proportional to the transcripts length

writeLines("####### reading transcriptome FASTA")
fastaFile="/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/gencode.v32.transcripts.fa.gz"
fasta = readDNAStringSet(fastaFile)
fastaFile2="/usr/share/sequencing/references/rnaseq/Homo_sapiens/ENSEMBL/Homo_sapiens.GRCh38.cdna.all.fa"
fasta2 = readDNAStringSet(fastaFile2)

readspertx = round(30 * width(fasta) / 100)

## then imagine we have 2 groups of cases and controls
## in cases, I select 200 random genes involved with immunoglobulins
writeLines("##### preparing the differentially expressed genes")
genes<- names(fasta)
genesGRCh38<-names(fasta2)
#############changed immunoglobulin to 400 genes because many are lost
immunoglobulin <- sample(which(grepl("immunoglobulin", genesGRCh38)), 400)
immunoglobulin <- genesGRCh38[immunoglobulin]

## and in controls I select 150 random genes involved with potassium
potassium <- sample(which(grepl("potassium", genesGRCh38)), 150)
potassium <- genesGRCh38[potassium]

############################then convert GRCH38 genes to gencode-there is some loss between the two
transcriptnamesgenes<-str_extract_all(genes, "ENST.............")
transcriptnamespotassium<-str_extract_all(potassium, "ENST.............")
transcriptnamesimmunoglobulin<-str_extract_all(immunoglobulin, "ENST.............")

############added in some output tables to compare to later results and for reference

write.table(transcriptnamesimmunoglobulin, "immunoglobulintranscripts.txt", sep="\t", row.names=FALSE, col.names=FALSE)

write.table(transcriptnamespotassium, "potassiumtranscripts.txt", sep="\t", row.names=FALSE, col.names=FALSE)

write.table(transcriptnamesgenes, "potassiumtranscripts.txt", sep="\t", row.names=FALSE, col.names=FALSE)

####match to get positions in gencode fasta file
immunomatched<-match(transcriptnamesimmunoglobulin, transcriptnamesgenes)
immunopositionsgencode<-immunomatched[!is.na(immunomatched)]

######get and write out length to see how many genes are not converted
numberimmunogenes<-length(immunopositionsgencode)
numberimmunogenes<-as.character(numberimmunogenes)
writeLines("length of gencode immunoglobulin gene list")
writeLines(numberimmunogenes)

############same as above for potassium genes
potmatched<-match(transcriptnamespotassium, transcriptnamesgenes)
potassiumpositionsgencode<-potmatched[!is.na(potmatched)]
numberpotgenes<-length(potassiumpositionsgencode)
numberpotgenes<-as.character(numberpotgenes)
writeLines("length of gencode potassium gene list")
writeLines(numberpotgenes)

writeLines("##### preparing fold changes matrix")
## I then create a baseline fold matrix with 1
fold_changes <- matrix(
  rep(1,2*length(genes)),
  nrow = length(genes)
)

## and substitute in the fold matrix the selected genes
## in cases and controls
## with random fold changes between 2 and 4
fold_changes[immunopositionsgencode,1] <- sample(2:4,length(immunopositionsgencode), replace = T)
fold_changes[potassiumpositionsgencode,2] <- sample(2:4,length(potassiumpositionsgencode), replace = T)

writeLines("#### simulating the reads")
## at this point I can simulate the experiment
simulate_experiment(fasta = fastaFile,
                    outdir = "reads",
                    num_reps = c(20,20),
                    reads_per_transcript=readspertx,
                    fold_changes=fold_changes,
                    size = 0.2,
                    error_model = "illumina5",
                    gzip=TRUE,
                    cores=16
                    )
