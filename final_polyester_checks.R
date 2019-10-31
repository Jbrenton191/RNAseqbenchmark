#######checking for consistency
###chagned gencode to extracted gencode from original for downstream analysis-the original genecode to genes
###had some genes removed for some reason.
###made no apparent difference
###also plugged in the positions used in the fold change matrix into the results matrices from DESeq
###and got a match for those I did with normal analysis
#####sz02 protocol
library(EnsDb.Hsapiens.v86)
library(Biostrings)
library(stringr)
library(dplyr)
library(tximport)
library(plyr)

#setwd("/home/AD/jbrenton/fastashuffle/output_15_10_19_sz03/alignments")
#####copied over files
#setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
##dir<-getwd()
###dir<-paste0(dir, "/salmon_quants_15_10_19_sz03")
###dir<-paste0(dir, "/salmon_quants_15_10_19_s02/alignments")
setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s02/alignments")
#setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s03/alignments")
dir<-getwd()

c(seq(01,40))->x
samples<-c(1:length(x))
header<-c(1:length(x))
for (i in x) {
  headername<-paste0("sample_", i, "_salmon_quant")
  header[i]<-headername
  #samples[i] <- file.path(dir, headername, "quant.sf")
  if (file.exists(file.path(dir, headername, "quant.sf"))==TRUE) {
    samples[i] <- file.path(dir, headername, "quant.sf")
  }
  
}

names(samples) <- paste0(header)

c(seq(1,9))->y
samples2<-c(1:length(y))
header2<-c(1:length(y))
for (i in 1:length(y)) {
  headername<-paste0("sample_0", i, "_salmon_quant")
  header2[i]<-headername
  if (file.exists(file.path(dir, headername, "quant.sf"))) {
    samples2[i] <- file.path(dir, headername, "quant.sf")
  }
}
names(samples2) <- paste0(header2)

samples[1:9]<-samples2[1:9]
names(samples[1:9])<-names(samples2[1:9])

###failed simulation fasta numbers for sz03
#y<-samples[c(21,27,02,28,34)]

###failed simulation fasta numbers for sz02
y<-samples[c(7,14,3,32,35,39)]
samples<-setdiff(samples,y)
###for some reason removes all names as well
#y<-c(21,27,02,28,34)
y<-c(7,14,3,32,35,39)
header<-setdiff(header,header[y])
names(samples) <- header
#names(samples[1:9]) <- paste0(header2[1:9])


all(file.exists(samples))

print("############### got all folders, now importing")
####removed gene version with sed in bash to make new file
###genecode_txid_to_geneid<-read.delim(file ="/usr/share/sequencing/internships/JonathanBrenton/polyester_simulation/simulate/polyester_tests_13.10.19_excit+inhib_genes/genecode_txid_to_geneid.txt",sep=" ", header=FALSE)
read.csv("extracted_gencode.csv")->extracted_gencode
#colnames(extracted_gencode)<-c("tx_id", "gene_id","gene_name")
#genecode_txid_to_geneid$tx_id<-sub("\\..+", "", genecode_txid_to_geneid$tx_id)
#genecode_txid_to_geneid$gene_id<-sub("\\..+", "", genecode_txid_to_geneid$gene_id)

#> length(unique(extracted_gencode$Gene_Name))
#[1] 59368

txi.salmon_sz02 <- tximport(samples, type = "salmon", tx2gene = extracted_gencode, ignoreTxVersion=TRUE)

#reading in files with read_tsv
#1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 
#removing duplicated transcript rows from tx2gene
#transcripts missing from tx2gene: 53
#summarizing abundance
#summarizing counts
#summarizing length

####1000 genes went where/missing?
#dim(res_sz02)
#[1] 58381     7

condition<-c(1:length(samples))
for (i in 1:length(samples)){
  if(i<20){
    condition[i]<-"A"
  } else{
    condition[i]<-"B"
  }
}

library(DESeq2)
sampleMetadata <- data.frame(
  samples = names(samples),
  condition = condition
)
dds_sz02 <- DESeqDataSetFromTximport(txi.salmon_sz02,
                                     colData = sampleMetadata,
                                     design = ~ condition)
keep <- rowSums(counts(dds_sz02)) >= 10
dds_sz02 <- dds_sz02[keep,]

dds_sz02 <- DESeq(dds_sz02)
res_sz02 <- results(dds_sz02)
res_sz02$gene <- row.names(res_sz02)
resOrdered_sz02 <- res_sz02[order(res_sz02$pvalue),]
resOrdered_sz02$gene <- row.names(resOrdered_sz02)
resOrdered_sz02 <- as.data.frame(resOrdered_sz02)

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(genefilter)
library(geneplotter)
library(topGO)


sigGenes <- rownames(subset(res_sz02, padj < 0.05 & abs(log2FoldChange > 1)))

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys=rownames(res_sz02), 
                              columns=c("GENENAME"),
                              keytype="SYMBOL")

anSig <- as.data.frame(subset(anno, SYMBOL %in% sigGenes))

sample_n(anSig, 5)

#####make a gene list again for some reason?? now make a background gene list that have a similar mean expression 
###level to DE genes just they aren't differentially expressed, but why not whole universe???

overallBaseMean <- as.matrix(res_sz02[, "baseMean", drop = F])

sig_idx <- match(anSig$SYMBOL, rownames(overallBaseMean))

backG <- c()



for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

backG <- setdiff(backG,  anSig$SYMBOL)
length(backG)

##Plotting the density of the average expressions, shows that the background matching has worked reasonably well.
res_sz02<-na.omit(res_sz02)

multidensity( list( 
  all= log2(res_sz02[,"baseMean"]) ,
  foreground =log2(res_sz02[anSig$SYMBOL, "baseMean"]), 
  background =log2(res_sz02[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

onts = c( "MF", "BP", "CC" )

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(anSig$SYMBOL,  backG) 
inSelection =  geneIDs %in% anSig$SYMBOL 
####don't really understand how there two logical vectors interact!
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tab = as.list(onts)
names(tab) = onts
for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "symbol" )
  
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim, 
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = length(tgd@graph@nodes))
  
}

seq<-c()
gotest<-c()
for (i in 1:3){
  seq<-rep(names(tab)[i],dim(tab[[i]][1])[1])
  gotest<-append(gotest,seq)
}

gotest<-as.data.frame(gotest)
topGOResults <- rbind.fill(tab)
topGOResults <- cbind(topGOResults, gotest)
#topGOResults<-left_join(topGOResults,gotest)
write.csv(topGOResults, file = "topGOResults.csv")

###################graphs

setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
####for sz 02
unique(Excitatory_genes$Gene_Name)->excit_gene_names
excit_gene_names<-as.data.frame(excit_gene_names)
names(resOrdered_sz02)[7]<-"Gene_Name"
dplyr::left_join(resOrdered_sz02, extracted_gencode, by='Gene_Name')->linked_sz02
subset(linked_sz02, padj < 0.05 & abs(log2FoldChange > 1))->sigOrderedgenes_sz02
#names(Excitatory_genes)[1]<-"tx_id"
Excitatory_genes<-dplyr::left_join(Excitatory_genes, extracted_gencode, by='Gene_Name')
names(excit_gene_names)[1]<-"Gene_Name"
sigOrderedgenes_sz02$gene_name %in% excit_gene_names$Gene_Name
#right_join(sigOrderedgenes,Excitatory_genes, by="gene_name")->test
names(sigOrderedgenes_sz02)[names(sigOrderedgenes_sz02)=="gene_name"]<-"Gene_Name"
right_join(sigOrderedgenes_sz02,excit_gene_names, by="Gene_Name")->test_sz02

#names(sigOrderedgenes)[names(sigOrderedgenes)=="gene_name"]<-"Gene_Name"
res_sz02$colour="black"
#read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
#names(Excitatory_genes)[1]<-"tx_id"
right_join(extracted_gencode, Excitatory_genes, by="Gene_Name")->excit_merged
#unique(na.omit(excit_merged))
na.omit(excit_merged)->excit_merged_no_na
unique(excit_merged_no_na$Gene_Name)->excit_gene_names
excit_matchvector<-na.omit(match(excit_gene_names,res_sz02$gene))
####NAs found in this
res_sz02$colour[excit_matchvector]="red"
#png(filename = "Full_DESeq Results Plot-Polyester Size factor=0.3.png", width = 1024, height = 768)
png(filename = "Full_DESeq Results Plot-Polyester Size factor=0.2.png", width = 1024, height = 768)
plot(res_sz02$log2FoldChange, -log10(res_sz02$padj), col=res_sz02$colour, ylim=c(-0.1,10), xlim=c(-10,10), main= "Full DESeq Results Plot: Polyester Size factor=0.2")
#plot(res$log2FoldChange, -log10(res$padj), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10), main= "Full DESeq Results Plot: Polyester Size factor=0.3")
dev.off()
x<-recordPlot()
#plot(res$log2FoldChange[excit_matchvector], -log10(res$padj[excit_matchvector]), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10))

png(filename = "Full_DESeq Results Plot-Excit_genes-Polyester Size factor=0.3.png", width = 1024, height = 768)
plot(res_sz02$log2FoldChange[excit_matchvector], -log10(res_sz02$padj[excit_matchvector]), col='red', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")
#png(filename = "Full_DESeq Results Plot-Excit_genes-Polyester Size factor=0.2.png", width = 1024, height = 768)
#plot(res$log2FoldChange[excit_matchvector], -log10(res$padj[excit_matchvector]), col='red', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")

y<-recordPlot()
dev.off()

###inhib

read.table("Inhibitory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Inhibitory_genes
unique(Inhibitory_genes$Gene_Name)->inhib_gene_names
inhib_gene_names<-as.data.frame(inhib_gene_names)
resOrdered_sz02<-as.data.frame(resOrdered_sz02)
dplyr::left_join(resOrdered_sz02, genecode_txid_to_geneid, by='gene_id')->linked_sz02

names(sigOrderedgenes_sz02)[9]<-"Gene_Name"

subset(linked_sz02, padj < 0.01 & abs(log2FoldChange > 1))->sigOrderedgenessz02
sigOrderedgenes_sz02$gene_name %in% inhib_gene_names   #####all false
right_join(sigOrderedgenes_sz02,Inhibitory_genes, by="Gene_Name")->test_sz02_inhib
na.omit(sigOrderedgenes_sz02) #####leaves 0 genes

names(Inhibitory_genes)[1]<-"tx_id"
left_join(genecode_txid_to_geneid, Inhibitory_genes, by="tx_id")->inhib_merged
na.omit(inhib_merged)->inhib_merged_no_na
unique(inhib_merged_no_na$gene_id)->inhib_gene_names
inhib_matchvector<-na.omit(match(inhib_gene_names,res_sz02$gene))
res_sz02$colour[inhib_matchvector]<-"purple"
#plot(res$log2FoldChange, -log10(res$padj), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10))
#png(filename = "Full_DESeq Results Plot-Inhib_genes-Polyester Size factor=0.3.png", width = 1024, height = 768)
png(filename = "Full_DESeq Results Plot-Inhib_genes-Polyester Size factor=0.2.png", width = 1024, height = 768)

plot(res_sz02$log2FoldChange[inhib_matchvector], -log10(res_sz02$padj[inhib_matchvector]), col='purple', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Inhib_genes-Polyester Size factor=0.2")
dev.off()
#z<-recordPlot()

###further checks

sigGenes <- rownames(subset(res, pvalue < 0.01 & abs(log2FoldChange > 1)))

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys=rownames(res), 
                              columns=c("SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

sample_n(anSig, 5)

anSig<-na.omit(anSig)


setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
unique(Excitatory_genes$Gene_Name)->excit_gene_names
excit_gene_names<-as.vector(excit_gene_names)
read.table("Inhibitory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Inhibitory_genes
unique(Inhibitory_genes$Gene_Name)->inhib_gene_names
inhib_gene_names<-as.vector(inhib_gene_names)
all_changed_genes<-append(excit_gene_names, inhib_gene_names)

all_changed_genes[all_changed_genes %in% anSig$SYMBOL]
####[1] "NRROS"

######final checks
####made csv of positions used in fold change matrices, copied from server to local
setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/")
read.csv("excit_positions.csv", header=FALSE)->excit_positions
read.csv("inhib_positions.csv", header=FALSE)->inhib_positions
excit_pos_in_gencode<-extracted_gencode[excit_positions[,1],]
inhib_pos_in_gencode<-extracted_gencode[inhib_positions[,1],]

factor(Excitatory_genes$Gene_Name %in% excit_pos_in_gencode$Gene_Name)
###Levels: TRUE

factor(Inhibitory_genes$Gene_Name %in% inhib_pos_in_gencode$Gene_Name)
##Levels: TRUE

###checking the two are exactly the same
#sort(excit_gene_names, decreasing=TRUE)->t1
#sort(unique(excit_pos_in_gencode$Gene_Name), decreasing=TRUE)->t2
#as.vector(t1)->t1
#as.vector(t2)->t2
#identical(t1,t2)
# [1] TRUE

###checking that these genes exist in the results file
factor(excit_pos_in_gencode$Gene_Name %in% resOrdered_sz02$gene)->t1
count(t1)
#x freq
#1 FALSE    4
#2  TRUE  871

#number of inhib genes present
inhib_pos_in_gencode$Gene_Name %in% resOrdered_sz02$gene->i1
count(i1)
##x freq
##1 FALSE    8
##2  TRUE 3473

unique(excit_pos_in_gencode$Gene_Name) %in% resOrdered_sz02$gene->e1
unique(excit_pos_in_gencode$Gene_Name) %in% resOrdered_sz02$gene->e1
#> e1
#[1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[15]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[29]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[43]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE
#[57]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[71]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[85]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[99]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[113]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#[127]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
#> unique(excit_pos_in_gencode$Gene_Name)[47]
#[1] DUX4L5
#59368 Levels: 5_8S_rRNA 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZEF1
#> unique(excit_pos_in_gencode$Gene_Name)[49]
#[1] DUX4L1
#59368 Levels: 5_8S_rRNA 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZEF1
#> unique(excit_pos_in_gencode$Gene_Name)[50]
#[1] DUX4L3
#59368 Levels: 5_8S_rRNA 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZEF1
#> unique(excit_pos_in_gencode$Gene_Name)[51]
#[1] DUX4L2
#59368 Levels: 5_8S_rRNA 5S_rRNA 7SK A1BG A1BG-AS1 A1CF A2M A2M-AS1 A2ML1 ... ZZEF1

####DUX4L4 works fine

#unique(inhib_pos_in_gencode$Gene_Name) %in% resOrdered_sz02$gene->i1
#i1
#> which(i1==FALSE)
#[1] 83 84 86 87 88 89

#unique(inhib_pos_in_gencode$Gene_Name)[which(i1==FALSE)]
#[1] USP17L25 USP17L26 USP17L27 USP17L28 USP17L29 USP17L30
