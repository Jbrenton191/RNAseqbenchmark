#####sz03 protocol
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
#setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s02/alignments")
setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s03/alignments")
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
y<-samples[c(21,27,02,28,34)]

###failed simulation fasta numbers for sz02
#y<-samples[c(7,14,3,32,35,39)]
samples<-setdiff(samples,y)
###for some reason removes all names as well
y<-c(21,27,02,28,34)
#y<-c(7,14,3,32,35,39)
header<-setdiff(header,header[y])
names(samples) <- header
#names(samples[1:9]) <- paste0(header2[1:9])


all(file.exists(samples))

print("############### got all folders, now importing")
####removed gene version with sed in bash to make new file
###genecode_txid_to_geneid<-read.delim(file ="/usr/share/sequencing/internships/JonathanBrenton/polyester_simulation/simulate/polyester_tests_13.10.19_excit+inhib_genes/genecode_txid_to_geneid.txt",sep=" ", header=FALSE)
genecode_txid_to_geneid<-read.delim(file ="C:/Users/JBrenton/Polyester_analysis_30_09_19/deseq_work_06.10.19/genecode_txid_to_geneid.txt", sep=" ", header=FALSE)
colnames(genecode_txid_to_geneid)<-c("tx_id", "gene_id","gene_name")
genecode_txid_to_geneid$tx_id<-sub("\\..+", "", genecode_txid_to_geneid$tx_id)
genecode_txid_to_geneid$gene_id<-sub("\\..+", "", genecode_txid_to_geneid$gene_id)

txi.salmon <- tximport(samples, type = "salmon", tx2gene = genecode_txid_to_geneid, ignoreTxVersion=TRUE)

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
dds <- DESeqDataSetFromTximport(txi.salmon,
                                colData = sampleMetadata,
                                design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds)
res$gene <- row.names(res)
resOrdered <- res[order(res$pvalue),]
resOrdered$gene <- row.names(resOrdered)
resOrdered <- as.data.frame(resOrdered)

###########graphing

setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
####for sz 03
unique(Excitatory_genes$Gene_Name)->excit_gene_names
excit_gene_names<-as.data.frame(excit_gene_names)
names(resOrdered)[7]<-"gene_id"
dplyr::left_join(resOrdered, genecode_txid_to_geneid, by='gene_id')->linked
subset(linked, padj < 0.01 & abs(log2FoldChange > 1))->sigOrderedgenes
names(Excitatory_genes)[1]<-"tx_id"
Excitatory_genes<-dplyr::left_join(Excitatory_genes, genecode_txid_to_geneid, by='tx_id')
names(excit_gene_names)[1]<-"Gene_Name"
sigOrderedgenes$gene_name %in% excit_gene_names$Gene_Name
#right_join(sigOrderedgenes,Excitatory_genes, by="gene_name")->test
names(sigOrderedgenes)[names(sigOrderedgenes)=="gene_name"]<-"Gene_Name"
right_join(sigOrderedgenes,excit_gene_names, by="Gene_Name")->test

#names(sigOrderedgenes)[names(sigOrderedgenes)=="gene_name"]<-"Gene_Name"
res$colour="black"
#read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
#names(Excitatory_genes)[1]<-"tx_id"
left_join(genecode_txid_to_geneid, Excitatory_genes, by="tx_id")->excit_merged
#unique(na.omit(excit_merged))
na.omit(excit_merged)->excit_merged_no_na
unique(excit_merged_no_na$gene_id.x)->excit_gene_names
excit_matchvector<-na.omit(match(excit_gene_names,res$gene))
res$colour[excit_matchvector]="red"
png(filename = "Full_DESeq Results Plot-Polyester Size factor=0.3.png", width = 1024, height = 768)
#png(filename = "Full_DESeq Results Plot-Polyester Size factor=0.2.png", width = 1024, height = 768)
#plot(res$log2FoldChange, -log10(res$padj), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10), main= "Full DESeq Results Plot: Polyester Size factor=0.2")
plot(res$log2FoldChange, -log10(res$padj), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10), main= "Full DESeq Results Plot: Polyester Size factor=0.3")
dev.off()
x<-recordPlot()
#plot(res$log2FoldChange[excit_matchvector], -log10(res$padj[excit_matchvector]), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10))

png(filename = "Full_DESeq Results Plot-Excit_genes-Polyester Size factor=0.3.png", width = 1024, height = 768)
plot(res$log2FoldChange[excit_matchvector], -log10(res$padj[excit_matchvector]), col='red', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.3")
#png(filename = "Full_DESeq Results Plot-Excit_genes-Polyester Size factor=0.2.png", width = 1024, height = 768)
#plot(res$log2FoldChange[excit_matchvector], -log10(res$padj[excit_matchvector]), col='red', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")

y<-recordPlot()
dev.off()

###inhib

read.table("Inhibitory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Inhibitory_genes
unique(Inhibitory_genes$Gene_Name)->inhib_gene_names
inhib_gene_names<-as.data.frame(inhib_gene_names)
resOrdered<-as.data.frame(resOrdered)
dplyr::left_join(resOrdered, genecode_txid_to_geneid, by='gene_id')->linked

names(sigOrderedgenes)[9]<-"Gene_Name"

subset(linked, padj < 0.01 & abs(log2FoldChange > 1))->sigOrderedgenes
sigOrderedgenes$gene_name %in% inhib_gene_names   #####all false
right_join(sigOrderedgenes,Inhibitory_genes, by="Gene_Name")->test
na.omit(test) #####leaves 0 genes

names(Inhibitory_genes)[1]<-"tx_id"
left_join(genecode_txid_to_geneid, Inhibitory_genes, by="tx_id")->inhib_merged
na.omit(inhib_merged)->inhib_merged_no_na
unique(inhib_merged_no_na$gene_id)->inhib_gene_names
inhib_matchvector<-na.omit(match(inhib_gene_names,res$gene))
res$colour[inhib_matchvector]<-"purple"
#plot(res$log2FoldChange, -log10(res$padj), col=res$colour, ylim=c(-0.1,10), xlim=c(-10,10))
png(filename = "Full_DESeq Results Plot-Inhib_genes-Polyester Size factor=0.3.png", width = 1024, height = 768)
#png(filename = "Full_DESeq Results Plot-Inhib_genes-Polyester Size factor=0.2.png", width = 1024, height = 768)

plot(res$log2FoldChange[inhib_matchvector], -log10(res$padj[inhib_matchvector]), col='purple', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Inhib_genes-Polyester Size factor=0.3")
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(genefilter)
library(geneplotter)
library(topGO)


sigGenes <- rownames(subset(res, padj < 0.01 & abs(log2FoldChange > 1)))

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys=rownames(res), 
                              columns=c("SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

sample_n(anSig, 5)

#####make a gene list again for some reason?? now make a background gene list that have a similar mean expression 
###level to DE genes just they aren't differentially expressed, but why not whole universe???

overallBaseMean <- as.matrix(res[, "baseMean", drop = F])

sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))

backG <- c()



for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

backG <- setdiff(backG,  anSig$ENSEMBL)
length(backG)

##Plotting the density of the average expressions, shows that the background matching has worked reasonably well.
res<-na.omit(res)

multidensity( list( 
  all= log2(res[,"baseMean"]) ,
  foreground =log2(res[anSig$ENSEMBL, "baseMean"]), 
  background =log2(res[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis")

onts = c( "MF", "BP", "CC" )

geneIDs = rownames(overallBaseMean)
inUniverse = geneIDs %in% c(anSig$ENSEMBL,  backG) 
inSelection =  geneIDs %in% anSig$ENSEMBL 
####don't really understand how there two logical vectors interact!
alg <- factor( as.integer( inSelection[inUniverse] ) )
names(alg) <- geneIDs[inUniverse]

tab = as.list(onts)
names(tab) = onts
for(i in 1:3){
  
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  
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




######Further checks

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
##[1] "NRROS"




