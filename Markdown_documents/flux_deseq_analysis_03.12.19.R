library(EnsDb.Hsapiens.v86)
library(Biostrings)
library(stringr)
library(dplyr)
library(tximport)
library(plyr)
library(ggplot2)


#setwd("/home/AD/jbrenton/fastashuffle/output_15_10_19_sz03/alignments")
#####copied over files
##dir<-getwd()
###dir<-paste0(dir, "/salmon_quants_15_10_19_sz03")
###dir<-paste0(dir, "/salmon_quants_15_10_19_s02/alignments")
#setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s02/alignments")
#setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/salmon_quants_15_10_19_s03/alignments")

read.table("C:/Users/JBrenton/Polyester_analysis_30_09_19/gencode_v32_tx2gene2name.txt", header=FALSE, sep=" ")->extracted_gencode

colnames(extracted_gencode)<-c("tx_id", "gene_id","gene_name")
##remove version info
extracted_gencode$tx_id<-sub("\\..+", "", extracted_gencode$tx_id)
extracted_gencode$tx_id<-sub("^>", "", extracted_gencode$tx_id)

extracted_gencode$gene_id<-sub("\\..+", "", extracted_gencode$gene_id)


setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19/flux_salmon_quants_12_12_19/alignments/")

dir<-getwd()

c(seq(01,10))->x
samples<-c(1:length(x))
header<-c(1:length(x))
for (i in x) {
  headername<-paste0("sample_", i, "_shuffled.fasta_salmon_quant")
 
  ##samples[i] <- file.path(dir, headername, "quant.sf")
  if (file.exists(file.path(dir, headername, "quant.sf"))==TRUE) {
    samples[i] <- file.path(dir, headername, "quant.sf")
    header[i]<-paste0("sample_", i)
  }
  
}


names(samples) <- paste0(header)

#c(seq(1,9))->y
#samples2<-c(1:length(y))
#header2<-c(1:length(y))
#for (i in 1:length(y)) {
 # headername<-paste0("sample_0", i, "_salmon_quant")
  #header2[i]<-headername
  #if (file.exists(file.path(dir, headername, "quant.sf"))) {
#    samples2[i] <- file.path(dir, headername, "quant.sf")
#  }
 # }
#names(samples2) <- paste0(header2)

##samples[1:9]<-samples2[1:9]
##names(samples[1:9])<-names(samples2[1:9])

all(file.exists(samples))

print("############### got all folders, now importing")

txi.salmon_flux <- tximport(samples, type = "salmon", tx2gene = extracted_gencode, ignoreTxVersion=TRUE)

condition<-c(1:length(samples))
for (i in 1:length(samples)){
  if(i<6){
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
dds_flux <- DESeqDataSetFromTximport(txi.salmon_flux,
                                      colData = sampleMetadata,
                                      design = ~ condition)
keep <- rowSums(counts(dds_flux)) >= 10
dds_flux <- dds_flux[keep,]

dds_flux <- DESeq(dds_flux)
#res_flux <- results(dds_flux)
#res <- results(dds, name="condition_treated_vs_untreated")
res_flux <- results(dds_flux, contrast=c("condition","A","B"))

res_flux$gene <- row.names(res_flux)
resOrdered_flux <- res_flux[order(res_flux$pvalue),]
resOrdered_flux$gene <- row.names(resOrdered_flux)
resOrdered_flux <- as.data.frame(resOrdered_flux)
resOrdered_flux$colour <- "black"

immuno_genes<-read.table(file="C:/Users/JBrenton/Polyester_analysis_30_09_19/gene_list.txt", header = FALSE)
upreg_genenames<-unique(extracted_gencode$gene_id[extracted_gencode$gene_name %in% immuno_genes$V1])

resOrdered_flux$colour[which(resOrdered_flux$gene %in% upreg_genenames)] <- "red"


plot(res_flux$log2FoldChange, -log10(res_flux$padj), col='red', ylim=c(-0.1,10), xlim=c(-10,10), main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")
plot(resOrdered_flux$log2FoldChange, -log10(resOrdered_flux$padj), col=resOrdered_flux$colour, main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")
plot(log10(resOrdered_flux$baseMean), resOrdered_flux$log2FoldChange, col=resOrdered_flux$colour, main= "DESeq Results Plot: Excit_genes-Polyester Size factor=0.2")

setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
png("Flux_simulator_output.png")
ggplot(resOrdered_flux, aes(log2FoldChange, -log10(padj),))+
  geom_point(colour=resOrdered_flux$colour)
dev.off()
  
  ggplot(resOrdered_flux, aes(log10(baseMean), log2FoldChange))+
  geom_point(colour=resOrdered_flux$colour)+
  ggtitle("Flux Simulator MA plot: Log 10 Scale")+
  xlab("Mean reads per gene")+
  ylab("Fold Change: Log 2 Scale")+
  ylim(-15, 15)

  
  
  #theme(legend.position = "bottom", text = "red dots indicate fold change genes")
  
#  + labs(title =" New title", x = "New x", y = "New y") 
  
  

  
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(genefilter)
library(geneplotter)
library(topGO)
library(AnnotationDbi)

sigGenes <- rownames(subset(res_flux, padj < 0.001 & abs(log2FoldChange > 2)))

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              keys=rownames(res_flux), 
              columns=c("SYMBOL","GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

sample_n(anSig, 5)



#####make a gene list again for some reason?? now make a background gene list that have a similar mean expression 
###level to DE genes just they aren't differentially expressed, but why not whole universe???

overallBaseMean <- as.matrix(res_flux[, "baseMean", drop = F])

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


############################ doing topGO analysis

multidensity( list( 
  all= log2(res_flux[,"baseMean"]) ,
  foreground =log2(res_flux[anSig$ENSEMBL, "baseMean"]), 
  background =log2(res_flux[backG, "baseMean"])), 
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


####extra go working out-to get pasted go term as variable and have attached names and IDs
#assign(paste(go_id), AnnotationDbi::select(org.Hs.eg.db, 
#                                           keys=go_id, 
#                                           columns=c("SYMBOL","GENENAME","ENSEMBL"),
#                                           keytype="GO"))
#get(paste(go_id))
