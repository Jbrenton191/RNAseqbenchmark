library(edgeR)
library(org.Hs.eg.db)
gene_list<-rownames(txi.salmon_sz02$counts)
gene_list<-as.data.frame(gene_list)
###need txi.salmon and conditon to be set from previous DESeq workflow
y <- DGEList(counts=txi.salmon$counts, group = condition, genes=rownames(txi.salmon$counts))
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)


idfound <- y$genes[,1] %in% mappedRkeys(org.Hs.egENSEMBL2EG)
y<-y[idfound,]
egENSEMBL2EG <- toTable(org.Hs.egENSEMBL2EG)
m <- match(y$genes[,1], egENSEMBL2EG$ensembl_id)
y$genes$EntrezGene <- egENSEMBL2EG$gene_id[m]
egSYMBOL <- toTable(org.Hs.egSYMBOL)

m <- match(y$genes$EntrezGene, egSYMBOL$gene_id)
y$genes$Symbol <- egSYMBOL$symbol[m]
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$Symbol)
y<-y[!d,]
nrow(y)
rownames(y)<-y$genes$Symbol

et <- exactTest(y)
topTags(et)
subset(et$table, PValue < 0.01 & abs(logFC) > 1)->EdgeR_sigGenes


setwd("C:/Users/JBrenton/Polyester_analysis_30_09_19")
read.table("Excitatory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Excitatory_genes
unique(Excitatory_genes$Gene_Name)->excit_gene_names
excit_gene_names<-as.vector(excit_gene_names)
read.table("Inhibitory_genes_transcripts_upregulated.csv", header = TRUE, sep=" ")->Inhibitory_genes
unique(Inhibitory_genes$Gene_Name)->inhib_gene_names
inhib_gene_names<-as.vector(inhib_gene_names)
all_changed_genes<-append(excit_gene_names, inhib_gene_names)

all_changed_genes[all_changed_genes %in% rownames(EdgeR_sigGenes)]
[1] "FAXC"   "MPLKIP" "SPPL2A" 

##rownames(EdgeR_sigGenes)[rownames(EdgeR_sigGenes) %in% all_changed_genes]
