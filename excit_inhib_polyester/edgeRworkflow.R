####EdgeR workflow

###used counts imported from txi.salmon from DESeq2 workflow
library(edgeR)
library(org.Hs.eg.db)
gene_list<-rownames(txi.salmon_sz02$counts)
gene_list<-as.data.frame(gene_list)
y <- DGEList(counts=txi.salmon_sz02$counts, group = condition, genes=rownames(txi.salmon_sz02$counts))
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
topTags(et)
plotMDS(y)


mappedkeys(org.Hs.egENSEMBL2EG)
###########y$genes doesn't work below because it's a dataframe and not a character vector
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
####removes lowly expressed genes:below
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
design<-model.matrix(~condition)
#y <- estimateDisp(y, design = condition)

library(statmod)
y <- estimateDisp(y, design = design, robust = TRUE)
plotBCV(y)
fit <- glmFit(y, design)

