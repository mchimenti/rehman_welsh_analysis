## Analysis of Rehman/Welsh human cultured epithelium normal and CF
## Date: 9.12.2018
## Author: Michael Chimenti
## Organism: hg38 / human 
## Aligners: hisat2 / salmon
## Design: CF/normal + cytokine treatments
## Reps: 6

##########
## Imports
##########

#source("https://bioconductor.org/biocLite.R")
#biocLite("DEGreport")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway and gene clusters
library('DEGreport')
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/welsh_lab/project_rehman_sept2018/") 

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 
#######################################
samples <- read.table("samples.csv", sep=',', header=TRUE)
rownames(samples) <- samples$sample
samples$group <- paste0(samples$condition, samples$treatment)
##samples 15 and 21 not found in the data
samples2 <- filter(samples, !sample %in% c("sample15","sample21"))

files <- file.path(getwd(), samples2$sample, 'salmon', 'quant.sf')
names(files) <- samples2$sample

tx2gene <- read_csv(file.path(getwd(), "tx2gene.csv"), col_names = FALSE)
tx2gene$X1 <- tx2gene$X1 %>%
  strsplit(split = '.', fixed = TRUE) %>%
  sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples2,
                                   design = ~ batch + group)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi <- DESeq(ddsTxi)

anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)

rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)
plotMA(object = ddsTxi, alpha = 0.05)
plotPCA(object = rldTxi, intgroup = 'group')

#looking at PCs 3 and 4
rld_mat <- assay(rldTxi)
pca <- prcomp(t(rld_mat))
df <- cbind(samples2, pca$x)
ggplot(df) + geom_point(aes(x=PC3,y=PC4, color = condition))

# drop samples?  

#ddsTxi <- ddsTxi[ , ddsTxi$sample != 's5']
#ddsTxi$sample <- droplevels(ddsTxi$sample)
#ddsTxi <- DESeq(ddsTxi)

##DE testing 
library("IHW")

###
res_wt_tnf <- results(ddsTxi, contrast = c("group","wttnf","wtnone"), filterFun = ihw)
res_wt_tnf <- na.omit(res_wt_tnf)  #drop NA rows
res_wt_tnf_sig <- res_wt_tnf[res_wt_tnf$padj < 0.05 & res_wt_tnf$baseMean > 5.0,]
res_wt_tnf_ord <- res_wt_tnf_sig[order(res_wt_tnf_sig$padj),]
res_wt_tnf_ord$ext_gene <- anno[row.names(res_wt_tnf_ord), "gene_name"]

png("volcano_nonCF_DEgenes_il17tnf_fdr_5percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_wt_tnf_ord, main = "Volcano: DE genes in nonCF cells w/ IL17 treatment", lfcthresh=1.5, sigthresh=0.05, textcx=.4, xlim=c(-4, 8), ylim = c(4,150))
dev.off()

degPlot(dds = ddsTxi, res = res_wt_tnf_ord, n = 9, xs = 'condition', group = 'treatment')
mycols <- c("baseMean", "log2FoldChange", "padj", "ext_gene")
write.csv(x = res_wt_tnf_ord[,mycols], file = "DEgenes_wt_tnf_vs_none_FDR_5percent.csv")

###
res_cf_tnf <- results(ddsTxi, contrast = c("group","cftnf","cfnone"), filterFun = ihw)
res_cf_tnf <- na.omit(res_cf_tnf)  #drop NA rows
res_cf_tnf_sig <- res_cf_tnf[res_cf_tnf$padj < 0.05 & res_cf_tnf$baseMean > 5.0,]
res_cf_tnf_ord <- res_cf_tnf_sig[order(res_cf_tnf_sig$padj),]
res_cf_tnf_ord$ext_gene <- anno[row.names(res_cf_tnf_ord), "gene_name"]

png("volcano_CF_DEgenes_il17tnf_FDR_5percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_cf_tnf_ord, main = "Volcano: DE genes in CF cells w/ IL17 treatment", lfcthresh=1.5, sigthresh=0.05, textcx=.4, xlim=c(-4, 8), ylim = c(4,100))
dev.off()

degPlot(dds = ddsTxi, res = res_cf_tnf_ord, n = 9, xs = "condition", group = "treatment")
write.csv(x = res_cf_tnf_ord[,mycols], file = "DEgenes_CF_tnf_vs_none_FDR5percent.csv")

###
res_cf_ncf <- results(ddsTxi, contrast = c("group","cfnone","wtnone"), filterFun = ihw)
res_cf_ncf <- na.omit(res_cf_ncf)  #drop NA rows
res_cf_ncf_sig <- res_cf_ncf[res_cf_ncf$padj < 0.05 & res_cf_ncf$baseMean > 5.0,]
res_cf_ncf_ord <- res_cf_ncf_sig[order(res_cf_ncf_sig$padj),]
res_cf_ncf_ord$ext_gene <- anno[row.names(res_cf_ncf_ord), "gene_name"]

png("volcano_CF_nCF_DEgenes_FDR_5percent.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_cf_ncf_ord, main = "Volcano: DE Genes btw CF and nonCF cells, FDR < 0.05", lfcthresh=1.5, sigthresh=0.05, textcx=.4, xlim=c(-6, 6), ylim = c(3,15))
dev.off()

degPlot(dds = ddsTxi, res = res_cf_ncf_ord, n = 9, xs = "condition", group = "treatment")
write.csv(x = res_cf_ncf_ord[,mycols], file = "DEgenes_noTNF_CF_vs_nonCF_FDR5percent.csv")

## get unique response of CF cells to treatment (i.e., genotype / treatment interaction)

ddsTxi2 <- DESeqDataSetFromTximport(txi,
                                   colData = samples2,
                                   design = ~ batch + condition + treatment + condition:treatment)

ddsTxi2 <- ddsTxi2[ rowSums(counts(ddsTxi2)) > 5, ]
ddsTxi2 <- DESeq(ddsTxi2)

ddsTxi3 <- DESeq(ddsTxi2, test = "LRT", reduced = ~ batch + condition + treatment)

res_lrt_int <- results(ddsTxi3)
res_lrt_int <- na.omit(res_lrt_int)  #drop NA rows
res_lrt_int_sig <- res_lrt_int[res_lrt_int$padj < 0.1 & res_lrt_int$baseMean > 5.0,]
res_lrt_int_ord <- res_lrt_int_sig[order(res_lrt_int_sig$padj),]
res_lrt_int_ord$ext_gene <- anno[row.names(res_lrt_int_ord), "gene_name"]

degPlot(dds = ddsTxi3, res = res_lrt_int_ord, n = 9, xs = "condition", group = "treatment")

png("test.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_lrt_int_ord, main = "Volcano:", lfcthresh=1.5, sigthresh=0.05, textcx=.4, xlim=c(-6, 6), ylim = c(3,15))
dev.off()
