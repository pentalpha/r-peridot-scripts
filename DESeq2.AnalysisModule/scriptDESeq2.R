args = commandArgs(trailingOnly = F)

localDir <- args[length(args)-3]

localDir

inputFilesDir <- args[length(args)-2]

inputFilesDir

outputFilesDir <- args[length(args)-1]

outputFilesDir

notFirstRun <- args[length(args)]

notFirstRun

localDir
setwd(localDir)
options(bitmapType='cairo')
#Get file config
FileConfigPath = paste(localDir, "config.txt", sep = "/")

#Read file config
FileConfig = read.table(FileConfigPath, header = TRUE, row.names = 1, sep = "|")

if(notFirstRun == "0"){
  peridotConditions = read.table(paste(inputFilesDir, "condition-input.tsv", sep = "/"), header=TRUE, row.names=1)
  peridotConditions

  #Read Path file
  inputTableFile = paste(inputFilesDir, "rna-seq-input.tsv", sep = "/")

  #Read file
  peridotCountTable = read.table(inputTableFile, header=TRUE, row.names=1 )

  #Ignore samples with "not-use" indicated
  #first, remove they from the conditions table
  peridotConditions <- subset(peridotConditions, condition != "not-use")
  #then, remove from the counts table
  for(i in colnames(peridotCountTable)){
    iContainsNotUse = length(grep("not.use", as.name(i))) > 0
    if(iContainsNotUse){
      #erases the column
      peridotCountTable[, i] = NULL
    }
  }
  #Finally, drop unused levels (not-use levels)
  peridotConditions = droplevels(peridotConditions)
  peridotConditions

  countData <- data.matrix(peridotCountTable)

  head(countData)

  condFac <- as.factor(peridotConditions$condition)

  condFac

  library(DESeq2)

  dds = DESeqDataSetFromMatrix(countData = countData, DataFrame(condFac), ~condFac)

  dds

  dds <- DESeq(dds)

  res <- results(dds)

  res
}else{
  load(file = "DESeq2.RData")
}

adaptativeRowMeans <- function(columns){
  if(is.null(dim(columns))){
    return(columns)
  }else{
    return(rowMeans(columns))
  }
}

if(notFirstRun == "0"){
  resFinal <- data.frame(baseMean = res$baseMean, row.names =  rownames(res))

  baseMeanPerLvl <- sapply( levels(condFac), function(lvl) adaptativeRowMeans( counts(dds,normalized=TRUE)[,condFac == lvl] ) )

  colnames(baseMeanPerLvl) <- c("baseMeanB", "baseMeanA")

  resFinal$baseMeanA <- baseMeanPerLvl[,1]

  resFinal$baseMeanB <- baseMeanPerLvl[,2]

  resFinal$foldChange <- 2^res$log2FoldChange

  resFinal$log2FoldChange <- res$log2FoldChange

  resFinal$pvalue <- res$pvalue

  resFinal$padj <- res$padj
}

head(resFinal)

#Histogram PValue and FDR Function
peridotPlotHist <- function(res){
  p1 <- with(resFinal, hist(pvalue, breaks=100, plot = F))
  p2 <- with(resFinal, hist(padj, breaks=100, plot = F))
  plot( p1, col="skyblue", main="Histogram", xlab = "Values")  # first histogram
  plot( p2, col=scales::alpha('red',.5), add=T)
  legend('topleft', c("PValue", "FDR(padjust)"), fill = c("skyblue", scales::alpha('red',.5)), bty = 'o', border = NA, cex = 0.8, bg = "white")
}

#MA Plot Function
peridotPlotMA <- function(res, config){
  with(resFinal, plot(log(baseMean), log2FoldChange, pch=20, main="MA Plot"))

  if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr & pvalue<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
    with(subset(resFinal, pvalue<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR(padj) < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr & pvalue<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0){
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(resFinal, pvalue<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("FDR(padj) < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }  
}

#Volcano Plot Function
peridotPlotVolcano <- function(res, config){
  with(resFinal, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
  
  if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr & pvalue<FileConfig$pValue), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
    with(subset(resFinal, pvalue<FileConfig$pValue), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR(padj) < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr & pvalue<FileConfig$pValue), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    legend('topleft', paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0){
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(resFinal, pvalue<FileConfig$pValue), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    legend('topleft', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(resFinal, padj<FileConfig$fdr), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
    legend('topleft', paste("FDR(padj) < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }
}

png(filename = paste(outputFilesDir, "2-histogram.png", sep = "/"), width=600, height=600)

peridotPlotHist(resFinal)

dev.off()

png(filename = paste(outputFilesDir, "3-MAPlot.png", sep = "/"), width=600, height=600)

peridotPlotMA(resFinal, FileConfig)

dev.off()

png(filename = paste(outputFilesDir, "4-volcanoPlot.png", sep = "/"), width=600, height=600)

peridotPlotVolcano(resFinal, FileConfig)

dev.off()

#Subset com PValue < FileConfig$pValue, Fold Change < FileConfig$log2FoldChange e FDR < FileConfig$fdr
if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & pvalue <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & pvalue <FileConfig$pValue))
}else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & padj < FileConfig$fdr))
}else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(pvalue <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$FoldChange != 0){
  resSig = subset(resFinal, abs(log2FoldChange) > log2(FileConfig$FoldChange))
}else if(FileConfig$pValue != 0){
  resSig = subset(resFinal, pvalue <FileConfig$pValue)
}else if(FileConfig$fdr != 0){
  resSig = subset(resFinal, padj < FileConfig$fdr)
}else{
  resSig = resFinal
}

##Remove NA
resSig = na.omit(resSig)

if(FileConfig$tops > 0){
  topRes = head(resSig, n = FileConfig$tops)

  write.table(topRes, paste(outputFilesDir, "/TopResults.tsv", sep = ""), sep = "\t")
}

##Create files csv

write.table(resSig, paste(outputFilesDir, "/1-res.tsv", sep = ""), sep = "\t")


pdf(file = paste(outputFilesDir, "5-plots.pdf", sep = "/"))

#Histogram
peridotPlotHist(resFinal)

#MA Plot
peridotPlotMA(resFinal, FileConfig)

#Volcano Plot
peridotPlotVolcano(resFinal, FileConfig)

dev.off()

save(resFinal, file = "DESeq2.RData")
