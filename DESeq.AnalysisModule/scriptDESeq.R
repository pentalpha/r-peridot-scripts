args = commandArgs(trailingOnly = F)

localDir <- args[length(args)-3]

localDir

inputFilesDir <- args[length(args)-2]

inputFilesDir

outputFilesDir <- args[length(args)-1]

outputFilesDir

notFirstRun <- args[length(args)]

notFirstRun

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

  library(DESeq)

  peridotCountTable[] = lapply(peridotCountTable, function(x) as.integer(x))

  #Create count data set
  cds = newCountDataSet(peridotCountTable, peridotConditions$condition)

  #Estimate the size factors for a count data set
  cds = estimateSizeFactors(cds)

  #The sizeFactors vector assigns to each column of the count data a value, the size factor
  sizeFactors(cds)

  #The counts slot holds the count data as a matrix of non-negative integer count values
  head(counts(cds, normalized=TRUE))

  #Obtains dispersion estimates for a count data set
  cds = estimateDispersions( cds , fitType = "local")

  #Compactly display the internal structure
  str( fitInfo(cds) )

  head( fData(cds) )

  #Tests for differences between the base means of two conditions (i.e., for differential expression in the case of RNA-Seq)
  res = nbinomTest( cds, levels(peridotConditions$condition)[2], levels(peridotConditions$condition)[1] )

  head(res)

}else{
  load(file = "DESeq.RData")
}

#Histogram PValue and FDR Function
peridotPlotHist <- function(res){
  p1 <- with(res, hist(pval, breaks=100, plot = F))
  p2 <- with(res, hist(padj, breaks=100, plot = F))
  plot( p1, col="skyblue", main="Histogram", xlab = "Values")  # first histogram
  plot( p2, col=scales::alpha('red',.5), add=T)
  legend('topleft', c("PValue", "FDR(padjust)"), fill = c("skyblue", scales::alpha('red',.5)), bty = 'o', border = NA, cex = 0.8, bg = "white")
}

#MA Plot Function
peridotPlotMA <- function(res, config){
  with(res, plot(log(baseMean), log2FoldChange, pch=20, main="MA Plot"))
  
  if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr & pval<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
    with(subset(res, pval<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR(padj) < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr & pval<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0){
    abline(h=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('bottomright', paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(res, pval<FileConfig$pValue), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr), points(log(baseMean), log2FoldChange, pch=20, col="red"))
    legend('bottomright', paste("FDR(padj) < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }  
}

#Volcano Plot Function
peridotPlotVolcano <- function(res, config){
  with(res, plot(log2FoldChange, -log10(pval), pch=20, main="Volcano plot"))
  
  if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr & pval<FileConfig$pValue), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
    with(subset(res, pval<FileConfig$pValue), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR(padj) < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr & pval<FileConfig$pValue), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    legend('topleft', paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$FoldChange != 0){
    abline(v=c(log2(FileConfig$FoldChange), log2(FileConfig$FoldChange)*(-1)), col="blue")
    legend('topleft', paste("log2FoldChange = mod(", round(log2(FileConfig$FoldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(res, pval<FileConfig$pValue), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    legend('topleft', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(res, padj<FileConfig$fdr), points(log2FoldChange, -log10(pval), pch=20, col="red"))
    legend('topleft', paste("FDR(padj) < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }
}

#Call Histogram Plot
png(filename = paste(outputFilesDir, "histogram.png", sep = "/"), width=600, height=600)

peridotPlotHist(res)

dev.off()

#Call MA Plot
png(filename = paste(outputFilesDir, "MAPlot.png", sep = "/"), width=600, height=600)

peridotPlotMA(res, FileConfig)

dev.off()

#Call Volcano Plot
png(filename = paste(outputFilesDir, "volcanoPlot.png", sep = "/"), width=600, height=600)

peridotPlotVolcano(res, FileConfig)

dev.off()

resAux = res

rownames(resAux) = resAux$id
resAux$id = NULL

#Subset com PValue < FileConfig$pValue, Fold Change < FileConfig$log2FoldChange e FDR < FileConfig$fdr
if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & pval <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$FoldChange != 0 & FileConfig$pValue != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & pval <FileConfig$pValue))
}else if(FileConfig$FoldChange != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$FoldChange) & padj < FileConfig$fdr))
}else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(pval <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$FoldChange != 0){
  resSig = subset(resAux, abs(log2FoldChange) > log2(FileConfig$log2FoldChange))
}else if(FileConfig$pValue != 0){
  resSig = subset(resAux, pval <FileConfig$pValue)
}else if(FileConfig$fdr != 0){
  resSig = subset(resAux, padj < FileConfig$fdr)
}else{
  resSig = resAux
}

##Remove NA
resSig = na.omit(resSig)

if(FileConfig$tops > 0){
  topRes = head(resSig, n = FileConfig$tops)

  write.table(topRes, paste(outputFilesDir, "/TopResults.tsv", sep = ""), sep = "\t")
}

##Create files csv
write.table(resSig, paste(outputFilesDir, "/res.tsv", sep = ""), sep = "\t")

pdf(file = paste(outputFilesDir, "plots.pdf", sep = "/"))

#Histogram of PValue
peridotPlotHist(res)

#MA Plot
peridotPlotMA(res, FileConfig)

#Volcano Plot
peridotPlotVolcano(res, FileConfig)

dev.off()

#Save RData
save(res, file = "DESeq.RData")
