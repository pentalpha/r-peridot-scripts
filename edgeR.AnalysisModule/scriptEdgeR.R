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

  #Load edgeR
  library(edgeR)

  edesign <- model.matrix(~peridotConditions$condition)
  peridotCountTable = as.matrix(peridotCountTable)

  #Create DGEList
  e <- DGEList(counts=peridotCountTable)

  #Calculate nomalization factors
  e <- calcNormFactors(e)

  normCounts <- e$counts/e$samples$norm.factors

  BaseMeanVect = rowMeans(normCounts)

  #Estimates a common negative binomial dispersion parameter for a DGE dataset
  e <- estimateGLMCommonDisp(e, edesign)

  #Estimates the abundace-disersion trend
  e <- estimateGLMTrendedDisp(e, edesign)

  #Compute an empirical Bayes estimate of the negative binomial dispersion parameter
  e <- estimateGLMTagwiseDisp(e, edesign)

  ## Fit the model, testing the coefficient for the treated vs untreated comparison
  efit <- glmFit(e, edesign)

  efit <- glmLRT(efit)#, coef="conditiontreated")

  ## Make a table of results
  etable <- topTags(efit, n=nrow(e))$table

  ## Create a column FoldChange
  etable$FoldChange = 2^etable$logFC

  ## List of conditions A
  factA = peridotConditions$condition==levels(peridotConditions$condition)[2]

  ## Columns with condition A
  colA = normCounts[,factA]

  ## List of conditions B
  factB = peridotConditions$condition==levels(peridotConditions$condition)[1]

  ## Columns with condition B
  colB = normCounts[,factB]

  ## Means of condition A
  baseMeanA = colA
  if(!(is.null(dim(colA)))){
    baseMeanA = rowMeans(colA)
  }

  ## Means of condition B
  baseMeanB = colB
  if(!(is.null(dim(colB)))){
    baseMeanB = rowMeans(colB)
  }

  ## Create data frame of baseMeans
  BaseMeandf = as.data.frame(BaseMeanVect)

  colnames(BaseMeandf) = "baseMean"

  BaseMeandf$baseMeanA = baseMeanA

  BaseMeandf$baseMeanB = baseMeanB

  df = merge(BaseMeandf, etable, by = 'row.names')

  res = df[c("Row.names", "baseMean", "baseMeanA", "baseMeanB", "FoldChange", "logFC", "logCPM", "LR", "PValue", "FDR")]

  rownames(res) = res$Row.names

  res$Row.names = NULL

  res$logCPM = NULL

  res$LR = NULL
}else{
  load(file = "edgeR.RData")
}

#Histogram PValue and FDR Function
peridotPlotHist <- function(res){
  p1 <- with(res, hist(PValue, breaks=100, plot = F))
  p2 <- with(res, hist(FDR, breaks=100, plot = F))
  plot( p1, col="skyblue", main="Histogram", xlab = "Values")  # first histogram
  plot( p2, col=scales::alpha('red',.5), add=T)
  legend('topleft', c("PValue", "FDR"), fill = c("skyblue", scales::alpha('red',.5)), bty = 'o', border = NA, cex = 0.8, bg = "white")
}

#MA Plot Funcition
peridotPlotMA <- function(res, config){
  with(res, plot(log(baseMean), logFC, pch=20, main="MA Plot"))
  
  if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr & PValue<FileConfig$pValue), points(log(baseMean), logFC, pch=20, col="red"))
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
    with(subset(res, PValue<FileConfig$pValue), points(log(baseMean), logFC, pch=20, col="red"))
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr), points(log(baseMean), logFC, pch=20, col="red"))
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('bottomright', c(paste("FDR < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr & PValue<FileConfig$pValue), points(log(baseMean), logFC, pch=20, col="red"))
    legend('bottomright', paste("FDR < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0){
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('bottomright', paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(res, PValue<FileConfig$pValue), points(log(baseMean), logFC, pch=20, col="red"))
    legend('bottomright', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr), points(log(baseMean), logFC, pch=20, col="red"))
    legend('bottomright', paste("FDR < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }
}

#Volcano Plot Function
peridotPlotVolcano <- function(res, config){
  with(res, plot(logFC, -log10(PValue), pch=20, main="Volcano plot"))
  
  if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr & PValue<FileConfig$pValue), points(logFC, -log10(PValue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
    with(subset(res, PValue<FileConfig$pValue), points(logFC, -log10(PValue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('topleft', c(paste("pvalue < ", FileConfig$pValue, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr), points(logFC, -log10(PValue), pch=20, col="red"))
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('topleft', c(paste("FDR < ", FileConfig$fdr, sep = ""), paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")), col = c("red", "blue"), bty = 'o', pch = c(16, NA), lty = c(NA, 1), bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr & PValue<FileConfig$pValue), points(logFC, -log10(PValue), pch=20, col="red"))
    legend('topleft', paste("FDR < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$foldChange != 0){
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    legend('topleft', paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = ""), col = "blue", bty = 'o', lty = 1, bg = "white", cex = 0.8)
  }else if(FileConfig$pValue != 0){
    with(subset(res, PValue<FileConfig$pValue), points(logFC, -log10(PValue), pch=20, col="red"))
    legend('topleft', paste("pvalue < ", FileConfig$pValue, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }else if(FileConfig$fdr != 0){
    with(subset(res, FDR<FileConfig$fdr), points(logFC, -log10(PValue), pch=20, col="red"))
    legend('topleft', paste("FDR < ", FileConfig$fdr, sep = ""), col = "red", bty = 'o', pch = 16, bg = "white", cex = 0.8)
  }
}

#Call Histogram Plot
png(filename = paste(outputFilesDir, "2-histogram.png", sep = "/"), width=600, height=600)

peridotPlotHist(res)

dev.off()

#Call MA Plot
png(filename = paste(outputFilesDir, "3-MAPlot.png", sep = "/"), width=600, height=600)

peridotPlotMA(res, FileConfig)

dev.off()

#Call Volcano Plot
png(filename = paste(outputFilesDir, "4-volcanoPlot.png", sep = "/"), width=600, height=600)

peridotPlotVolcano(res, FileConfig)

dev.off()

#Subset com PValue < FileConfig$pValue, Fold Change < FileConfig$foldChange e FDR < FileConfig$fdr
if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(res, abs(logFC) > log2(FileConfig$foldChange) & PValue < FileConfig$pValue & FDR < FileConfig$fdr)
}else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
  resSig = subset(res, abs(logFC) > log2(FileConfig$foldChange) & PValue <FileConfig$pValue)
}else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
  resSig = subset(res, abs(logFC) > log2(FileConfig$foldChange) & FDR < FileConfig$fdr)
}else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(res, PValue <FileConfig$pValue & FDR < FileConfig$fdr)
}else if(FileConfig$foldChange != 0){
  resSig = subset(res, abs(logFC) > log2(FileConfig$foldChange))
}else if(FileConfig$pValue != 0){
  resSig = subset(res, PValue <FileConfig$pValue)
}else if(FileConfig$fdr != 0){
  resSig = subset(res, FDR < FileConfig$fdr)
}else{
  resSig = res
}

##Remove NA
resSig = na.omit(resSig)
#

head(resSig)

if(FileConfig$tops > 0){
  topRes = head(resSig, n = FileConfig$tops)

  write.table(topFDR, paste(outputFilesDir, "TopResults.tsv", sep = "/"), sep = "\t")
}

##Create files csv
#if(length(resSig$FDR > 0)){
  write.table(resSig, paste(outputFilesDir, "1-res.tsv", sep = "/"), sep = "\t")
#}
#

pdf(file = paste(outputFilesDir, "5-plots.pdf", sep = "/"))

#Histogram PValue and FDR
peridotPlotHist(res)

#MA Plot
peridotPlotMA(res, FileConfig)

#Volcano Plot
peridotPlotVolcano(res, FileConfig)

dev.off()

save(etable, res, file = "edgeR.RData")
