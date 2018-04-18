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
  plot( p1, col="skyblue", main = "R-Peridot: Histogram", xlab = "Values")  # first histogram
  plot( p2, col=scales::alpha('red',.5), add=T)
  legend('topleft', c("PValue", "FDR(padjust)"), fill = c("skyblue", scales::alpha('red',.5)), bty = 'o', border = NA, cex = 0.8, bg = "white")
}

#Colors and Texts for R-Peridot Plots
xlim = ylim = c(-1, 1) * quantile(abs(resFinal$log2FoldChange[is.finite(resFinal$log2FoldChange)]), 
                           probs = 0.99) * 1.1

colFP = ifelse(resFinal$padj >= FileConfig$fdr & resFinal$pvalue >= FileConfig$pValue, "gray32", "red3")

colF = ifelse(resFinal$padj >= FileConfig$fdr, "gray32", "red3")

colP = ifelse(resFinal$pvalue >= FileConfig$pValue, "gray32", "red3")

col = "gray32"

textLegFP = c(paste("FDR(padj) < ", FileConfig$fdr, " & pvalue < ", FileConfig$pValue, sep = ""))

textLegP = c(paste("pvalue < ", FileConfig$pValue, sep = ""))

textLegF = c(paste("FDR(padj) < ", FileConfig$fdr, sep = ""))

textLine = paste("log2FoldChange = mod(", round(log2(FileConfig$foldChange), 3), ")", sep = "")

pardefault = par(no.readonly = T)

#Set Legend outside of plot
peridotPar <- function(){
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
}

#MA Plot Function
peridotPlotMA <- function(res, config){
  main = "R-Peridot: MA Plot"
  
  par(oma = c(3,1,1,1))
  
  if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colFP))
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegFP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colP))
    
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottom', c(textLegP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colF))
    
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottom', c(textLegF, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colFP))
    
    peridotPar()
    
    legend('bottom', textLegFP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = col))
    
    abline(h=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottom', textLine, col = "blue", bty = 'n', lty = 1, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colP))
    
    peridotPar()
    
    legend('bottom', textLegP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$fdr != 0){
    with(res, plot(log(baseMean), pmax(ylim[1], pmin(ylim[2], res$log2FoldChange)), 
                   pch = ifelse(res$log2FoldChange < ylim[1], 6, ifelse(res$log2FoldChange > ylim[2], 2, 20)), 
                   main = main, ylab = expression(log[2] ~ Fold ~ Change), ylim = ylim, col = colF))
    peridotPar()
    
    legend('bottom', textLegF, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }
  par(pardefault)
}

#Volcano Plot Function
peridotPlotVolcano <- function(res, config){
  main = "R-Peridot: Volcano Plot"
  
  par(oma = c(3,1,1,1))
  
  if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colFP))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegFP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colP))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colF))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegF, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colFP))
    peridotPar()
    
    legend('bottomright', textLegFP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = col))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', textLine, col = "blue", bty = 'n', lty = 1, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colP))
    peridotPar()
    
    legend('bottomright', textLegP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pvalue), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colF))
    peridotPar()
    
    legend('bottomright', textLegF, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }
  par(pardefault)
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
if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & pvalue <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & pvalue <FileConfig$pValue))
}else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & padj < FileConfig$fdr))
}else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resFinal, c(pvalue <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$foldChange != 0){
  resSig = subset(resFinal, abs(log2FoldChange) > log2(FileConfig$foldChange))
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
