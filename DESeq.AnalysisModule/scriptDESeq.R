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
  res = nbinomTest( cds, levels(peridotConditions$condition)[1], levels(peridotConditions$condition)[2] )

  head(res)

}else{
  load(file = "DESeq.RData")
}

#Histogram PValue and FDR Function
peridotPlotHist <- function(res){
  p1 <- with(res, hist(pval, breaks=100, plot = F))
  p2 <- with(res, hist(padj, breaks=100, plot = F))
  plot( p1, col="skyblue", main="R-Peridot: Histogram", xlab = "Values")  # first histogram
  plot( p2, col=scales::alpha('red',.5), add=T)
  legend('topleft', c("PValue", "FDR(padjust)"), fill = c("skyblue", scales::alpha('red',.5)), bty = 'o', border = NA, cex = 0.8, bg = "white")
}


#Colors and Texts for R-Peridot Plots
xlim = c(-1, 1) * quantile(abs(res$log2FoldChange[is.finite(res$log2FoldChange)]), 
                           probs = 0.99) * 1.1

ylim = c(-1, 1) * quantile(abs(res$log2FoldChange[is.finite(res$log2FoldChange)]), 
                           probs = 0.99) * 1.1
colFP = ifelse(res$padj >= FileConfig$fdr & res$pval >= FileConfig$pValue, "gray32", "red3")

colF = ifelse(res$padj >= FileConfig$fdr, "gray32", "red3")

colP = ifelse(res$pval >= FileConfig$pValue, "gray32", "red3")

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
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colFP))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegFP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colP))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegP, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colF))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', c(textLegF, textLine), col = c("red", "blue"), bty = 'n', pch = c(20, NA), 
           lty = c(NA, 1), bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colFP))
    peridotPar()
    
    legend('bottomright', textLegFP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$foldChange != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = col))
    
    abline(v=c(log2(FileConfig$foldChange), log2(FileConfig$foldChange)*(-1)), col="blue")
    
    peridotPar()
    
    legend('bottomright', textLine, col = "blue", bty = 'n', lty = 1, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$pValue != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colP))
    peridotPar()
    
    legend('bottomright', textLegP, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }else if(FileConfig$fdr != 0){
    with(res, plot(pmax(xlim[1], pmin(xlim[2], res$log2FoldChange)), -log10(pval), 
                   pch = ifelse(res$log2FoldChange < xlim[1], 6, ifelse(res$log2FoldChange > xlim[2], 2, 20)), 
                   main = main, xlab = expression(log[2] ~ Fold ~ Change), xlim = xlim, col = colF))
    peridotPar()
    
    legend('bottomright', textLegF, col = "red", bty = 'n', pch = 20, 
           bg = "white", cex = 0.8, xpd = T, inset = c(0,0.04))
  }
  par(pardefault)
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

resAux = res

rownames(resAux) = resAux$id
resAux$id = NULL

#Subset com PValue < FileConfig$pValue, Fold Change < FileConfig$log2FoldChange e FDR < FileConfig$fdr
if(FileConfig$foldChange != 0 & FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & pval <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$foldChange != 0 & FileConfig$pValue != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & pval <FileConfig$pValue))
}else if(FileConfig$foldChange != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(abs(log2FoldChange) > log2(FileConfig$foldChange) & padj < FileConfig$fdr))
}else if(FileConfig$pValue != 0 & FileConfig$fdr != 0){
  resSig = subset(resAux, c(pval <FileConfig$pValue & padj < FileConfig$fdr))
}else if(FileConfig$foldChange != 0){
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
write.table(resSig, paste(outputFilesDir, "/1-res.tsv", sep = ""), sep = "\t")

pdf(file = paste(outputFilesDir, "5-plots.pdf", sep = "/"))

#Histogram of PValue
peridotPlotHist(res)

#MA Plot
peridotPlotMA(res, FileConfig)

#Volcano Plot
peridotPlotVolcano(res, FileConfig)

dev.off()

#Save RData
save(res, file = "DESeq.RData")
