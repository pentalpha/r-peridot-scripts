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
#Get directory
FileConfig = getwd()

#Read file example
peridotCountTable = read.table(paste(inputFilesDir, "rna-seq-input.tsv", sep = "/"), header=TRUE, row.names=1 )

geneNames = rownames(peridotCountTable)

peridotConditions = read.table(paste(inputFilesDir, "condition-input.tsv", sep = "/"), header=TRUE, row.names=1)
peridotConditions

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

rld = log2(peridotCountTable+0.99)

cor = cor(as.matrix(rld))

# Normalizar o dado de entrada #
peridotCountTable = as.data.frame(lapply(peridotCountTable, function(x) (x/sum(x))*1000000))

rownames(peridotCountTable) = geneNames

#############################

# Gerar boxplot #

lev = levels(peridotConditions$condition)

color.code<-colorRampPalette(c('blue','red'), space="rgb")(length(lev))

peridotCountTableNA = as.data.frame(apply(peridotCountTable, c(1, 2), function(x){
  #print(x)
  if(x == 0){
    x = NA
  }else{
    x = x
  }
}))

png(filename = paste(outputFilesDir, "G-BoxPlot.png", sep = "/"), width=600, height=600)

boxplot(log2(peridotCountTableNA), outline = F, col=color.code[peridotConditions$condition], main = "Boxplot", las=2)

dev.off()

###############################

library(pvclust)

library(gplots)

library(RColorBrewer)

# Abrir o arquivo de miRNAs achados nos pacotes do R-peridot #
intersectFile = paste(inputFilesDir, "VennDiagram.PostAnalysisModule/1-Intersect.tsv", sep = "/")

inter = read.table(file = intersectFile, header = F, sep = "\t")

################################

# Calcular o PCA #
tperidot = t(peridotCountTable)

pca = prcomp(tperidot)

png(filename = paste(outputFilesDir, "E-PCA.png", sep = "/"), width=600, height=600)

plot(pca$x[,1], pca$x[,2], xlab="PCA 1", ylab="PCA 2",type="p", pch=19, col=color.code[peridotConditions$condition] , cex=1.0, xlim=c(min(pca$x[,1])*1.1, max(pca$x[,1])*1.1), ylim=c(min(pca$x[,2]*1.1), max(pca$x[,2])*1.1), main = "PCA")

text(pca$x[,1] -4,pca$x[,2]-1, rownames(tperidot),cex=0.7, pos = 1)

dev.off()

################################

# Salvar arquivos normalizados #

rownames(peridotCountTable) = geneNames

write.table(peridotCountTable, paste(outputFilesDir, "F-NormalizedCounts.tsv", sep = "/"), sep = "\t", row.names = T)

###############################

length(inter[,1])

if(length(inter[,1]) > 6){
  # Calcular dendrograma e heatmap só funciona com mais de 6 samples encontrados #

  #peridotCountTable[] = lapply(peridotCountTable, function(x) as.integer(x))

  subInter = intersect(rownames(peridotCountTable), inter[,1])

  d = as.matrix(peridotCountTable[subInter,])

  d.pv = pvclust(d, nboot = 1000, parallel = TRUE, method.hclust = "complete", method.dist = "euclidean")

  png(filename = paste(outputFilesDir, "D-Dendrogram.png", sep = "/"), width=600, height=600)

  plot(d.pv)

  dev.off()
  
  z = t(scale(t(d)))
  
  hclustfunc = function(x) hclust(x, method = "complete")
  
  distfunc = function(x) dist(x, method = "maximum")
  
  cols = c("dodgerblue3", "firebrick3")[peridotConditions$condition]
  
  d2 = as.matrix(apply(d, c(1,2), function(x){
    if(x == 0){
      x = 0.01
    }else{
      x = x
    }
  }))
  
  peridotPar <- function(){
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  }
  
  
  png(filename = paste(outputFilesDir, "A-HeatMapScale.png", sep = "/"), width=600, height=600)

  heatmap.2(z, hclustfun = hclustfunc, distfun = distfunc, dendrogram = "both", key = T, keysize = 1.4,
            key.par=list(cex=0.5), col = greenred(200), scale = "row", trace = "none", cexCol = 0.7, srtCol = 90,
            density.info = 'histogram', main = "R-Peridot: HeatMap", labRow = "", margins = c(10,7),
            ColSideColors = cols)
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )
  
  dev.off()
  
  png(filename = paste(outputFilesDir, "B-HeatMapCor.png", sep = "/"), width=600, height=600)
  
  heatmap.2(cor, symm = T, col = colorRampPalette(c("darkblue", "white"))(100),
            labCol = colnames(cor), labRow = colnames(cor),
            distfun = function(c) as.dist(1 - c), trace = "none", Colv = T,
            cexRow = 0.9, cexCol = 0.9, key = F, font = 2, RowSideColors = cols,
            ColSideColors = cols, main = "R-Peridot: Samples correlation Heatmap", margins = c(10,8))
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )
  
  dev.off()
  
  png(filename = paste(outputFilesDir, "C-HeatMapLog2.png", sep = "/"), width=600, height=600)
  
  sig.dat = rld[subInter,]
  
  annC = data.frame(condition=peridotConditions)
  
  rownames(annC) = colnames(sig.dat)
  
  heatmap.2(as.matrix(sig.dat), scale = "row", trace = "none", margins = c(10, 8), 
            col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(256), 
            cexRow = 0.9, cexCol = 0.9, key = T, keysize = 1.4, key.par=list(cex=0.5), 
            offsetRow = T, offsetCol = T,reorderfun = function(d,w) reorder(d,w, agglo.FUN = mean), 
            main = expression(R-Peridot: ~ log[2] ~ (Count ~ reads + 0.99) ~ Heatmap), lhei = c(2,7), lwid = c(2,7), ColSideColors = cols)
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )

  dev.off()

  PDFheight = nrow(d)/8
  PDFwidth = ncol(d)/8
  if(PDFheight <= 8) {PDFheight = 10}
  if(PDFwidth <= 8) {PDFwidth = 10}

  pdf(file = paste(outputFilesDir, "3-boxplot-pca-dendrogram.pdf", sep = "/"))

  par(cex.axis=0.8)

  boxplot(log2(peridotCountTableNA), outline = F, col=color.code[peridotConditions$condition], main = "Boxplot", las=2)

  plot(pca$x[,1], pca$x[,2], xlab="PCA 1", ylab="PCA 2",type="p", pch=19, col=color.code[peridotConditions$condition] , cex=1.0, xlim=c(min(pca$x[,1])*1.1, max(pca$x[,1])*1.1), ylim=c(min(pca$x[,2]*1.1), max(pca$x[,2])*1.1), main = "PCA")

  text(pca$x[,1] -4,pca$x[,2]-1, rownames(tperidot),cex=0.7, pos = 1)

  plot(d.pv, cex = 0.8, print.pv = F, main = "Dendrogram")

  dev.off()

  pdf(file = paste(outputFilesDir, "1-HeatMapCor.pdf", sep = "/"), height = PDFheight, width = PDFwidth)
  
  #### HeatMap 1 ####
  
  heatmap.2(cor, symm = T, col = colorRampPalette(c("darkblue", "white"))(100),
            labCol = colnames(cor), labRow = colnames(cor),
            distfun = function(c) as.dist(1 - c), trace = "none", Colv = T,
            cexRow = 0.9, cexCol = 0.9, key = F, font = 2, RowSideColors = cols,
            ColSideColors = cols, main = "R-Peridot: Samples correlation Heatmap", margins = c(10,8))
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )
  
  dev.off()

  pdf(file = paste(outputFilesDir, "2-HeatMaps.pdf", sep = "/"), height = PDFheight, width = PDFwidth)

  #### HeatMap 2 ####

  heatmap.2(as.matrix(sig.dat), scale = "row", trace = "none", margins = c(10, 8), 
            col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(256), 
            cexRow = 0.9, cexCol = 0.9, key = T, keysize = 1.4, key.par=list(cex=0.5), 
            offsetRow = T, offsetCol = T,reorderfun = function(d,w) reorder(d,w, agglo.FUN = mean), 
            main = expression(R-Peridot: ~ log[2] ~ (Count ~ reads + 0.99) ~ Heatmap), lhei = c(2,PDFheight/2), lwid = c(2,PDFwidth/2), ColSideColors = cols)
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )
  
  #### HeatMap 3 ####
  heatmap.2(z, hclustfun = hclustfunc, distfun = distfunc, dendrogram = "both", key = T, keysize = 1.4, 
            key.par=list(mar=c(3,1,3,1)), col = greenred(200), scale = "row", trace = "none", cexRow = 0.75, cexCol = 0.9, srtCol = 90, 
            density.info = 'histogram', main = "R-Peridot: Scale(Count Reads) Heatmap", margins = c(10, 8),
            lhei = c(2,PDFheight/2), lwid = c(2,PDFwidth/2), ColSideColors = cols)
  
  peridotPar()
  
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(peridotConditions$condition), # category labels
         col = levels(as.factor(cols)),  # color key
         lty= 1,             # line style
         lwd = 8            # line width
  )
  
  dev.off()

}else{
  stop("Number of results less than 6")
}
