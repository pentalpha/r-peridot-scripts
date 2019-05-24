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
#Get directory
FileConfig = getwd()

intersectFilePath = paste(inputFilesDir, "VennDiagram.PostAnalysisModule/1-Intersect.tsv", sep = "/")
interUniverse = read.table(intersectFilePath, header=TRUE, row.names=1)

peridotCountTable = read.table(paste(inputFilesDir, "rna-seq-input.tsv", sep = "/"), header=TRUE, row.names=1)

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

# Normalizar o dado de entrada #
peridotCountTable = as.data.frame(lapply(peridotCountTable, function(x) (x/sum(x))*10000000))

rownames(peridotCountTable) = geneNames

# Plot normalized counts for each differentially expressed gene
dir.create(file.path(outputFilesDir, 'countPlots'))
apply(X = as.data.frame(interUniverse[,1]), MARGIN = 1, FUN = function(x){
  png(filename = paste(outputFilesDir, "/countPlots/", x, ".png", sep = ""), width = 600, height = 600)
  
  plot(x = as.integer(peridotConditions$condition) + runif(ncol(peridotCountTable),-.05,.05), y = peridotCountTable[rownames(peridotCountTable) == x,], main = x, xlab = "group", ylab = "normalized count", xlim=c(.5,max(as.integer(peridotConditions$condition))+.5), log="y", xaxt = "n")
  axis(1, at=seq_along(levels(as.factor(peridotConditions$condition))), levels(as.factor(peridotConditions$condition)))
  
  dev.off()  
})
