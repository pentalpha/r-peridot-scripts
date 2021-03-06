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
FileConfigPath = paste(localDir, "config.txt", sep = "/")
FileConfig = read.table(FileConfigPath, header = TRUE, row.names = 1, sep = "|")
minPackages = FileConfig$minimumPackagesForConsensus

setDESeq = NULL
setEBSeq = NULL
setEdgeR = NULL
setDESeq2 = NULL
setsSeq = NULL

if(file.exists(paste(inputFilesDir, "/DESeq.AnalysisModule/1-res.tsv", sep = ""))){
  setDESeq = read.csv(paste(inputFilesDir, "/DESeq.AnalysisModule/1-res.tsv", sep = ""), sep = "\t")
}

if(file.exists(paste(inputFilesDir, "/EBSeq.AnalysisModule/1-res.tsv", sep = ""))){
  setEBSeq = read.csv(paste(inputFilesDir, "/EBSeq.AnalysisModule/1-res.tsv", sep = ""), sep = "\t")
}

if(file.exists(paste(inputFilesDir, "/edgeR.AnalysisModule/1-res.tsv", sep = ""))){
  setEdgeR = read.csv(paste(inputFilesDir, "/edgeR.AnalysisModule/1-res.tsv", sep = ""), sep = "\t")
}

if(file.exists(paste(inputFilesDir, "/DESeq2.AnalysisModule/1-res.tsv", sep = ""))){
  setDESeq2 = read.csv(paste(inputFilesDir, "/DESeq2.AnalysisModule/1-res.tsv", sep = ""), sep = "\t")
}

if(file.exists(paste(inputFilesDir, "/sSeq.AnalysisModule/1-res.tsv", sep = ""))){
  setsSeq = read.csv(paste(inputFilesDir, "/sSeq.AnalysisModule/1-res.tsv", sep = ""), sep = "\t")
}

library(limma)

set1<-row.names(setDESeq)
set2<-row.names(setEBSeq)
set3<-row.names(setEdgeR)
set4<-row.names(setDESeq2)
set5<-row.names(setsSeq)

ncol = 0
sets <- as.character()
names <- as.character()
alias <- as.character()

if(!is.null(set1)){
  sets <- c(sets, set1)
  ncol <- ncol + 1
  colset1 = ncol
  names <- c(names, "DESeq")
  alias <- c(alias, "DESeq")
}
if(!is.null(set2)){
  sets <- c(sets, set2)
  ncol <- ncol + 1
  colset2 = ncol
  names <- c(names, "EBSeq")
  alias <- c(alias, "EBSeq")
}
if(!is.null(set3)){
  sets <- c(sets, set3)
  ncol <- ncol + 1
  colset3 = ncol
  names <- c(names, "edgeR")
  alias <- c(alias, "edgeR")
}
if(!is.null(set4)){
  sets <- c(sets, set4)
  ncol <- ncol + 1
  colset4 = ncol
  names <- c(names, "DESeq2")
  alias <- c(alias, "DESeq2")
}
if(!is.null(set5)){
  sets <- c(sets, set5)
  ncol <- ncol + 1
  colset5 = ncol
  names <- c(names, "sSeq")
  alias <- c(alias, "sSeq")
}

# What are the possible letters in the universe?
universe <- sort(unique(sets))

# Generate a matrix, with the sets in columns and possible letters on rows
Counts <- matrix(0, nrow=length(universe), ncol=ncol)

# Populate the said matrix
for (i in 1:length(universe)) {
  if(!is.null(set1)){
    Counts[i,colset1] <- universe[i] %in% set1
  }
  if(!is.null(set2)){
    Counts[i,colset2] <- universe[i] %in% set2
  }
  if(!is.null(set3)){
    Counts[i,colset3] <- universe[i] %in% set3
  }
  if(!is.null(set4)){
    Counts[i,colset4] <- universe[i] %in% set4
  }
  if(!is.null(set5)){
    Counts[i,colset5] <- universe[i] %in% set5
  }
}

# Name the columns with the sample names
#colnames(Counts) <- c("set1","set2","set3")
colnames(Counts) <- names

interPacks = apply(Counts, MARGIN = 1, FUN = function(x){
  paste(alias[x == 1], collapse = ",")
})

#Remover o 0 de valores de fora do universo amostral
vCounts = vennCounts(Counts)
vCounts[1, length(colnames(vCounts))] = NA

# Specify the colors for the sets
cols<-c("Red", "Green", "Blue", "Black", "Pink")

png(filename = paste(outputFilesDir, "/2-vennDiagramPlot.png", sep = ""), width=600, height=600)

vennDiagram(vCounts, circle.col=cols)

dev.off()

pdf(file = paste(outputFilesDir, "/3-plots.pdf", sep = ""))

vennDiagram(vCounts, circle.col=cols)

dev.off()

##### END PLOT #####

##### INTERSCTS #####
listSets <- list()

if(!is.null(set1)){
  listSets <- append(listSets, list(set1))
}
if(!is.null(set2)){
  listSets <- append(listSets, list(set2))
}
if(!is.null(set3)){
  listSets <- append(listSets, list(set3))
}
if(!is.null(set4)){
  listSets <- append(listSets, list(set4))
}
if(!is.null(set5)){
  listSets <- append(listSets, list(set5))
}

intersect2 <- function(...) {
  args <- list(...)
  nargs <- length(args)
  if(nargs <= 1) {
    if(nargs == 1 && is.list(args[[1]])) {
      do.call("intersect2", args[[1]])
    } else {
      stop("cannot evaluate intersection fewer than 2 arguments")
    }
  } else if(nargs == 2) {
    intersect(args[[1]], args[[2]])
  } else {
    intersect(args[[1]], intersect2(args[-1]))
  }
}

if(length(listSets) > 1){
  interSets = intersect2(listSets)
}else{
  interSets = listSets
}

#write.table(interSets, paste(outputFilesDir, "1-Intersect.tsv", sep = "/"), sep = "\t", row.names = F, col.names = F)

interUniverse = data.frame(list = universe, packages = interPacks)

ord = apply(interUniverse, 1, function(x){
  nchar(x[2])
})

interUniverse = interUniverse[order(ord, decreasing = T),]

colnames(interUniverse) <- c('id', 'packs')
colnames(interUniverse)

write.table(interUniverse, paste(outputFilesDir, "1-Intersect.tsv", sep = "/"), sep = "\t", row.names = F, col.names = T)


enough_packages = apply(interUniverse, 1, function(x){
  length(unlist(strsplit(x[2],","))) >= minPackages
})
consensus <- interUniverse[enough_packages, ]
if (length(consensus[,1]) > 0){
  write.table(consensus, paste(outputFilesDir, "4-Consensus.tsv", sep = "/"), sep = "\t", row.names = F, col.names = T)
}else{
  "No consensus to write."
}
