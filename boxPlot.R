args = commandArgs(trailingOnly = F)

countReads <- args[length(args)-2]

countReads

conditions <- args[length(args)-1]

conditions

saveTo <- args[length(args)]

saveTo

options(bitmapType='cairo')

#Read file example
pasillaCountTable = read.table(countReads, header=TRUE, row.names=1 )

geneNames = rownames(pasillaCountTable)

head(geneNames)

pasillaDesign = read.table(conditions, header=TRUE, row.names=1)
pasillaDesign

#Ignore samples with "not-use" indicated
#first, remove they from the conditions table
pasillaDesign <- subset(pasillaDesign, condition != "not-use")
#then, remove from the counts table
for(i in colnames(pasillaCountTable)){
  iContainsNotUse = length(grep("not.use", as.name(i))) > 0
  if(iContainsNotUse){
    #erases the column
    pasillaCountTable[, i] = NULL
  }
}
#Finally, drop unused levels (not-use levels)
pasillaDesign = droplevels(pasillaDesign)
pasillaDesign

#Normalization
pasillaCountTable = as.data.frame(lapply(pasillaCountTable, function(x) (x/sum(x))*10000000))

names = colnames(pasillaCountTable)

row.names(pasillaCountTable) = geneNames

head(pasillaCountTable)

lev = levels(pasillaDesign$condition)

color.code<-colorRampPalette(c('blue','red'), space="rgb")(length(lev))

png(filename = saveTo, width=600, height=600)

boxplot(pasillaCountTable, outline = F, col = color.code[pasillaDesign$condition])

dev.off()