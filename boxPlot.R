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

sampleNames = rownames(pasillaCountTable)

head(sampleNames)

pasillaDesign = read.table(conditions, header=TRUE, row.names=1)
pasillaDesign

#Normalization
pasillaCountTable = as.data.frame(lapply(pasillaCountTable, function(x) (x/sum(x))*10000000))

names = colnames(pasillaCountTable)

row.names(pasillaCountTable) = sampleNames

head(pasillaCountTable)

lev = levels(pasillaDesign$condition)

color.code<-colorRampPalette(c('blue','yellow'), space="rgb")(length(lev))

conditionColor<-function(condition, conditions, cols){
  finalColor <- cols[1]
  for (i in 1:length(conditions)){
    if (condition == conditions[i]){
      finalColor <- cols[i]
    }
  }
  finalColor
}
colors <- color.code[pasillaDesign$condition]
colors
levelColors = c()
levelIndex = 1
for(cond in lev){
  color <- conditionColor(cond, pasillaDesign$condition, colors)
  levelColors <- c(levelColors, color)
  print(color)
}

levelColors

png(filename = saveTo, width = 600, height = 600)
par(xpd = T, mar = par()$mar + c(0.5,0,-3,7))
boxplot(pasillaCountTable, outline = F, col = colors, axes = FALSE, axisnames = FALSE)

getPos<-function(pos, bounds, coord){
  xmin <- bounds[1]
  xmax <- bounds[2]
  ymin <- bounds[3]
  ymax <- bounds[4]
  xsize <- xmax - xmin
  ysize <- ymax - ymin

  position <- xmin + (xsize * pos)
  if (coord == "y"){
    position <- ymin + (ysize * pos)
  }
  position
}

text(seq(1, length(pasillaDesign$condition), by=1), getPos(0.027,par("usr"),"y"), 
  labels = row.names(pasillaDesign), srt = 45, adj = c(1.1,1.1), xpd = TRUE, cex=.9)
text(getPos(0.5,par("usr"),"x"),getPos(1.01,par("usr"),"y"), labels = "(rows with only zeros removed)", adj=c(0.5,0.5))
axis(2)

ncols <- ceiling(length(levelColors) / 3)
ncols
par("usr")
x <- getPos(1.013,par("usr"),"x")
y <- getPos(1.0,par("usr"),"y")
legend(x, y, 
  title="Conditions", lev, fill=levelColors, ncol=1,cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()