args = commandArgs(trailingOnly = F)

localDir <- args[length(args)-3]

localDir

inputFilesDir <- args[length(args)-2]

inputFilesDir

outputFilesDir <- args[length(args)-1]

outputFilesDir

notFirstRun <- args[length(args)]

notFirstRun

options(bitmapType='cairo')
setwd(localDir)

#Get directory
paramFile = paste(localDir, "config.txt", sep = "/")
params = read.table(paramFile, header = TRUE, row.names = 1, sep = "|")
params

genelistinput = paste(inputFilesDir, "VennDiagram.PostAnalysisModule/1-Intersect.tsv", sep = "/");

genelist = inter = read.table(file = genelistinput, header = F, sep = "\t")

library(clusterProfiler)
library(ggplot2)

refOrganism = params$referenceOrganism
stopifnot(refOrganism == "Human" || refOrganism == "Mouse" || refOrganism == "Fly")

orgDBName<-NULL
universe<-NULL

if(refOrganism == "Human"){
  require(org.Hs.eg.db)
  universe = org.Hs.egENSEMBL
  orgDBName = "org.Hs.eg.db"
  keggOrg = "hsa"
}else if(refOrganism == "Mouse"){
  require(org.Mm.eg.db)
  universe = org.Mm.egENSEMBL
  orgDBName = "org.Mm.eg.db"
  keggOrg = "mmu"
}else if(refOrganism == "Fly"){
  require(org.Dm.eg.db)
  universe = org.Dm.egENSEMBL
  orgDBName = "org.Dm.eg.db"
  keggOrg = "dme"
}

eg = bitr(gsub("\\..*","",genelist[,1]), fromType= as.character(params$idType), toType = c("SYMBOL","ENTREZID", "UNIPROT"), OrgDb = orgDBName, drop = T)

head(eg$ENTREZID)

mappedGenes = mappedkeys(universe)

ego = enrichGO(gene = eg$ENTREZID, universe = mappedGenes, OrgDb = orgDBName, ont = "MF", readable = T, pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff = 1)
ego@result = subset(ego@result, (pvalue < params$pValue & qvalue < params$fdr))
head(ego@result)

ego2 = enrichGO(gene = eg$ENTREZID, universe = mappedGenes, OrgDb = orgDBName, ont = "CC", readable = T, pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff = 1)
ego2@result = subset(ego2@result, (pvalue < params$pValue & qvalue < params$fdr))
head(ego2@result)

ego3 = enrichGO(gene = eg$ENTREZID, universe = mappedGenes, OrgDb = orgDBName, ont = "BP", readable = T, pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff = 1)
ego3@result = subset(ego3@result, (pvalue < params$pValue & qvalue < params$fdr))
head(ego3@result)

egos = c(ego,ego2,ego3)
lenCategoryEgos = c()
for (enGo in egos){
  lenCategoryEgos = c(lenCategoryEgos, length(enGo@result$ID))
}
lensDescriptionChar = c()
i = 1
for (enGo in egos){
  applyed = sapply(enGo@result$Description, nchar)
  newLen = 0
  if (length(applyed) > 0){
    newLen = max(applyed)
  }
  lensDescriptionChar = c(lensDescriptionChar, newLen)
  i = i + 1
}

height = max(lenCategoryEgos, 48, na.rm = T)
width = max(lensDescriptionChar, 48, na.rm = T)

if(height/6 <= 8){
  height = 10
}else if(height/6 > 50){
  height = 50
}else{
  height = height/6
}

if(width/6 <= 8){
  width = 10
}else if(width/6 > 50){
  width = 50
}else{
  width = width/6
}

titles = c("Ontology = MF", "Ontology = CC", "Ontology = BP")
pngs = c("enrichGOMF.png", "enrichGOCC.png", "enrichGOBP.jpg")

i = 1
for (enGo in egos){
  if(lenCategoryEgos[i] > 0){
    dotplot(enGo, title = titles[i], showCategory = lenCategoryEgos[i], colorBy = "pvalue")

    ggsave(filename = paste(outputFilesDir, pngs[i], sep = "/"), width = width, height = height)
  }
  i = i + 1
}

#0,077898551 inches per char
pdf(file = paste(outputFilesDir, "enrich.pdf", sep = "/"), width = width, height = height)
i = 1
for (enGo in egos){
  if(lenCategoryEgos[i] > 0){
    print(dotplot(enGo, title = titles[i], showCategory = lenCategoryEgos[i], colorBy = "pvalue"))
  }
  i = i + 1
}

dev.off()

##### KEGG #####
eg2np <- bitr_kegg(eg$ENTREZID, fromType='kegg', toType='ncbi-geneid', organism=keggOrg)

kegg <- enrichKEGG(eg2np$kegg, organism = keggOrg, pvalueCutoff = 1, qvalueCutoff = 1)

kegg@result = subset(kegg@result, (pvalue < params$pValue & qvalue < params$fdr))

if(length(kegg@result$ID) > 0){
  dotplot(kegg, title = "KEGG Ontology", showCategory = length(kegg@result$ID), colorBy = "pvalue")
  
  ggsave(filename = paste(outputFilesDir, "KEGG_GO.jpg", sep = "/"), width = width, height = height)
}else{
  print("No enriched terms found for KEGG Ontology")
}

urlPath = vector()

urlPath = sapply(kegg@result$ID, function(x){
  paste0("http://www.kegg.jp/kegg-bin/show_pathway?", x, '/', kegg[x, "geneID"])
})

write.table(x = data.frame(urlPath), file = paste(outputFilesDir, "urlPaths.tsv", sep = "/"), sep = "\t", quote = F)

tables <- c("enrichGOMF.tsv", "enrichGOCC.tsv", "enrichGOBP.tsv")
i = 1
for (enGo in egos){
  if(lenCategoryEgos[i] > 0){
    write.table(x = enGo@result, file = paste(outputFilesDir, tables[i], sep = "/"), sep = "\t", quote = F)
  }else{
    print(paste("No enriched terms found for ", titles[i]))
  }
  i = i + 1
}