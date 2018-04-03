args = commandArgs(trailingOnly = F)

localDir <- args[length(args)-3]

localDir

inputFilesDir <- args[length(args)-2]

inputFilesDir

outputFilesDir <- args[length(args)-1]

outputFilesDir = getwd()

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
head(ego)

ego2 = enrichGO(gene = eg$ENTREZID, universe = mappedGenes, OrgDb = orgDBName, ont = "CC", readable = T, pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff = 1)
ego2@result = subset(ego2@result, (pvalue < params$pValue & qvalue < params$fdr))
head(ego2)

ego3 = enrichGO(gene = eg$ENTREZID, universe = mappedGenes, OrgDb = orgDBName, ont = "BP", readable = T, pAdjustMethod = "BH", pvalueCutoff  = 1, qvalueCutoff = 1)
ego3@result = subset(ego3@result, (pvalue < params$pValue & qvalue < params$fdr))
head(ego3)

#Length of the number of Category
lenCategoryEgo = length(ego@result$ID)
lenCategoryEgo2 = length(ego2@result$ID)
lenCategoryEgo3 = length(ego3@result$ID)

#Length of the number of chars of description
lenDescriptionCharEgo = nchar(max(ego@result$Description))
lenDescriptionCharEgo2 = nchar(max(ego2@result$Description))
lenDescriptionCharEgo3 = nchar(max(ego3@result$Description))

height = max(lenCategoryEgo, lenCategoryEgo2, lenCategoryEgo3, 48, na.rm = T)
width = max(lenDescriptionCharEgo, lenDescriptionCharEgo2, lenDescriptionCharEgo3, 48, na.rm = T)

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

if(length(ego@result$ID) > 0){
  dotplot(ego, title = "Ontology = MF", showCategory = lenCategoryEgo, colorBy = "pvalue")

  ggsave(filename = paste(outputFilesDir, "enrichGOMF.png", sep = "/"), width = width, height = height)
}

if(length(ego2@result$ID) > 0){
  dotplot(ego2, title = "Ontology = CC", showCategory = 100, colorBy = "pvalue")
  
  ggsave(filename = paste(outputFilesDir, "enrichGOCC.png", sep = "/"), width = width, height = height)
}

if(length(ego3@result$ID) > 0){
  dotplot(ego3, title = "Ontology = BP", showCategory = length(ego3@result$ID), colorBy = "pvalue")

  ggsave(filename = paste(outputFilesDir, "enrichGOBP.jpg", sep = "/"), width = width, height = height)
}

#0,077898551 inches per char
pdf(file = paste(outputFilesDir, "enrich.pdf", sep = "/"), width = width, height = height)
if(length(ego@result$ID) > 0){
  dotplot(ego, title = "Ontology = MF", showCategory = length(ego@result$ID), colorBy = "pvalue")
}

if(length(ego2@result$ID) > 0){
  dotplot(ego2, title = "Ontology = CC", showCategory = length(ego2@result$ID), colorBy = "pvalue")
}

if(length(ego3@result$ID) > 0){
  dotplot(ego3, title = "Ontology = BP", showCategory = length(ego3@result$ID), colorBy = "pvalue")
}

dev.off()

##### KEGG #####
eg2np <- bitr_kegg(eg$ENTREZID, fromType='kegg', toType='ncbi-geneid', organism=keggOrg)

kk <- enrichKEGG(eg2np$kegg, organism = keggOrg, pvalueCutoff = 1, qvalueCutoff = 1)

kk@result = subset(kk@result, (pvalue < params$pValue & qvalue < params$fdr))

if(length(kk@result$ID) > 0){
  dotplot(kk, title = "KEGG Ontology", showCategory = length(kk@result$ID), colorBy = "pvalue")
  
  ggsave(filename = paste(outputFilesDir, "KEGG_GO.jpg", sep = "/"), width = width, height = height)
}

urlPath = vector()

urlPath = sapply(kk@result$ID, function(x){
  paste0("http://www.kegg.jp/kegg-bin/show_pathway?", x, '/', kk[x, "geneID"])
})

write.table(x = data.frame(urlPath), file = paste(outputFilesDir, "urlPaths.tsv", sep = "/"), sep = "\t", quote = F)

write.table(x = ego@result[,c(2,8)], file = paste(outputFilesDir, "enrichGOMF.tsv", sep = "/"), sep = "\t", quote = F)

write.table(x = ego2@result[,c(2,8)], file = paste(outputFilesDir, "enrichGOCC.tsv", sep = "/"), sep = "\t", quote = F)

write.table(x = ego3@result[,c(2,8)], file = paste(outputFilesDir, "enrichGOBP.tsv", sep = "/"), sep = "\t", quote = F)
