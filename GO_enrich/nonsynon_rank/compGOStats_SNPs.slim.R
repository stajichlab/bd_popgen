# read in the table
nonSynSNPtable <-read.table("49strains_Nonsynonymous_UM142_anno_top10quan.txt",header=T,sep="\t",stringsAsFactors=F, quote="")
# get ride of genes which empy lines
genes <- subset(nonSynSNPtable$Broad_Gene_ID,nonSynSNPtable$Broad_Gene_ID != 'NA')

#get the list of all genes
allgenes <- read.csv("BD_trans_to_gene.tab",
	 header=F,stringsAsFactors=F,sep=" ",quote="")
# just want the names of the universe of possible genes
universe <- unique(allgenes$V2)


#universe
#genes

# problem matching mode of this before
mode(universe)
mode(genes)

library("AnnotationDbi")
# define the table of GO terms for genes from this tab delimited file
godat <- read.table("JEL423.IPR.GO_slim",header=F,sep="\t");
goframeData <- data.frame(godat$V1, godat$V2, godat$V3)
goFrame <- GOFrame(goframeData,organism="Batrachochytrium dendrobatidis")
goAllFrame=GOAllFrame(goFrame)

library("GSEABase")
library("GOstats")

gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "MF",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")


OverMF <- hyperGTest(params)
summary(OverMF)
OverMF
write.csv(summary(OverMF),"OverMF_enrich.slim.csv");

paramsCC <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "CC",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")

OverCC <- hyperGTest(paramsCC)
summary(OverCC)
OverCC
write.csv(summary(OverCC),"SNP.OverCC_enrich.slim.csv");
paramsBP <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
	  geneIds = genes,
	  universeGeneIds = universe,
	  ontology = "BP",
	  pvalueCutoff = 0.05,
	  conditional = FALSE,
	  testDirection = "over")

OverBP <- hyperGTest(paramsBP)
summary(OverBP)
OverBP
write.csv(summary(OverBP),"SNP.OverBP_enrich.slim.csv");
# you can look for underrepresented too with the 'testDirection="under"'

paramsCC <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
          geneIds = genes,
          universeGeneIds = universe,
          ontology = "CC",
          pvalueCutoff = 0.05,
          conditional = FALSE,
          testDirection = "under")

UnderCC <- hyperGTest(paramsCC)
summary(UnderCC)
UnderCC
write.csv(summary(UnderCC),"SNP.UnderCC_enrich.slim.csv");

paramsBP <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
          geneIds = genes,
          universeGeneIds = universe,
          ontology = "BP",
          pvalueCutoff = 0.05,
          conditional = FALSE,
          testDirection = "under")

UnderBP <- hyperGTest(paramsBP)
summary(UnderBP)
UnderBP
write.csv(summary(UnderBP),"SNP.UnderBP_enrich.slim.csv");

params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
          geneIds = genes,
          universeGeneIds = universe,
          ontology = "MF",
          pvalueCutoff = 0.05,
          conditional = FALSE,
          testDirection = "under")


Under <- hyperGTest(params)
summary(UnderMF)
UnderMF
write.csv(summary(UnderMF),"SNP.UnderMF_enrich.slim.csv");
