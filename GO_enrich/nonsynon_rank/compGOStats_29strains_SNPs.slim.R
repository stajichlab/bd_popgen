library("GSEABase")
library("GOstats")
library("AnnotationDbi")

# read in the table
nonSynSNPtable <-read.table("29strains_dnds_UM142_over1.txt",header=F,sep="\t",stringsAsFactors=F, quote="")
# get ride of genes which empy lines
genes <- nonSynSNPtable$V1

#get the list of all genes
allgenes <- read.csv("BD_trans_to_gene.tab",
	 header=F,stringsAsFactors=F,sep=" ",quote="")
# just want the names of the universe of possible genes
universe <- allgenes$V2

#universe
#genes

# problem matching mode of this before
mode(universe)
mode(genes)

# define the table of GO terms for genes from this tab delimited file
godat <- read.table("JEL423.INTERPRO.go_slim",header=F,sep="\t");
goframeData <- data.frame(godat$V1, godat$V2, godat$V3)
goFrame <- GOFrame(goframeData,organism="Batrachochytrium dendrobatidis")
goAllFrame=GOAllFrame(goFrame)


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
write.csv(summary(OverMF),"OverMF_enrich_29strains.slim.csv");

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
write.csv(summary(OverCC),"OverCC_enrich_29strains.slim.csv");
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
write.csv(summary(OverBP),"OverBP_enrich_29strains.slim.csv");
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
write.csv(summary(UnderCC),"UnderCC_enrich_29strains.slim.csv");

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
write.csv(summary(UnderBP),"UnderBP_enrich_29strains.slim.csv");

params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
          geneIds = genes,
          universeGeneIds = universe,
          ontology = "MF",
          pvalueCutoff = 0.05,
          conditional = FALSE,
          testDirection = "under")


UnderMF <- hyperGTest(params)
summary(UnderMF)
write.csv(summary(UnderMF),"UnderMF_enrich_29strains.slim.csv");
