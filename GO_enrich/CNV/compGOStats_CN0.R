# read in the table
CN0 <-read.table("CN0_genes.tab",header=T,sep="\t",stringsAsFactors=F, quote="")
genes <- CN0$GENE

allgenes <- read.csv("BD_trans_to_gene.tab",
	 header=F,stringsAsFactors=F,sep=" ",quote="")
universe <- allgenes$V2

#universe
#genes

# problem matching mode of this before
mode(universe)
mode(genes)

library("AnnotationDbi")
godat <- read.table("JEL423.IPR.GO",header=F);
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
write.csv(summary(OverMF),"CN0.OverMF_enrich.csv");

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
write.csv(summary(OverCC),"CN0.OverCC_enrich.csv");

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
write.csv(summary(OverBP),"CN0.OverBP_enrich.csv");
