# read in the table
CN0 <-read.table("loh_10k_strict_filter25.lst",header=F,sep="\t",stringsAsFactors=F, quote="")
genes <- CN0$V1

allgenes <- read.csv("BD_trans_to_gene.tab",
	 header=F,stringsAsFactors=F,sep=" ",quote="")
universe <- unique(allgenes$V2)

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


Over <- hyperGTest(params)
summary(Over)
Over
write.csv(summary(Over),"OverMF_enrich.csv");
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
write.csv(summary(OverCC),"OverCC_enrich.csv");

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
write.csv(summary(OverBP),"OverBP_enrich.csv");


params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
          geneSetCollection=gsc,
          geneIds = genes,
          universeGeneIds = universe,
          ontology = "MF",
          pvalueCutoff = 0.05,
          conditional = FALSE,
          testDirection = "under")


Under <- hyperGTest(params)
summary(Under)
Under
write.csv(summary(Under),"UnderMF_enrich.csv");


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
write.csv(summary(UnderCC),"UnderCC_enrich.csv");

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
write.csv(summary(UnderBP),"UnderBP_enrich.csv");
