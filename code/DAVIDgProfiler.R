# DoriTool is the most complete bioinformatic tool offering functional 'in silico' annotation of variants
# Copyright (C) 2017 Isabel Martín-Antoniano, Lola Alonso,Miguel Madrid; Evangelina López de Maturana; Núria Malats. Genetic and Molecular Epidemiology Group of Spanish National Cancer Research Centre (CNIO)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


##################### gProfileR  ########################################
suppressPackageStartupMessages(library("gProfileR"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("FGNet"))
library("gProfileR")
library("org.Hs.eg.db")
library("AnnotationDbi")
library("FGNet")

setwd(Sys.getenv("OUTPUTDIR"))
# load genelist
genelist.data <- read.table("./VEP/genelist.txt") #Symbolgene list
genelist.data <- unlist(as.list(levels(genelist.data$V1)))
geneLabels <- unlist(as.list(org.Hs.egSYMBOL))

Result = tryCatch({
   gprofilerRESULTS <- gprofiler(genelist.data, png_fn="gprofilerResults.png", significant=T, correction_method="fdr", organism="hsapiens")
}, warning=function(w) {}, error=function(e) {})

# Load gene symbols list and ENSG (Ensemble) genes
ensemblGeneList <- scan("./VEP/genelistENS.txt", what="", sep="\n")

# import Variants per Gene file: list of genes (Ensembl) with the number of variants to create colorcode.
variants_per_gene <- read.table(file="./VEP/variants_per_gene2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
variants_per_gene_tmp = variants_per_gene

# Load gene annotation for symbols
geneLabels <- unlist(as.list(org.Hs.egENSEMBL))

# Read and create a data frame for symbol genelist.txt
geneList <- geneLabels[which(geneLabels %in% ensemblGeneList)]

# Recode the variants_per_gene list to colour genes by colorcode
# 1 variant per gene --> (NA = "blue")
# 2 variants per gene --> (-1 = "green")
# 3 variants per gene --> (0 = "yellow")
# >3 variants per gene --> (1 = "red")
variants_per_gene$variants[variants_per_gene$variants == 1] = NA
variants_per_gene$variants[variants_per_gene$variants == 2]= -1
variants_per_gene$variants[variants_per_gene$variants == 3] = 0
variants_per_gene$variants[variants_per_gene$variants > 3] = 1
variants_per_gene

list_variantspergene <- setNames(variants_per_gene$variants, as.list(rownames(variants_per_gene)))

###################### Recover gene descriptions and annotations #######################
ensg_symbol <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="SYMBOL", keytype="ENSEMBL")))
colnames(ensg_symbol) = "symbol"
ensg_genename <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="GENENAME", keytype="ENSEMBL")))
colnames(ensg_genename) = "genename"
ensg_OMIM <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="OMIM", keytype="ENSEMBL")))
colnames(ensg_OMIM) = "OMIM"
ensg_GO <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="GO", keytype="ENSEMBL")))
colnames(ensg_GO) = "GO"
ensg_PFAM <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="PFAM", keytype="ENSEMBL")))
colnames(ensg_PFAM) = "PFAM"
ensg_pathway <- as.data.frame(as.character(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="PATH", keytype="ENSEMBL")))
colnames(ensg_pathway) = "pathway"
ensg_pathway$pathway <- gsub("^","hsa",ensg_pathway$pathway)
ENSGenes <- cbind(ensg_symbol, ensg_genename, ensg_OMIM, ensg_GO, ensg_PFAM, ensg_pathway, variants_per_gene_tmp)

rm(ensg_symbol, ensg_genename, ensg_OMIM, ensg_GO, ensg_PFAM, ensg_pathway, variants_per_gene_tmp)

# org.Hs.egENSEMBL is an R object that provides mappings between entrez gene identifiers and gene abbreviations.
ENSEMBL2EG <- org.Hs.egENSEMBL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_ENSG <- mappedkeys(ENSEMBL2EG)
# Convert to a list
ENSEMBL2EG <- as.list(ENSEMBL2EG[mapped_ENSG])
ENSGenes[geneList,"ENTREZ"]=names(geneList)

write.table(ENSGenes, file="./genes_annotations.tsv",sep="\t", quote=FALSE, col.names=NA)

########### FGNet ############
#Build and visualize functional gene and term networks from clustering of enrichment analyses in multiple annotation spaces, using DAVID

geneListSymbol = geneList
geneListSymbol <- as.character(ENSGenes[geneList,"symbol"])
geneListSymbol <- setNames(geneListSymbol,names(geneList))

list_variantspergene <- setNames(list_variantspergene,ENSGenes[names(list_variantspergene),"symbol" ] )

Result = tryCatch ({
   feaResults_DAVID <- fea_david(names(geneListSymbol), geneIdType="ENTREZ_GENE_ID", geneLabels=geneListSymbol, annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), argsWS=c(overlap=1L, initialSeed=1L, finalSeed=1L, linkage=0.5, kappa=35L), jobName="DAVID", email=NULL, downloadFile=TRUE)
}, warning = function(w) {}, error = function(e) {} )

Result = tryCatch ({
   FGNet_report(feaResults_DAVID, geneExpr=list_variantspergene, plotKeggPw=TRUE, plotGoTree=FALSE)
}, warning = function(w) {}, error = function(e) {} )

# Color Kegg pathways from FeaResults_DAVID and annotation ID
dir.create("./KEGG")
setwd("./KEGG")

# 1- annotation
Result = tryCatch ({
  for (keggNumber in levels(as.factor(ENSGenes$pathway))) {
      plotKegg(keggNumber, geneExpr=list_variantspergene, geneIDtype="ENSEMBL", colType= "discrete")
  }
}, warning = function(w) {}, error = function(e) {} )

# 2- DAVID
dir.create("./KEGGDAVID")
setwd("./KEGGDAVID")

keggIDs=list()
all_keggIDs=list()

Result = tryCatch ( {
  for (i in 1:length(getTerms(feaResults_DAVID, returnValue="KEGG"))) {
    all_keggIDs = getTerms(feaResults_DAVID, returnValue="KEGG")[[i]]
    keggIDs=c(keggIDs, all_keggIDs[all_keggIDs != "named list()"] )
  }
  }, warning = function(w) {}, error = function(e) {} )

Result = tryCatch ( {
for (keggID in keggIDs) {
  plotKegg(keggID, geneExpr=list_variantspergene, geneIDtype="SYMBOL", colType= "discrete")
}
}, warning = function(w) {}, error = function(e) {} )
