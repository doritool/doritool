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


#####################GPROFILER########################################
#install.packages("gProfileR")
#upload library
library(gProfileR)
library("org.Hs.eg.db")
#correction_method = "fdr"
#upload genelist and read

#setwd("/home/bioinfo/Escritorio/doritool_pruebas/output")
initialDirectory = "../output"

genelist.data <- read.table("./VEP/genelist.txt") #Symbolgene list

genelist.data <- unlist(as.list(levels(genelist.data$V1)))

geneLabels <- unlist(as.list(org.Hs.egSYMBOL))


Result = tryCatch({
gprofilerRESULTS <- gprofiler(genelist.data, png_fn = "gprofilerResults.png" , significant = T, correction_method = "fdr", organism = "hsapiens")
}, warning = function(w) {}, error = function(e) {})


########### FGnet ############
#Build and visualize functional gene and term networks from clustering of enrichment analyses in multiple annotation spaces.

#upload library
library(FGNet)
#library("AnnotationDbi")
# Load gene symbols list and ENSG (Ensemble) genes
# geneList <- scan("./genes/genelist.txt", what="", sep="\n")
ensemblGeneList <- scan("./VEP/genelistENS.txt", what="", sep="\n")
ensemblGeneList

# import Variants per Gene file: Input: list of genes (Ensembl) with the number of variants and colorcode.
variants_per_gene <- read.table(file="./VEP/variants_per_gene2.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, row.names=1)
#variants_per_gene

# Load gene annotation for symbols
geneLabels <- unlist(as.list(org.Hs.egENSEMBL))

# Read and create a data frame for symbol genelist.txt
geneList <- geneLabels[which(geneLabels %in% ensemblGeneList)]
#geneList

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
#list_variantspergene


###################### (1) Recover gene descriptions and annotations #######################
#setwd(paste0(initialDirectory)
ensg_symbol <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="SYMBOL", keytype="ENSEMBL"))
colnames(ensg_symbol) = "symbol"
ensg_genename <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="GENENAME", keytype="ENSEMBL"))
colnames(ensg_genename) = "genename"
ensg_OMIM <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="OMIM", keytype="ENSEMBL"))
colnames(ensg_OMIM) = "OMIM"
ensg_GO <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="GO", keytype="ENSEMBL"))
colnames(ensg_GO) = "GO"
ensg_PFAM <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="PFAM", keytype="ENSEMBL"))
colnames(ensg_PFAM) = "PFAM"
ensg_pathway <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensemblGeneList, column="PATH", keytype="ENSEMBL"))
colnames(ensg_pathway) = "pathway"
ENSGenes <- cbind(ensg_symbol, ensg_genename, ensg_OMIM, ensg_GO, ensg_PFAM, ensg_pathway, variants_per_gene)

rm(ensg_symbol, ensg_genename, ensg_OMIM, ensg_GO, ensg_PFAM, ensg_pathway)

# org.Hs.egENSEMBL is an R object that provides mappings between entrez gene identifiers and gene abbreviations.
# Get the entrez gene identifiers that are mapped to an ENSEMBL gene
ENSEMBL2EG <- org.Hs.egENSEMBL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_ENSG <- mappedkeys(ENSEMBL2EG)
# Convert to a list
ENSEMBL2EG <- as.list(ENSEMBL2EG[mapped_ENSG])
#ENSEMBL2EG[ensemblGeneList]

ENSGenes[geneList,"ENTREZ"]=names(geneList)

write.table(ENSGenes, file="./genes_annotations.tsv",sep="\t", quote=FALSE, col.names=NA)
#ls()



######### 1.1 David ############################


#source("https://bioconductor.org/biocLite.R")
#biocLite("RDAVIDWebService")

# work with symbol genelist
#Funtional gene network. This tryCatch will perform the fea_david analysis the same way, because the parameters are the same (overlap = 1L, initialSeed = 1L, finalSeed = 1L). At the beginning, in order not to have a very high number of KEGG pathways, we had a more restrictive condition, where overlap = 2L, initialSeed = 2L, finalSeed = 2L. In this case, performs first the restrictive condition, and only if it fails, it launches the permisive one (1L). But when we launched only 1 variant, if its genes do not clusterize with others, no KEGG is given. Since we want to exhaustively annotate the pathways in which a gene participates, we have tuned the parameters to be the same (permisive), but keep the tryCatch original function in case we want to go back to the original behaviour.

geneListSymbol = geneList
geneListSymbol <- as.character(ENSGenes[geneList,"symbol"])
geneListSymbol <- setNames(geneListSymbol,names(geneList))

list_variantspergene <- setNames(list_variantspergene,ENSGenes[names(list_variantspergene),"symbol" ] )

#Result = tryCatch ({
feaResults_DAVID <- fea_david(names(geneListSymbol), geneIdType="ENTREZ_GENE_ID", geneLabels=geneListSymbol, annotations=c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO"), argsWS=c(overlap=1L, initialSeed=1L, finalSeed=1L, linkage=0.5, kappa=35L), jobName="DAVID", email=NULL, downloadFile=TRUE)
#}, warning = function(w) {}, error = function(e) {} )

#Result = tryCatch ({
FGNet_report(feaResults_DAVID, geneExpr=list_variantspergene, plotKeggPw=TRUE, plotGoTree=FALSE)
#}, warning = function(w) {}, error = function(e) {} )

#Color Kegg pathways from FeaResults_DAVID and annotation ID

dir.create("./KEGG")
setwd("./KEGG")

#1 annotation

Result = tryCatch ({
levels(ENSGenes$pathway)
for (keggNumber in levels(ENSGenes$pathway)) {
  plotKegg(paste0("hsa",keggNumber), geneExpr=list_variantspergene, geneIDtype="ENSEMBL", colType= "discrete")
}
}, warning = function(w) {}, error = function(e) {} )

#2 david

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



# Add all terms in annotation.tsv file and incorporate them into completeReport.tsv.

setwd("../../")

ENSGenes$pathway <- sapply(ENSGenes$pathway, function(path) {ifelse(is.na(path), NA, paste0("hsa:", path))})
ens <- ENSGenes[, c("symbol", "pathway")]
ens <- ens[!is.na(ens$symbol),]
rownames(ens) <- NULL

DAVID_formatted <- read.csv(file="./DAVID/DAVID_formatted.txt", sep="\t", header=TRUE)

# Take all KEGG pathways from DAVID_formatted txt
d <- DAVID_formatted[grep("KEGG", DAVID_formatted$Terms), c("Terms", "Genes")]
d$Terms <-sub("KEGG:","",d$Terms)
d$Genes <- as.character(d$Genes)
ens$symbol <- as.character(ens$symbol)
df <- sapply(ens$symbol, function(symbol){sapply(d$Genes, function(genes){grep(symbol, genes)})})
df_t <- apply(df, 1, function(x){x==1})
#d$Terms[!is.na(df_t[1,])]
david_pathways <- sapply(rownames(df_t), function(symbol){paste(d$Terms[!is.na(df_t[symbol,])], collapse=";")})
david_pathways <- data.frame(david_pathways)
david_pathways$symbol <- rownames(david_pathways)
rownames(david_pathways) <- NULL
all_pathways <- merge(ens, david_pathways, by="symbol")
all_pathways$final_pathways <- paste(all_pathways$pathway,all_pathways$david_pathways, sep=";")
all_pathways$final_pathways <- sub("NA;", "", all_pathways$final_pathways)
all_pathways$final_pathways <- sub("^$", NA, all_pathways$final_pathways)
all_pathways$final_pathways <- sub(";$", "", all_pathways$final_pathways)
all_pathways$pathway <- NULL
all_pathways$david_pathways <- NULL
final_table <- merge(ENSGenes, all_pathways, by="symbol")
final_table$genename <- NULL
final_table$GO <- NULL
final_table$pathway <- NULL
final_table$ENTREZ <- NULL
final_table$variants <- NULL

write.table(final_table, file="./genes_annotations2.tsv",sep="\t", na = "-", quote=FALSE, col.names=NA)

#list_variantspergene
##keggIDs_DAVID = gsub("hsa","",as.vector(unlist(keggIDs_DAVID)))
##keggIDs_annotation = levels(ENSGenes$pathway)
##keggIDs = unique(c(keggIDs_DAVID,keggIDs_annotation))

##AnnotationDbi::as.list(KEGGPATHID2NAME)[keggIDs] # KEGG IDS NAMES, hay que juntarla
##feaResults_DAVID$geneTermSets[feaResults_DAVID$geneTermSets$Category=="KEGG_PATHWAY",]   # KEGG IDS NAMES, hay que juntarla



#keggIDs = paste0("hsa",keggIDs)
#Result = tryCatch ( {
#    for (keggID in keggIDs) {
#      plotKegg(keggID, geneExpr=list_variantspergene, geneIDtype="SYMBOL", colType= "discrete")
#    }
#  }, warning = function(w) {},error = function(e) {} )
#keggIDs
#library("KEGG.db")

