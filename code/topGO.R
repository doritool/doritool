#!/usr/bin/R

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



#################### Install packages if needed ##################
# in the terminal (needed for RCurl and FGNet respectively)
# sudo apt-get install libcurl4-gnutls-dev
# sudo apt-get install libxml2-dev

## in case libcurl4-gnutls-dev gives problems, use libcurl4-openssl-dev

# from CRAN
#install.packages("RCurl")   
#install.packages("gProfileR")
#install.packages("XML")

# from Bioconductor
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
#biocLite("FGNet")
#biocLite("AnnotationDbi")
#biocLite("topGO")
#biocLite("KEGGprofile")

# This script takes 5 minutes to run completely: Rscript topGO_DAVID_gProfiler.

# Load required libraries but suppressing the Package Startup Messages (this is in case you want to embed this inside an automated script) to avoid flooding the terminal with non interesting messages.
suppressPackageStartupMessages(library("FGNet"))
suppressPackageStartupMessages(library("AnnotationDbi"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("topGO"))
suppressPackageStartupMessages(library("KEGGprofile"))
suppressPackageStartupMessages(library("gProfileR"))
library("FGNet")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("topGO")	
library("KEGGprofile")
library("gProfileR")

rm(list=ls())

#1.FGNet   ############################################################
#Perform a Functional Enrichment Analysis (FEA) on a list of genes or expression set
#https://www.bioconductor.org/packages/devel/bioc/vignettes/FGNet/inst/doc/FGNet.html

#1.1 TopGO ############################################################
#TopGO package provides tools for testing GO terms while accounting for the topology of the GO graph.

#Read and create a data frame for ENS gene list
#result=tryCatch, allows to continue working although the outcome is an error. 
#numeric. Minimum size of GO terms. TopGo authors recommend 5-10 for more stable results, 1 for no prune. pValThr = 0.01

ensembldata <- scan("../output/VEP/genelistENS.txt", what="", sep="\n") #ENS gene list   

Result <- tryCatch({
  	feaResults_topGO <- fea_topGO(ensembldata, geneIdType = "ENSEMBL", annotations = c("GO_BP","GO_MF","GO_CC"), organism ="Hs", jobName= "topGO", nodeSize=1, pValThr = 0.01)
	}, error = function(e) 
	{
  "Error in match.names(clabs, names(xi)) : 
    names do not match previous names"
	})

#1.2 ensg_genenames description    
#labels <- as.vector(unlist(as.list(org.Hs.egENSEMBL)))
#Write table ensg_genenames (description of each gene)
ensg_genenames <- as.data.frame(mapIds(org.Hs.eg.db, keys=ensembldata, column="GENENAME", keytype="ENSEMBL"))

write.table(ensg_genenames, file = "./topGO/ensg_genenames.txt",sep="\t", quote=FALSE, col.names = FALSE )
