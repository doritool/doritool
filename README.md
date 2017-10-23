# DoriTool - https://doritool.github.io/

DoriTool modules are written in Perl, Bash, and R scripts, and run on any UNIX-like operating system.
To use  DoriTool, simply download the execution script of a docker container, which automatically download the necessary container image with the packages and annotation files. Alternatively, the user can install its dependencies manually (See in detail the [Dockerfile](https://github.com/doritool/doritool/blob/master/Dockerfile).) The full source code and the container is freely available on the GitHub repository (https://github.com/doritool/doritool).
DoriTool website (https://doritool.github.io/) includes general information about the purpose of the tool, instances and explanations about its uses, as well as the link to connect to the GitHub repository in order to download the tool directly.
The DoriTool is available under a GNU GPLv3 license. (https://www.gnu.org/licenses/gpl-3.0.txt).

## Steps
1. Install Docker from the official web site
  - [Linux](https://docs.docker.com/engine/installation/#supported-platforms)
  - [Mac](https://docs.docker.com/docker-for-mac/install/#download-docker-for-mac)
  - [Window](https://docs.docker.com/docker-for-windows/install/) (In case it is needed, use _[docker toolbox](https://www.docker.com/products/docker-toolbox)_ )
       - Notice: Activate BIOS VT-X/AMD-v if is not enabled. Enabling it in the BIOS is mandatory

2. Download the DoriTool repository
    - `git clone https://github.com/doritool/doritool.git`

     or the zip file

    - <https://github.com/doritool/doritool/archive/master.zip>

    Notice:

    For window's users download the DoriTool repository using cmd (click windows R and write cmd, the shell will be opened). Or directly from the docker toolbox Shell.

3. See the next section for learning how to run the _doritool_ script (<font color="red">the first time, the docker image will be downloaded, be patient</font>)

# Quick start

1) Go to GitHub Repository at <https://github.com/doritool/doritool> and read the [Readme.md](https://github.com/doritool/doritool/blob/master/README.md) file. <font color="red">DoriTool uses _GRCh37_ human assembly</font>.

**Read [TROUBLESHOOTING.md](https://github.com/doritool/doritool/blob/master/TROUBLESHOOTING.md) file in case you are not able to run DoriTool**

2) Input data a mutation/variant call format file (VCF) or an rs identifier SNP list.

3) Perform the functional in silico analysis from the shell as follows bellow

- Linux users

    `./doritool`

- Windows users

    `doritool.bat`

Run one of these scripts with the next parameters

3.1 Variant annotation analyslis (without Linkage Disequilibrium and eQTLs).

There are two ways to specify the input file (`-- file`, `-i file` )

`-i web.rs`

`-- web.rs`

3.2. Variant annotation analysis with linkage Disequilibrium (--LD , -l). Notice the cutoff must be specified. LD needs input rs.

`-l 0.90 -i web.rs`

`--LD 0.90 -i web.rs`

3.3. Variant annotation analysis with eQTLs (`--GTEx`, `-e`). Notice that the specific human tissue required must be downloaded previously from https://www.gtexportal.org/home/

    --GTEx Brain_Caudate_basal_ganglia_Analysis.nominal.filtered.txt -i web.rs


___

    -e Brain_Caudate_basal_ganglia_Analysis.nominal.filtered.txt -i web.rs

3.4.Variant annotation analysis with Linkage Disequilibrium and eQTLs

    -l 0.90 -e Brain_Nucleus_accumbens_basal_ganglia_Analysis.nominal.filtered.txt -i web.rs.

# Features considered for DoriTool
Annotation levels and the bioinformatics tools used in DoriTool are shown in the table below (Table1)

#### Table1
Annotation Level | Tool | Input | Use | Script language
---|---|---|---|---
Variant annotation |VEP | SNPs (VCF, CSV) | Variant effect prediction | Perl
| | GTEx | SNPs (LD proxies) | eQTLs | Bash
| |	LD proxy | SNPs rs id	| variants in LD with any overlapping existing variants from the Ensembl variation databases | Perl
Gene annotation	| g:Profiler | Symbol Genes | Functional Analysis | R
| | org<span>.</span>Hs<span>.</span>eg<span>.</span>db |ENS Genes | Annotation | R
| | TopGO (FGNet) | ENS Genes | Functional Annotation | R
Pathway level	| org<span>.</span>Hs<span>.</span>eg<span>.</span>db | ENS Genes	| Annotation | R
| | DAVID(FGNet) | Symbol Genes | Cluster </br> Colour pathways	| R
| | KEGGprofile | Kegg ID | Colour Pathways	| R
Network level	| DAVID (FGNet) | Symbol Genes | Cluster </br> Colour pathways | R


### Mutation/variant level
Variant annotation. Input variants are functionally annotated using VEP, which is an open source, free to use, actively maintained and developed toolset (https://github.com/Ensembl/ensembl-vep) for the analysis, annotation, and prioritization of genomic variants in both coding and non-coding regions. VEP allows determining the effect of the variants (SNPs) on genes, transcripts and proteins, as well as regulatory regions, non-genic variants, and including transcription factor binding sites (TFBS), using regularly updated data files that are distributed by Ensembl and its output follows a standard form (VCF). The plugins and default parameters of VEP considered in DoriTool are shown in the following table.
#### Table2
Plugin options | Features	| Description
---|---|---
--symbol | | Adds the gene symbol to the output
--plugin UpDownDistance | 10000, 5000 | Changes the distance to call upstream and downstream consequence types
--plugin NearestGene | Limit=1 </br> --Max_range=1000000 | Finds the nearest gene to an intergenic variant
--plugin Condel* | | SNP effect prediction

*Condel score (a consensus score considering SIFT and PolyPhen-2, which ranges 0-1, being 0 neutral and 1 deleterious).


DoriTool also locates each variant in its corresponding cytoband, retrieving the information from the database downloaded from UCSC hg19 (http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz).
Expression QTLs. The user has the option to provide a database including tissue-specific significant SNP- gene pairs (with the same format as the ones GTEx provides) to obtain the effect size of the eQTL as well as the eGene.
The specific human tissue required must be downloaded previously from https://www.gtexportal.org/home/
Linkage Disequilibrium proxies. DoriTool also explores proxy and putatively functional SNPs for a query SNP in a selected 1000GP population, using the Ensembl REST API in an integrated Perl script. It computes and returns pairwise LD values 1) among the input variants and 2) between each input variant and all other variants in a surrounding window. The default parameters in DoriTool for window size and strength of LD (r2) are 500 Kb (SNPs for which LD with other variants in a 500 Kb window extracted in EUR population; includes CEU, IBS, TSI, FIN, GBR) and 0.90.

### Gene level
Gene annotation. DoriTool uses three different tools to perform the gene annotation, and therefore allowing interpreting and identifying the biological processes for the gene list tagged by the input list of mutations/variants.
1) org.Hs.eg.db (org.Hs.eg.db: Genome wide annotation for Human. R package version 3.4.1)
2) FGNet, which uses TopGO R package and
3) g:Profiler R package.

DoriTool retrieves the following gene information to annotate: Symbol, gene name, OMIM identifier (Online Mendelian Inheritance in Man) and focuses primarily on inherited or heritable genetic diseases), Gene Ontology identifiers, protein families (PFAM) and Pathway.
topGO functions are included in FGNet R package (https://bioconductor.org/packages/release/bioc/html/topGO.html) to facilitate semi-automated enrichment analysis for Gene Ontology (GO) terms, mapping the genes tagged by the input variants with the associated biological annotation terms (e.g. GO Terms), and then statistically examine the enrichment of gene members for each of the annotation terms on the basis of gene counts. The default parameters in DoriTool are nodeSize=1 (no prune) and pValThr=0.01.
g:Profiler is an open source, free to use, R package actively maintained on CRAN (https://CRAN.R-project.org/package=gProfileR) that performs Functional Enrichment Analysis, including transcription factor binding site predictions, Mendelian disease annotations, information about protein expression and complexes, statistically significant Gene Ontology terms, pathways and other gene function related terms. Multiple testing correction is performed by selecting Benjamini-Hochberg (“FDR”).

### Pathway level
Pathway annotation. It is performed using three R packages: FGNet, which uses DAVID, org.Hs.eg.db, which maps Entrez Gene identifiers to KEGG pathways https://www.bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf. Mappings were based on data provided by: KEGG GENOME ftp://ftp.genome.jp/pub/kegg/genomes.
For the DoriTool pipeline, we included a code to obtain coloured KEGG pathways considering the number of variants per gene by using KEGGprofile.
Colour legend for description of variants affecting genes in networks and pathways are shown below.
#### Table 3

Number of variants|	Colorcode FGNet	   |	KEGG colour    |  DAVID colour|
---|---|---|---
No variants	|	| | |
1	|	NA| <font color="blue">Blue</font> | <font color="green">Green</font>
2 |	-1|	<font color="green">Green</font> | White
3	|	0	|	<font color="#D0D00">Yellow</font> | <font color="red">Red</font>
|>3|	1	|	<font color="red">Red</font> | <font color="red">Red</font>

Moreover, the offline mode of DAVID is considered, since it does not need specific considerations for the installation and guarantees obtaining always results even if the server is down. However, the offline mode does not allow to modify the default settings (default cutoff Linkage: 0.5 and overlap=4, initialseed=4, finalseed=4), resulting in a lower number of KEGG pathways and GO terms, which is compensated with the org.Hs.eg.db annotation.

### Network level
DoriTool also provides coloured Functional Networks (see Table 3) of the list of genes tagged by the input variant list. These functional connections between the different genes were based on annotations (GO) and given by DAVID functions included in FGNet. Building functional networks provids an overview of the biological functions of the genes/terms, and permits links between genes, overlapping between clusters.
Genes in solid colour have a unique cluster (and their color corresponds to the one shown in the table 3), while genes with a white background are shared genes between clusters. The intersection network,  is a simplified functional network where all the genes that belong to only one metagroup are clustered into a single node. Clusters are represented by square boxes and genes are circles.

# Please cite us
When using Doritool for a publication please cite as

**Martín-Antoniano,I; Alonso,L; Madrid,M; López de Maturana, E; Malats,N**. _DoriTool: A bioinformatics integrative tool for post-association functional annotation_. Public Health Genomics. 2017;20(2):126-135. doi: 10.1159/000477561 (<https://doi.org/10.1159/000477561>).

