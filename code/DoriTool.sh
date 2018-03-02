#!/bin/bash

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


#set -e    #para interrumpir script en caso de salir errores o faltar archivo
# create work directory and paths
BASEDIR="$(dirname $(dirname "$0"))"
echo "BASEDIR = ${BASEDIR}"
export WORKDIR_CODE="$BASEDIR"/code
export WORKDIR_DATABASES="$BASEDIR"/databases
export INPUTDIR="$PWD"/input
export OUTPUTDIR="$PWD"/output

LDrequired=
GTExrequired=
function usage
{
    echo "usage: Please provide following arguments --input inputfile.rs (or VCF) --LD LD_value (optional) --GTEx filename_with_eQTLs (optional)"
}

while [ "$1" != "" ]; do
    case $1 in
	  -i | --input )		shift
					FILEINPUT=$1
					;;
        -l | --LD )           shift
					LDrequired=$1
                                ;;
        -e | --GTEx ) 		shift
					filenameGTEx=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
    esac
    shift
done

mkdir "$INPUTDIR"
mkdir "$OUTPUTDIR"
mkdir "$OUTPUTDIR"/VEP

# from user_input.csv we create input.vcf following this format
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT

firstline=$(head -1 $FILEINPUT)
thirdcol=$(tail -1 $FILEINPUT | cut -f3)

if [[ $FILEINPUT =~ .*\.(VCF|vcf) ]] ; then             # VCF standard
  cp "$PWD"/"$FILEINPUT" "$INPUTDIR"/"$FILEINPUT"
elif [[ $firstline =~ ^rs ]] ; then                     # rs list one by line
  cp "$PWD"/"$FILEINPUT" "$INPUTDIR"/"$FILEINPUT"
elif [[ $thirdcol =~ ^rs ]] ; then                      # user provides chr, pos, rs, ref, alt
   awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5"\t.\t.\t.\t."} ' "$PWD"/"$FILEINPUT" | sed '1d' > "$INPUTDIR"/"$FILEINPUT".vcf
  sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' "$INPUTDIR"/"$FILEINPUT".vcf
  FILEINPUT="${FILEINPUT}.vcf"
else                                                    # user provides chr, pos, ref, alt
  awk 'BEGIN{OFS="\t"} {print $1,$2,$2,$3,$4"\t.\t.\t.\t."} ' "$BASEDIR"/"$FILEINPUT" | sed '1d' > "$INPUTDIR"/"$FILEINPUT".vcf
  sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' "$INPUTDIR"/"$FILEINPUT".vcf
  FILEINPUT="${FILEINPUT}.vcf"
fi

################ Variant Effect Predictor ###################
vep --database --force_overwrite --port 3337 -i "$INPUTDIR"/"$FILEINPUT" --symbol --tab --distance 10000,5000  --regulatory --plugin NearestGene,limit=1,max_range=1000000 --plugin Condel -o "$OUTPUTDIR"/VEP/condeloutput.txt


if [[ $LDrequired != "" ]] ; then
    "$WORKDIR_CODE"/returnsLD_pairwise.pl "$INPUTDIR"/"$FILEINPUT"
    "$WORKDIR_CODE"/returnsLD_othervariantsEUR.pl "$INPUTDIR"/"$FILEINPUT" "$LDrequired"
fi

if [[ "$filenameGTEx" != "" ]] ; then
   vep --force_overwrite --port 3337 -i "$INPUTDIR"/"$FILEINPUT" --vcf -o "$OUTPUTDIR"/VEP/condeloutput.vcf
   awk 'BEGIN{FS="\t"} { if ( ! /^#/ ) {print $3"\t"$1"_"$2"_"$4"_"$5"_b37"} }' "$OUTPUTDIR"/VEP/condeloutput.vcf > "$OUTPUTDIR"/VEP/unique_variants.txt
   printf  "rs_id\tvariant_id\tgene_id\ttss_distance\tpval_nominal\tslope\tslope_se\tslope_fpkm\tslope_fpkm_se\tpval_nominal_threshold\tmin_pval_nominal\tpval_beta" > "$OUTPUTDIR"/variants_eQTLs.tsv
   for variant in `cut -f2 "$OUTPUTDIR"/VEP/unique_variants.txt`; do printf $variant"\n"; grep $variant -w $filenameGTEx ; done >> "$OUTPUTDIR"/variants_eQTLs.tsv
#    rm "$OUTPUTDIR"/VEP/condeloutput.vcf
fi

# Change condeloutput.txt_summary.html name into Summary_VEP.html
mv "$OUTPUTDIR"/VEP/condeloutput.txt_summary.html "$OUTPUTDIR"/VEP/Summary_VEP.html

############### CytoBand script ########
#1- clean condeloutput.txt file and within gene variable incorporate "intergenic" when "-" occur.
grep -v '^##' "$OUTPUTDIR"/VEP/condeloutput.txt | awk 'BEGIN{FS="\t";OFS="\t"} {if ($4 == "-") {$4="INTERGENIC"} {print $0} }' | sed -e 's/#//g' > "$OUTPUTDIR"/VEP/condeloutputclean0.txt

# Add final column with cytoband into condeloutputclean0.txt. We create condeloutputclean1.txt with the cytoband script. The condeloutputclean1.txt is created by the perl script.
echo "WORKDIR_CODE cyto = ${WORKDIR_CODE}"
"$WORKDIR_CODE"/cytoband_finder.pl "$OUTPUTDIR"/VEP/condeloutputclean0.txt "$WORKDIR_DATABASES"/cytoBand.tsv

#1.1 clean condeloutputclean1.txt and leave it without without header
sed '1d' "$OUTPUTDIR"/VEP/condeloutputclean1.txt > "$OUTPUTDIR"/VEP/condeloutputclean2.txt

#1.2a) genelist with the variant position and gene. Count the number of variants per gene. And Create a file "variants_per_gene1.txt" with the number and ENSG id that will be more robust for DAVID
cut -f2,4 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | sort -u | sort -k2 > "$OUTPUTDIR"/VEP/variants_per_gene1.txt
for gene in $(cut -f2 "$OUTPUTDIR"/VEP/variants_per_gene1.txt | sort -u | grep -w -v 'INTERGENIC' | grep -v 'Location')
do
  printf $gene"\t"; grep $gene -w -c "$OUTPUTDIR"/VEP/variants_per_gene1.txt
done > "$OUTPUTDIR"/VEP/variants_per_gene3.txt

#1.2b) genelist with the variant position and gene. Count the number of variants per gene. And Create a file "variants_per_gene4.txt" with the number and gene SYMBOL, that will be more significant for the user
cut -f2,18 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | sort -u | sort -k2 > "$OUTPUTDIR"/VEP/variants_per_gene4.txt
for gene in $(cut -f2 "$OUTPUTDIR"/VEP/variants_per_gene4.txt | sort -u | grep -w -v '-' | grep -v 'Location')
do
  printf $gene"\t"; grep $gene -c -w "$OUTPUTDIR"/VEP/variants_per_gene4.txt
done > "$OUTPUTDIR"/VEP/variants_per_gene.txt

#Incorporate headline into variants_per_gene1.txt
grep -v "Gene" "$OUTPUTDIR"/VEP/variants_per_gene3.txt > "$OUTPUTDIR"/VEP/variants_per_gene2.txt
sed -i '1i gene\tvariants' "$OUTPUTDIR"/VEP/variants_per_gene2.txt

#1.3 Create a gene list: ENS and Symbol list
cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | grep -v 'INTERGENIC' | sort -u > "$OUTPUTDIR"/VEP/genelistENS.txt
cut -f18 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | grep -v '^-' |sort -u > "$OUTPUTDIR"/VEP/genelist.txt

############# topGO.R and DAVIDgPROFILER.R scripts ##############
(cd "$OUTPUTDIR" && Rscript "$WORKDIR_CODE"/topGO.R)
cd "$OUTPUTDIR"
(cd "$OUTPUTDIR" && Rscript "$WORKDIR_CODE"/DAVIDgProfiler.R)

########## create Complete report.tsv #############################
# Match for every ENSG its ontologies in topGO.txt
for i in $(cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean1.txt ); do printf $i"\n"; grep $i -w "$OUTPUTDIR"/topGO/topGO.txt; done | cut -f2 | awk '/^ENSG/ {printf "\n"$0"\t"; next} {printf $0";"}' | sed 's/INTERGENIC;/\nINTERGENIC/g' | sed '1 s/Gene;/Gene\tGO/g' | sed 's/GO:[0-9]*://g' | sed 's/;/; /g' > "$OUTPUTDIR"/VEP/topaste0.txt

# Match for every SYMBOL (there are not ENSG available in DAVID) its KEGG pathways in DAVID_formatted.txt
for i in `cut -f18 "$OUTPUTDIR"/VEP/condeloutputclean1.txt`; do printf "CADENADECONTROL\n"; grep $i -w "$OUTPUTDIR"/DAVID/DAVID_formatted.txt | grep 'KEGG'; done | cut -f4 | tr '\n' ' ' | sed 's/CADENADECONTROL/\n/g' | sed 's/^ KEGG://' | sed 's/ KEGG:/;/g' | sed '1d' | sed '1d' > "$OUTPUTDIR"/VEP/topaste1.txt

# Match for every ENSG its genenames with the gene names from topGO result file ensg_genenames.txt
for i in $(cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean1.txt ); do printf $i"\t\n"; grep $i -w -m1 "$OUTPUTDIR"/topGO/ensg_genenames.txt; done | awk 'BEGIN{FS="\t"; OFS="\t"} /^INTERGENIC/ {printf $0"this region is not a gene\n"; next}{print $0} ' | cut -f2 | sed '/^$/d' > "$OUTPUTDIR"/topGO/topaste2.txt

# Match for every ENSG its pathways from annotation to write them in another column
for i in $(cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean1.txt ); do printf "CADENADECONTROL\n"; grep $i -w  "$OUTPUTDIR"/genes_annotations.tsv | grep 'hsa'; done | cut -f7 | tr '\n' ' ' | sed 's/CADENADECONTROL/\n/g' | sed 's/^ hsa/hsa/' | sed 's/ hsa/;hsa/g' > "$OUTPUTDIR"/VEP/topaste3.txt

#write "KEGG_DAVID" header in topaste1.txt, "KEGG_annot" in topaste3 and write header "description" in topaste2.txt
sed -i '1i KEGG_DAVID' "$OUTPUTDIR"/VEP/topaste1.txt
sed -i '1d;2d' "$OUTPUTDIR"/VEP/topaste3.txt
sed -i '1i KEGG_annot' "$OUTPUTDIR"/VEP/topaste3.txt
sed -i '1i Gene_Description' "$OUTPUTDIR"/topGO/topaste2.txt

# Paste KEGG_DAVID, KEGG_annot; complete Gene Name, GO, beside condeloutputclean1.txt
paste "$OUTPUTDIR"/VEP/condeloutputclean1.txt "$OUTPUTDIR"/VEP/topaste1.txt > "$OUTPUTDIR"/VEP/condeloutputclean3.txt
paste "$OUTPUTDIR"/VEP/condeloutputclean3.txt "$OUTPUTDIR"/VEP/topaste3.txt > "$OUTPUTDIR"/VEP/condeloutputclean4.txt
paste "$OUTPUTDIR"/VEP/condeloutputclean4.txt "$OUTPUTDIR"/topGO/topaste2.txt > "$OUTPUTDIR"/VEP/condeloutputclean5.txt
paste "$OUTPUTDIR"/VEP/condeloutputclean5.txt "$OUTPUTDIR"/VEP/topaste0.txt > "$OUTPUTDIR"/CompleteReport.tsv

# create SummaryReport with the most relevant fields
cut -f1-4,7,14,17,18,21,22,23,24,25,26,27,28,29,30,32 "$OUTPUTDIR"/CompleteReport.tsv | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' > "$OUTPUTDIR"/SummaryReport0.tsv
# we remove inline the header to avoid considering it in the next sort function
sed -i '1d' "$OUTPUTDIR"/SummaryReport0.tsv

# make the $IFS the newline character
BACKUPIFS=$IFS
IFS=$'\n'
for vargene in `cut -f3,4,5 "$OUTPUTDIR"/SummaryReport0.tsv | sort -u | sort -n`; do grep -w -m1 -P "$vargene\t" "$OUTPUTDIR"/SummaryReport0.tsv; done > "$OUTPUTDIR"/SummaryReport.tsv
# we add the header back again
IFS=$BACKUPIFS

sed -i '1i Uploaded_variation\tAllele\tLocation\tGene\tConsequence\tIMPACT\tFLAGS\tSYMBOL\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tNearestGene\tCondel\tcytoband\tKEGG_DAVID\tKEGG_annot\tGene_Description\tGO' "$OUTPUTDIR"/SummaryReport.tsv

# remove extra intermediary files
# rm "$OUTPUTDIR"/VEP/condeloutput.txt  # we keep this file because the VEP.html links directly to it
rm "$OUTPUTDIR"/VEP/condeloutputclean0.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean1.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean2.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean3.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean4.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean5.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene1.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene2.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene3.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene4.txt
rm "$OUTPUTDIR"/VEP/topaste0.txt
rm "$OUTPUTDIR"/VEP/topaste1.txt
rm "$OUTPUTDIR"/topGO/topaste2.txt
rm "$OUTPUTDIR"/VEP/topaste3.txt
rm "$OUTPUTDIR"/topGO/ensg_genenames.txt
rm "$OUTPUTDIR"/SummaryReport0.tsv
rm "$OUTPUTDIR"/VEP/genelist.txt
rm "$OUTPUTDIR"/VEP/genelistENS.txt

exit
