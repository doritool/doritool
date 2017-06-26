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
BASEDIR=${PWD}
WORKDIR="$BASEDIR"/code
INPUTDIR="$BASEDIR"/input
OUTPUTDIR="$BASEDIR"/output

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
#         * )                     usage
#                                 exit 1
    esac
    shift
done

mkdir "$INPUTDIR"
mkdir "$OUTPUTDIR"
mkdir "$OUTPUTDIR"/VEP

# from user_input.csv we create input.vcf following this format
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT

firstline=`head -1 $FILEINPUT `
thirdcol=`tail -1 $FILEINPUT | cut -f3`

if [[ $FILEINPUT =~ .*\.(VCF|vcf) ]] ; then             # VCF standard
  cp "$BASEDIR"/"$FILEINPUT" "$INPUTDIR"/"$FILEINPUT"
elif [[ $firstline =~ ^rs ]] ; then                     # rs list one by line
  cp "$BASEDIR"/"$FILEINPUT" "$INPUTDIR"/"$FILEINPUT"
elif [[ $thirdcol =~ ^rs ]] ; then                      # user provides chr, pos, rs, ref, alt
   awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5"\t.\t.\t.\t."} ' "$BASEDIR"/"$FILEINPUT" | sed '1d' > "$INPUTDIR"/"$FILEINPUT".vcf
  sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' "$INPUTDIR"/"$FILEINPUT".vcf
  FILEINPUT=`echo $FILEINPUT.vcf`
else                                                    # user provides chr, pos, ref, alt
  awk 'BEGIN{OFS="\t"} {print $1,$2,$2,$3,$4"\t.\t.\t.\t."} ' "$BASEDIR"/"$FILEINPUT" | sed '1d' > "$INPUTDIR"/"$FILEINPUT".vcf
  sed -i '1i #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT' "$INPUTDIR"/"$FILEINPUT".vcf
  FILEINPUT=`echo $FILEINPUT.vcf`
fi

################ Variant Effect Predictor ###################
vep --database --force_overwrite --port 3337 -i "$INPUTDIR"/"$FILEINPUT" --symbol --tab --plugin UpDownDistance,10000,5000 --regulatory --plugin NearestGene,limit=1,max_range=1000000 --plugin Condel -o "$OUTPUTDIR"/VEP/condeloutput.txt


if [[ $LDrequired != "" ]] ; then
    "$WORKDIR"/returnsLD_pairwise.pl "$INPUTDIR"/"$FILEINPUT"
    "$WORKDIR"/returnsLD_othervariantsEUR.pl "$INPUTDIR"/"$FILEINPUT" $LDrequired
fi

if [[ "$filenameGTEx" != "" ]] ; then
   vep --force_overwrite --port 3337 -i "$INPUTDIR"/"$FILEINPUT" --vcf -o "$OUTPUTDIR"/VEP/condeloutput.vcf
   awk 'BEGIN{FS="\t"} { if ( ! /^#/ ) {print $3"\t"$1"_"$2"_"$4"_"$5"_b37"} }' "$OUTPUTDIR"/VEP/condeloutput.vcf > "$OUTPUTDIR"/VEP/unique_variants.txt
   printf  "rs_id\tvariant_id\tgene_id\ttss_distance\tpval_nominal\tslope\tslope_se\tslope_fpkm\tslope_fpkm_se\tpval_nominal_threshold\tmin_pval_nominal\tpval_beta" > "$OUTPUTDIR"/variants_eQTLs.tsv
   for variant in `cut -f2 "$OUTPUTDIR"/VEP/unique_variants.txt`; do printf $variant"\n"; grep $variant -w $filenameGTEx ; done >> "$OUTPUTDIR"/variants_eQTLs.tsv
fi

# Change condeloutput.txt_summary.html name into Summary_VEP.html
mv "$OUTPUTDIR"/VEP/condeloutput.txt_summary.html "$OUTPUTDIR"/VEP/Summary_VEP.html


############### CytoBand script ########

#1- clean condeloutput file and  within gene variable incorporate "intergenic" when "-" occur.
grep -v '^##' "$OUTPUTDIR"/VEP/condeloutput.txt | awk 'BEGIN{FS="\t";OFS="\t"} {if ($4 == "-") {$4="INTERGENIC"} {print $0} }' | sed -e 's/#//g' > "$OUTPUTDIR"/VEP/condeloutputclean0.txt

# Add final column with cytoband into condeloutputclean0.txt . We create condeloutputclean.txt with the cytoband script.
"$WORKDIR"/cytoband_finder.pl "$OUTPUTDIR"/VEP/condeloutputclean0.txt "$WORKDIR"/../databases/cytoBand.tsv

#1.1 clean condeloutputclean0 without header

sed '1d' "$OUTPUTDIR"/VEP/condeloutputclean0.txt  > "$OUTPUTDIR"/VEP/condeloutputclean2.txt

#1.2a) genelist with the variant position and gene. Count the number of variants per gene. And Create a file "variants_per_gene1.txt" with the number and ENSG id that will be more robust for DAVID

cut -f2,4 "$OUTPUTDIR"/VEP/condeloutputclean0.txt | sort -u | grep -v '^#' | sort -k2 > "$OUTPUTDIR"/VEP/variants_per_gene1.txt
for gene in $(cut -f2 "$OUTPUTDIR"/VEP/variants_per_gene1.txt | sort -u | grep -w -v 'INTERGENIC' | grep -v 'Location')
do
  printf $gene"\t"; grep $gene -w -c "$OUTPUTDIR"/VEP/variants_per_gene1.txt
done > "$OUTPUTDIR"/VEP/variants_per_gene2_tmp.txt

#1.2b) genelist with the variant position and gene. Count the number of variants per gene. And Create a file "variants_per_gene3.txt" with the number and gene SYMBOL, that will be more significant for the user
cut -f2,18 "$OUTPUTDIR"/VEP/condeloutputclean0.txt | sort -u | grep -v '^#' | sort -k2 > "$OUTPUTDIR"/VEP/variants_per_gene3.txt
for gene in $(cut -f2 "$OUTPUTDIR"/VEP/variants_per_gene3.txt | sort -u | grep -w -v '-' | grep -v 'Location')
do
  printf $gene"\t"; grep $gene -c -w "$OUTPUTDIR"/VEP/variants_per_gene3.txt
done > "$OUTPUTDIR"/VEP/variants_per_gene.txt

#Incorporate headline into variants_per_gene1.txt
grep -v "Gene" "$OUTPUTDIR"/VEP/variants_per_gene2_tmp.txt > "$OUTPUTDIR"/VEP/variants_per_gene2.txt
sed -i '1i gene\tvariants' "$OUTPUTDIR"/VEP/variants_per_gene2.txt

#1.3 Create a gene list: ENS and Symbol list

cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | grep -v 'INTERGENIC' | sort -u > "$OUTPUTDIR"/VEP/genelistENS.txt
cut -f18 "$OUTPUTDIR"/VEP/condeloutputclean2.txt | grep -v '^-' |sort -u > "$OUTPUTDIR"/VEP/genelist.txt


#############TOPGO and DAVIDgPROFILER ######Rscript##############

(cd "$OUTPUTDIR" && Rscript "$WORKDIR"/topGO.R)
cd "$OUTPUTDIR"
(cd "$OUTPUTDIR" && Rscript "$WORKDIR"/DAVIDgProfiler.R)


##########Complete report.tsv#############################

# Match for every ENSG its ontologies in topGO.txt
for i in $(cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean0.txt ); do printf $i"\n"; grep $i -w "$OUTPUTDIR"/topGO/topGO.txt; done | cut -f2 | awk '/^ENSG/ {printf "\n"$0"\t"; next} {printf $0";"}' | sed 's/INTERGENIC;/\nINTERGENIC/g' | sed '1 s/Gene;/Gene\tGO/g' | sed 's/GO:[0-9]*://g' | sed 's/;/; /g' > "$OUTPUTDIR"/VEP/topaste0.txt

# Match for every SYMBOL (there are not ENSG available in DAVID) its KEGG pathways in DAVID_formatted.txt
for i in `cut -f18 "$OUTPUTDIR"/VEP/condeloutputclean0.txt`; do printf "CADENADECONTROL\n"; grep $i -w "$OUTPUTDIR"/DAVID/DAVID_formatted.txt | grep 'KEGG'; done | cut -f4 | tr '\n' ' ' | sed 's/CADENADECONTROL/\n/g' | sed 's/^ KEGG/KEGG/' | sed 's/ KEGG/; KEGG/g'  | sed '1d' | sed '1d'> "$OUTPUTDIR"/VEP/topaste1.txt


# Match for every ENSG its genenames with the gene names from topGO result file ensg_genenames.txt
for i in $(cut -f4 "$OUTPUTDIR"/VEP/condeloutputclean0.txt ); do printf $i"\t\n"; grep $i -w -m1 "$OUTPUTDIR"/topGO/ensg_genenames.txt; done | awk 'BEGIN{FS="\t"; OFS="\t"} /^INTERGENIC/ {printf $0"this region is not a gene\n"; next}{print $0} ' | cut -f2 | sed '/^$/d' > "$OUTPUTDIR"/topGO/topaste2.txt

#write "KEGG" header in topaste1.txt and write header "description" in topaste2.txt
sed -i '1i KEGG' "$OUTPUTDIR"/VEP/topaste1.txt
sed -i '1i Gene_Description' "$OUTPUTDIR"/topGO/topaste2.txt


# Paste KEGG,  the complete Gene Name , GO, beside condeloutputclean0
if [ -f "$OUTPUTDIR"/VEP/topaste1.txt ];
then
paste "$OUTPUTDIR"/VEP/condeloutputclean.txt "$OUTPUTDIR"/VEP/topaste1.txt > "$OUTPUTDIR"/VEP/condeloutputclean1.txt
fi

paste "$OUTPUTDIR"/VEP/condeloutputclean1.txt "$OUTPUTDIR"/topGO/topaste2.txt > "$OUTPUTDIR"/VEP/condeloutputclean2.txt
paste "$OUTPUTDIR"/VEP/condeloutputclean2.txt "$OUTPUTDIR"/VEP/topaste0.txt > "$OUTPUTDIR"/CompleteReport.tsv


# CompleteReport.tsv Summary
# offline
#cut -f1-4,7,14,17,18,21,22,23,24,26 "$OUTPUTDIR"/CompleteReport.tsv | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > "$OUTPUTDIR"/SummaryReport0.tsv

#online (for nearest gene)
cut -f1-4,7,14,17,18,21,22,23,24,25,26,27,28,29,31 "$OUTPUTDIR"/CompleteReport.tsv | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$3,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' > "$OUTPUTDIR"/SummaryReport0.tsv

BACKUPIFS=$IFS
IFS=$'\n'
for vargene in `cut -f3,4,5 SummaryReport0.tsv | sort -u | sort -n `; do grep -m1 -P "$vargene\t" SummaryReport0.tsv; done > "$OUTPUTDIR"/SummaryReport.tsv
IFS=$BACKUPIFS

# remove extra intermediary files


rm "$OUTPUTDIR"/VEP/condeloutput.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean0.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean1.txt
rm "$OUTPUTDIR"/VEP/condeloutputclean2.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene1.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene2.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene2_tmp.txt
rm "$OUTPUTDIR"/VEP/variants_per_gene3.txt
rm "$OUTPUTDIR"/VEP/topaste0.txt
rm "$OUTPUTDIR"/VEP/topaste1.txt
rm "$OUTPUTDIR"/VEP/condeloutput.vcf
rm "$OUTPUTDIR"/topGO/topaste2.txt
rm "$OUTPUTDIR"/topGO/ensg_genenames.txt
rm "$OUTPUTDIR"/SummaryReport0.tsv
rm "$OUTPUTDIR"/genes_annotations2.tsv
#rm "$OUTPUTDIR"/VEP/genelist.txt
#rm "$OUTPUTDIR"/VEP/genelistENS.txt

exit
