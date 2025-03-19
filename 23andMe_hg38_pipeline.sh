#!/bin/bash
START=$(date +%s.%N)

# setup parameters

YSEQID=${PWD##*/}
# YSEQID="1234" # (the above command simply gets the name of the last segment of the current working directory)

NUM_THREADS=$(getconf _NPROCESSORS_ONLN)
echo "We can use ${NUM_THREADS} threads."

REF_HG19="/genomes/0/refseq/hg19/hg19.fa"
REF_HG38="/genomes/0/refseq/hg38/hg38.fa"

BAMFILE="${YSEQID}_bwa-mem_hg38.bam"

BAMFILE_SORTED="${YSEQID}_bwa-mem_hg38_sorted.bam"
VCF_FILE="${YSEQID}_23andMe_hg38.vcf"

REF_23ANDME="23andMe_all_hg19_ref.tab"
REF_23ANDME_HG38="23andMe_all_hg38_ref.tab"


TEMPLATE_23ANDME="/genomes/0/refseq/hg19/23andMe_all_hg19_raw.tab"
if [ -f "${TEMPLATE_23ANDME}" ]; then
    echo "${TEMPLATE_23ANDME} already stored. Using it. Make sure the file is up to date."
    cp ${TEMPLATE_23ANDME} .
else 
    echo "${TEMPLATE_23ANDME} does not exist. Downloading ..."
    # Prepare newest 23andMe Reference from their API:
    wget -O 23andMe_all_hg19_raw.tab https://api.23andme.com/1/genome_snp_map/
fi
echo "23andMe SNP definitions available" 

echo "#CHROM	POS	ID" >23andMe_all_hg19_ref.tab

while IFS=$'\t' read -r index snp chromosome chromosome_position
do
	if [[ $index == \#* || $index =~ ^index ]]; then
		echo "skipping $index"
	else 
		CHROM=${chromosome//MT/M}
		echo "chr${CHROM}	${chromosome_position}	${snp}	--" >>23andMe_all_hg19_unsorted.tab
	fi
done < 23andMe_all_hg19_raw.tab


sort -t $'\t' -k1,2 -V 23andMe_all_hg19_unsorted.tab > 23andMe_all_hg19_sorted.tab
cat 23andMe_all_hg19_sorted.tab >> 23andMe_all_hg19_ref.tab

# Convert the 23andMe TAB file into a VCF
bcftools convert -c CHROM,POS,ID,AA -s SampleName -f ${REF_HG19} --tsv2vcf 23andMe_all_hg19_ref.tab -Ov -o 23andMe_all_hg19_ref.vcf


bgzip -c 23andMe_all_hg19_ref.tab > 23andMe_all_hg19_ref.tab.gz
tabix -s1 -b2 -e2 23andMe_all_hg19_ref.tab.gz


# CrossMap the 23andMe VCF to hg38
map_chain="/genomes/0/refseq/chain/hg19ToHg38.over.chain"
CrossMap.py vcf $map_chain 23andMe_all_hg19_ref.vcf ${REF_HG38} 23andMe_all_hg38_ref.vcf

# Sort again (due to inversions during hg19 > hg38 conversion)
bcftools sort -Oz 23andMe_all_hg38_ref.vcf -o 23andMe_all_hg38_ref_sorted.vcf.gz
tabix -p vcf 23andMe_all_hg38_ref_sorted.vcf.gz

# Convert back to TSV
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_all_hg38_ref.tab.gz
tabix -s1 -b2 -e2 23andMe_all_hg38_ref.tab.gz
echo "23andMe SNP definitions translated to hg38"

# Split in separate chr files for speedup of mpileup
bcftools query -r chr1  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr1_hg38_ref.tab.gz &
bcftools query -r chr2  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr2_hg38_ref.tab.gz &
bcftools query -r chr3  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr3_hg38_ref.tab.gz &
bcftools query -r chr4  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr4_hg38_ref.tab.gz &
bcftools query -r chr5  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr5_hg38_ref.tab.gz &
bcftools query -r chr6  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr6_hg38_ref.tab.gz &
bcftools query -r chr7  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr7_hg38_ref.tab.gz &
bcftools query -r chr8  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr8_hg38_ref.tab.gz &
bcftools query -r chr9  -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr9_hg38_ref.tab.gz &
bcftools query -r chr10 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr10_hg38_ref.tab.gz &
bcftools query -r chr11 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr11_hg38_ref.tab.gz &
bcftools query -r chr12 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr12_hg38_ref.tab.gz &
bcftools query -r chr13 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr13_hg38_ref.tab.gz &
bcftools query -r chr14 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr14_hg38_ref.tab.gz &
bcftools query -r chr15 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr15_hg38_ref.tab.gz &
bcftools query -r chr16 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr16_hg38_ref.tab.gz &
bcftools query -r chr17 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr17_hg38_ref.tab.gz &
bcftools query -r chr18 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr18_hg38_ref.tab.gz &
bcftools query -r chr19 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr19_hg38_ref.tab.gz &
bcftools query -r chr20 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr20_hg38_ref.tab.gz &
bcftools query -r chr21 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr21_hg38_ref.tab.gz &
bcftools query -r chr22 -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chr22_hg38_ref.tab.gz &
bcftools query -r chrX -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chrX_hg38_ref.tab.gz &
bcftools query -r chrY -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chrY_hg38_ref.tab.gz &
#bcftools query -r chrM -f'%CHROM\t%POS\t%REF,%ALT\n' 23andMe_all_hg38_ref_sorted.vcf.gz | bgzip -c > 23andMe_chrM_hg38_ref.tab.gz &
wait

tabix -s1 -b2 -e2 23andMe_chr1_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr2_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr3_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr4_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr5_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr6_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr7_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr8_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr9_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr10_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr11_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr12_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr13_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr14_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr15_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr16_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr17_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr18_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr19_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr20_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr21_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chr22_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chrX_hg38_ref.tab.gz &
tabix -s1 -b2 -e2 23andMe_chrY_hg38_ref.tab.gz &
#tabix -s1 -b2 -e2 23andMe_chrM_hg38_ref.tab.gz &
wait

# Generate 23andMe mockup file

# Parallel SNP calling by chromosome
echo "MpileUp started"

PARAMC=0

bcftools mpileup -R 23andMe_chr1_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr1_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr2_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr2_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr3_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr3_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr4_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr4_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr5_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr5_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr6_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr6_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr7_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr7_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr8_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr8_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr9_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr9_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr10_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr10_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr11_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr11_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr12_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr12_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr13_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr13_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr14_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr14_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr15_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr15_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr16_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr16_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr17_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr17_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr18_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr18_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr19_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr19_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr20_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr20_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr21_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr21_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chr22_hg38_ref.tab.gz -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chr22_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chrX_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrX_${VCF_FILE}.gz &
bcftools mpileup -R 23andMe_chrY_hg38_ref.tab.gz  -Ou -C ${PARAMC} -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrY_${VCF_FILE}.gz &
bcftools mpileup -r chrM                          -Ou -C 50 -f $REF_HG38 $BAMFILE_SORTED | bcftools call -O z --threads $NUM_THREADS -V indels -m -P 0 > chrM_${VCF_FILE}.gz &
wait

# Concatenate all chromosome VCFs to one big file
bcftools concat -O z chr[1-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chr[M,X-Y]_${VCF_FILE}.gz > ${YSEQID}_23andMe_called_hg38_vcf.gz
tabix -p vcf ${YSEQID}_23andMe_called_hg38_vcf.gz

# Delete no longer needed VCFs
rm -f chr[0-9]_${VCF_FILE}.gz chr[1-2][0-9]_${VCF_FILE}.gz chrX_${VCF_FILE}.gz

echo "MPileUp complete"

# Reverse translate hg38 > hg19
back_chain="/genomes/0/refseq/chain/hg38ToHg19.over.chain"
CrossMap.py vcf $back_chain ${YSEQID}_23andMe_called_hg38_vcf.gz ${REF_HG19} ${YSEQID}_23andMe_translated_hg19.vcf

echo "Translating SNP calls to hg19 complete"

# Sort again (due to inversions during hg38 > hg19 conversion) & 
bcftools sort -Oz ${YSEQID}_23andMe_translated_hg19.vcf -o ${YSEQID}_23andMe_sorted_hg19.vcf.gz
tabix -p vcf ${YSEQID}_23andMe_sorted_hg19.vcf.gz

echo "Sorting hg19 SNPs complete"

# Filter (-R) for tested chip positions only
bcftools view -Oz -R 23andMe_all_hg19_ref.tab.gz ${YSEQID}_23andMe_sorted_hg19.vcf.gz -o ${YSEQID}_23andMe_sorted_filtered_hg19.vcf.gz
tabix -p vcf ${YSEQID}_23andMe_sorted_filtered_hg19.vcf.gz

echo "Filtering hg19 SNPs complete"


# Annotate the SNP names (rs numbers)
bcftools annotate -O z -a 23andMe_all_hg19_ref.tab.gz -c CHROM,POS,ID ${YSEQID}_23andMe_sorted_filtered_hg19.vcf.gz > ${YSEQID}_23andMe_annotated_hg19.vcf.gz
tabix -p vcf ${YSEQID}_23andMe_annotated_hg19.vcf.gz

echo "SNP name annotation complete"

bcftools query -f '%ID\t%CHROM\t%POS[\t%TGT]\n' ${YSEQID}_23andMe_annotated_hg19.vcf.gz | \
    sed 's/chr//' | \
    sed 's/\tM\t/\tMT\t/' | \
    sed 's/\///' | \
    sed -E 's/^(.*)\tY\t(.*)\t(.).*$/\1\tY\t\2\t\3/' | \
    sed -E 's/^(.*)\tMT\t(.*)\t(.).*$/\1\tMT\t\2\t\3/' | \
    sed 's/\.\.$/--/' | \
    sed 's/\t$/\t--/' | \
    sed 's/TA$/AT/' | \
    sed 's/TC$/CT/' | \
    sed 's/TG$/GT/' | \
    sed 's/GA$/AG/' | \
    sed 's/GC$/CG/' | \
    sed 's/CA$/AC/' > ${YSEQID}_23andMe_all_hg19.tab

sort -t $'\t' -k2,3 -V ${YSEQID}_23andMe_all_hg19.tab > ${YSEQID}_23andMe_all_hg19_sorted.tab

# DATE=$(date +%Y-%m-%d_%H:%M:%S)  # Old legacy format
DATE=$(date +%a %b %d %H:%M:%S %Y)

# 23andMe header
echo "# This data file generated by 23andMe at: ${DATE}" > ${YSEQID}_23andMe_all_hg19.txt
echo '#' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# This file contains raw genotype data, including data that is not used in 23andMe reports.' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# This data has undergone a general quality review however only a subset of markers have been ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# individually validated for accuracy. As such, this data is suitable only for research, ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# educational, and informational use and not for medical or other use.' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# Below is a text version of your data.  Fields are TAB-separated' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# Each line corresponds to a single SNP.  For each SNP, we provide its identifier ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# (an rsid or an internal id), its location on the reference human genome, and the ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# genotype call oriented with respect to the plus strand on the human reference sequence.' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# We are using reference human assembly build 37 (also known as Annotation Release 104).' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# Note that it is possible that data downloaded at different times may be different due to ongoing ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# improvements in our ability to call genotypes. More information about these changes can be found at:' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# https://www.23andme.com/you/download/revisions/' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# ' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# More information on reference human assembly build 37 (aka Annotation Release 104):' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=9606' >> ${YSEQID}_23andMe_all_hg19.txt
echo '#' >> ${YSEQID}_23andMe_all_hg19.txt
echo '# rsid	chromosome	position	genotype' >> ${YSEQID}_23andMe_all_hg19.txt

cat ${YSEQID}_23andMe_all_hg19_sorted.tab >> ${YSEQID}_23andMe_all_hg19.txt

echo "${YSEQID}_23andMe_all_hg19.txt created"

zip ${YSEQID}_23andMe_all_hg19.zip ${YSEQID}_23andMe_all_hg19.txt

echo "${YSEQID}_23andMe_all_hg19.zip packed"

# Cleanup no longer needed files
echo "Cleaning up ..."

rm -f 23andMe_chr*ref.tab*
rm -f *${YSEQID}_23andMe_hg38.vcf.gz
rm -f ${YSEQID}_23andMe_all_hg19_sorted.tab
rm -f ${YSEQID}_23andMe_all_hg19.tab
rm -f ${YSEQID}_23andMe_all_hg19.txt
#rm -f ${YSEQID}_23andMe_annotated_hg19.vcf.gz
#rm -f ${YSEQID}_23andMe_annotated_hg19.vcf.gz.tbi
#rm -f ${YSEQID}_23andMe_called_hg38_vcf.gz
#rm -f ${YSEQID}_23andMe_called_hg38_vcf.gz.tbi
rm -f ${YSEQID}_23andMe_sorted_filtered_hg19.vcf.gz
rm -f ${YSEQID}_23andMe_sorted_filtered_hg19.vcf.gz.tbi
rm -f ${YSEQID}_23andMe_sorted_hg19.vcf.gz
rm -f ${YSEQID}_23andMe_sorted_hg19.vcf.gz.tbi
rm -f ${YSEQID}_23andMe_translated_hg19.vcf
rm -f ${YSEQID}_23andMe_translated_hg19.vcf.unmap
rm -f 23andMe_all*

END=$(date +%s.%N)

DIFF=$(echo "$END - $START" | bc)
echo "Finished in ${DIFF} Seconds"




