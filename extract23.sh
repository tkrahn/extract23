#!/bin/bash

# extract23 Version 1.3 

OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Fetch parameters from commandline

while getopts ":vt:o:r:b:" opt; do
    case "$opt" in
    v)  verbose=1
        ;;
    b)  BAMFILE_SORTED=$OPTARG
        ;;
    r)  REF=$OPTARG
        ;;
    t)  REF_23ANDME=$OPTARG
        ;;
    o)  OUTFILE=$OPTARG
        ;;
    *)
        echo "======================================== EXTRACT23 ======================================="
        echo "Usage: "
        echo "extract23.sh -b <sorted_hg19_bamfile.bam> -r <hg19_ref.fasta> [-t <23andMe_V3_hg19_ref.tab.gz>] [-o output.txt] [-v]"
        echo "  Parameters:"
        echo "  -b The whole genome hg19 referenced, sorted and indexed BAM file. REQUIRED!"
        echo "  -r The hg19 reference file (downloaded to your computer and indexed with samtools index <ref.fasta>. REQUIRED!"
        echo "  -t The bgziped and tabix -s1 -b2 -e2 indexed translation database. Defaults to 23andMe_V3_hg19_ref.tab.gz"
        echo "  -o Output file name without the .zip suffix. Defaults to 23andMe_V3_hg19.txt"
        echo "  -v Verbose output."
        exit 0
        ;;
    esac
done

shift $((OPTIND-1))


if [ -z "${BAMFILE_SORTED}" ]; then
    echo "BAM file required (option -b)"
    exit 0
fi

if [ -z "${REF}" ]; then
    echo "FASTA reference file required (option -r)"
    exit 0
fi

if [ -z "${REF_23ANDME}" ]; then
    echo "Using Pre-defined translation database 23andMe_V3_hg19_ref.tab.gz (option -t)"
    REF_23ANDME="23andMe_V3_hg19_ref.tab.gz"
fi

if [ -z "${OUTFILE}" ]; then
    echo "Using Pre-defined output file name 23andMe_V3_hg19.txt (option -o)"
    OUTFILE="23andMe_V3_hg19.txt"
fi

if [ -z "${verbose}" ]; then
    verbose=0
fi

# Generate 23andMe mockup file

if [ ${verbose} -gt 0 ]; then
    echo "Starting mpileup... Please be patient!"
fi

# The -l parameter requires mpileup to exactly pick the constellations at the SNP positions
# without having to screen through the whole genome base by base.
samtools mpileup -C 50 -v -l ${REF_23ANDME} -f ${REF} ${BAMFILE_SORTED} > 23andMe_raw.vcf.gz
tabix -p vcf 23andMe_raw.vcf.gz

if [ ${verbose} -gt 0 ]; then
    echo "Mpileup completed. Starting SNP calling..."
fi


# Now we call the SNPs from the raw mpileup data with the -m (mixed base) genotype caller
bcftools call -O z -V indels -m -P 0 23andMe_raw.vcf.gz > 23andMe_called.vcf.gz
tabix -p vcf 23andMe_called.vcf.gz

if [ ${verbose} -gt 0 ]; then
    echo "SNP calling completed. Starting annotation..."
fi


# Here we annotate the SNP names (rs numbers) to each SNP position
bcftools annotate -O z -a ${REF_23ANDME} -c CHROM,POS,ID 23andMe_called.vcf.gz > 23andMe_annotated.vcf.gz
tabix -p vcf 23andMe_annotated.vcf.gz

if [ ${verbose} -gt 0 ]; then
    echo "Annotation completed. Starting extraction from VCF ..."
fi


# Pick the data from the vcf and convert it into a tab delimited table.
# The loop brackets [] are required for the translated genotypes, even though we only have one sample.
# We have to reformat several things with the sed editor in a pipe to make the format compatible.
bcftools query -f '%ID\t%CHROM\t%POS[\t%TGT]\n' 23andMe_annotated.vcf.gz | \
    sed 's/chr//' | \
    sed 's/\tM\t/\tMT\t/g' | \
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
    sed 's/CA$/AC/' > 23andMe_V3_hg19.tab

if [ ${verbose} -gt 0 ]; then
    echo "Extraction from VCF completed. Sorting by chromosome and position ..."
fi

# Finally sort the table by chromosome and position
sort -t $'\t' -k2,3 -V 23andMe_V3_hg19.tab > 23andMe_V3_hg19_sorted.tab


# 23andMe header
# I have implemented exactly the original 23andMe header here to avoid compatibility issues.
# However I have no clue about copyright issues, so use at your own risk!
# Probably most tools will ignore the # comment lines, so I recommend to change the header and only use this original 
# 23andMe header when compatibility issues occur.

echo '# This data file generated by 23andMe at: Thu Dec 29 11:59:59 2016' > ${OUTFILE}
echo '#' >> ${OUTFILE}
echo '# This file contains raw genotype data, including data that is not used in 23andMe reports.' >> ${OUTFILE}
echo '# This data has undergone a general quality review however only a subset of markers have been ' >> ${OUTFILE}
echo '# individually validated for accuracy. As such, this data is suitable only for research, ' >> ${OUTFILE}
echo '# educational, and informational use and not for medical or other use.' >> ${OUTFILE}
echo '# ' >> ${OUTFILE}
echo '# Below is a text version of your data.  Fields are TAB-separated' >> ${OUTFILE}
echo '# Each line corresponds to a single SNP.  For each SNP, we provide its identifier ' >> ${OUTFILE}
echo '# (an rsid or an internal id), its location on the reference human genome, and the ' >> ${OUTFILE}
echo '# genotype call oriented with respect to the plus strand on the human reference sequence.' >> ${OUTFILE}
echo '# We are using reference human assembly build 37 (also known as Annotation Release 104).' >> ${OUTFILE}
echo '# Note that it is possible that data downloaded at different times may be different due to ongoing ' >> ${OUTFILE}
echo '# improvements in our ability to call genotypes. More information about these changes can be found at:' >> ${OUTFILE}
echo '# https://www.23andme.com/you/download/revisions/' >> ${OUTFILE}
echo '# ' >> ${OUTFILE}
echo '# More information on reference human assembly build 37 (aka Annotation Release 104):' >> ${OUTFILE}
echo '# http://www.ncbi.nlm.nih.gov/mapview/map_search.cgi?taxid=9606' >> ${OUTFILE}
echo '#' >> ${OUTFILE}
echo '# rsid	chromosome	position	genotype' >> ${OUTFILE}

# Append the genotype table to the header
cat 23andMe_V3_hg19_sorted.tab >> ${OUTFILE}

if [ ${verbose} -gt 0 ]; then
    echo "${OUTFILE} was created. Compressing ..."
fi

# Zip the file
zip ${OUTFILE}.zip ${OUTFILE}
echo "extract23: Output file ${OUTFILE}.zip was created."

# Delete all no longer needed files
rm -f 23andMe_raw.vcf.*
rm -f 23andMe_called.vcf.*
rm -f 23andMe_annotated.vcf.*
rm -f 23andMe_V3_hg19*.tab
rm -f ${OUTFILE}




