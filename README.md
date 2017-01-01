# extract23
Extract a simulated 23andMe (V3) style file from a Whole Genome BAM file

Many DNA and genealogy tools use a file upload for allele call (summary) files in the 23andMe format (http://fileformats.archiveteam.org/wiki/23andMe). Those tools can't cope with huge VCF files or even raw FastQ or BAM files. With this script I show how it is possible to extract exactly the SNPs listed in a template from a BAM file and reformat the data to a 23andMe style allele table that can be uploaded to various interpretation services such as Gedmatch or Promethease. Of course you can change the template to any other arbitrary format such as Ancestry or FTDNA. I have just started with 23andMe V3 since this is the most universal usable format.


# Requirements:

The extract23.sh script depends on bash, a modern UNIX/Linux shell and is only tested on Debian Linux. I will not work on Windows platforms without significant rewriting, but you may try to use it with Cygwin https://www.cygwin.com/.

The script makes heavy use of the htslib versions of samtools and bcftools and tabix http://www.htslib.org/ . At the time of writing I was using samtools Version: 1.3.1-42-g0a15035 and bcftools 1.3.1-201-g87456cf (using htslib 1.3.2-176-gafd9b56). Older versions before inclusion of htslib will likely not work because features like bcftools call were not available at that time. HTSlib itself depends on zlib. Since the 23andMe files usually come in a Zip format, you need to have zip installed to compress it.

SNP calling of huge BAM files requires more than average computing power, disk space and RAM. The minimum recommended setup is 4 Gbyte of RAM and 100 Gbyte free disk space. A 64 bit processor is recommended.

With the given template the whole genome BAM file must be in the hg19 format. This means that the reference FASTA files must have chr1, chr2, ... chr22, chrM, chrX and chrY. If the BAM file is in GRCh37 format (1, 2, 3 ... 22, M, X, Y) you may try to change the template by removing the leading 'chr' characters like this:

<code>
gunzip 23andMe_V3_hg19_ref.vcf.gz > 23andMe_V3_hg19_ref.vcf
cat 23andMe_V3_hg19_ref.vcf | sed 's/^chr//' > 23andMe_V3_GRCh37_ref.vcf
bgzip 23andMe_V3_GRCh37_ref.vcf > 23andMe_V3_GRCh37_ref.vcf.gz
tabix -s1 -b2 -e2 23andMe_V3_GRCh37_ref.vcf.gz
</code>

You can change

Of course you'll need to change the reference to the correct template in the script file itself.


# Installation/Usage:

Make sure htslib, samtools, bcftools, tabix, git, gzip and zip are installed and in your executable path.

<code>
git clone https://github.com/tkrahn/extract23
cd extract23
./extract23.sh /path/to/bamfile_in_hg19.bam 23andMe_V3_hg19_ref.vcf.gz
</code>

The BAM file and the BAI index must be in the same directory.

