# extract23

Extract a simulated 23andMe (V3) style file from a Whole Genome BAM file

Many DNA and genealogy tools use a file upload for allele call (summary) files in the 
23andMe format (http://fileformats.archiveteam.org/wiki/23andMe). Those tools can't cope 
with huge VCF files or even raw FastQ or BAM files. With this script I show how it is 
possible to extract exactly the SNPs listed in a template from a BAM file and reformat 
the data to a 23andMe style allele table that can be uploaded to various interpretation 
services such as Gedmatch or Promethease. Of course you can change the template to any 
other arbitrary format such as Ancestry or FTDNA. I have just started with 23andMe V3 
since this is the most universal usable format.


## Requirements:

The extract23.sh script depends on bash, a modern UNIX/Linux shell and is only tested 
on Debian Linux. I will not work on Windows platforms without significant rewriting, 
but you may try to use it with Cygwin https://www.cygwin.com/.

The script makes heavy use of the htslib versions of samtools and bcftools and tabix 
http://www.htslib.org/ . At the time of writing I was using 
samtools Version: 1.3.1-42-g0a15035 and bcftools 1.3.1-201-g87456cf 
(using htslib 1.3.2-176-gafd9b56). Older versions before inclusion of htslib will 
likely not work because features like bcftools call were not available at that time. 
HTSlib itself depends on zlib. Since the 23andMe files usually come in a Zip format, 
you need to have zip installed to compress it.

SNP calling of huge BAM files requires more than average computing power, disk space 
and RAM. The minimum recommended setup is 4 Gbyte of RAM and 100 Gbyte free disk 
space. A 64 bit processor is recommended.

Of course you need a hg19 referenced, sorted and indexed WGS BAM file. I haven't tried
an exome BAM file, but it might work, since the 23andMe SNPs are often in exome 
regions. This is certainly sufficient for single SNP diagnostics, but it may be 
problematic for segment analysis or phasing.

With the given template the whole genome BAM file must be in the hg19 format. This 
means that the reference FASTA files must have chr1, chr2, ... chr22, chrM, chrX and 
chrY. If the BAM file is in GRCh37 format (1, 2, 3 ... 22, M, X, Y) you may try to 
change the template by removing the leading 'chr' characters like this:

```bash
gunzip -c 23andMe_V3_hg19_ref.tab.gz > 23andMe_V3_hg19_ref.tab
cat 23andMe_V3_hg19_ref.tab | sed 's/^chr//' > 23andMe_V3_GRCh37_ref.tab
bgzip -c 23andMe_V3_GRCh37_ref.tab > 23andMe_V3_GRCh37_ref.tab.gz
tabix -s1 -b2 -e2 23andMe_V3_GRCh37_ref.tab.gz
```

Of course you'll need to change the correct template in the script 
call itself with the help of the -t parameter.


## Installation/Usage:

Make sure htslib, samtools, bcftools, tabix, git, gzip and zip are installed and 
in your executable path. You need to have the samtools indexed hg19 reference on
your computer and point to it with the -r parameter. The .fai file must be in the same
directory as the hg19_ref.fasta file.

```bash
git clone https://github.com/tkrahn/extract23
cd extract23
./extract23.sh -b /path/to/bamfile_in_hg19.bam -r /path/to/hg19_ref.fasta -v
```

The BAM file and the BAI index must be in the same directory.
Be patient, the run for a 30x WGS BAM file can take 15 minutes up to hours, depending 
on the speed of your processor.

```
======================================== EXRACT23 =======================================
Usage: 
extract23.sh -b <sorted_hg19_bamfile.bam> -r <hg19_ref.fasta> [-t <23andMe_V3_hg19_ref.tab.gz>] [-o output.txt] [-v]
  Parameters:
  -b The whole genome hg19 referenced, sorted and indexed BAM file. REQUIRED!
  -r The hg19 reference file (downloaded to your computer and indexed with samtools index <ref.fasta>. REQUIRED!
  -t The bgziped and tabix -s1 -b2 -e2 indexed translation database. Defaults to 23andMe_V3_hg19_ref.tab.gz
  -o Output file name without the .zip suffix. Defaults to 23andMe_V3_hg19.txt
  -v Verbose output.
```



## Issues

The formatting is very simplistic and currently completely ignores Indels. The script 
actually calls Indels like point mutations with the single base at the given location in 
the template. This should be improved with another lookup table to properly call
Indels and translate them to 23andMe's DD, DI, II format. However neither Gedmatch or 
Promethease care about incorrect Indels, so I didn't look into this further.

Another difference is that haplotype SNPs on the X, Y and MT chromosomes are still
called with two alleles. Not sure if there is any tool available that has problems with 
this. 

