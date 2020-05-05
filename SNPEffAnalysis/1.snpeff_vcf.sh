#!/bin/sh
#$ -cwd

/usr/bin/java -Xmx8g -jar /mnt/isilon/dbhi_bfx/bin/snpeff/4.3/snpEff/snpEff.jar -v hg19 in.vcf > out.vcf
