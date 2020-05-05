#Using Plink/Seq to complete burden test

#Create Pseq project
pseq veo-ibd new-project --resources /scr1/users/shraimr/hg19
pseq veo-ibd load-vcf --vcf veo.ibd.vcf
pseq veo-ibd load-pheno --file pop.pheno
pseq veo-ibd loc-summary > summary.txt

#Complete burden test
pseq veo-ibd assoc --phenotype phe1 --mask loc.group=refseq > burdentest.txt

#burden test with a minor allele frquency set 
pseq veo-ibd assoc --phenotype phe1 --mask loc.group=refseq maf=0-0.05 > burdentest.MAF0.05.txt