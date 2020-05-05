bcftools norm -m-any --check-ref w -f /mnt/isilon/dbhi_bfx/reference/human/g1k_v37/human_g1k_v37.fasta /scr1/users/shraimr/veo.remicade.062419.reheader.vcf | bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | bcftools view -t ^X,Y | bcftools norm -Ob --rm-dup both > /scr1/users/shraimr/veo.remicade.062419.bcf
bcftools index /scr1/users/shraimr/veo.remicade.062419.bcf
