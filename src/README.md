# Instruction

1. Download the relevant assembly .sdx and .seq file, place in the same directory
2. Run (assumes pigz installed):
   ./bin/bystro_to_vcf assembly.sdx <(pigz -d -c annotation.tsv.gz) annotation.sample_list | pigz -c > 1000g_sample_annotation.vcf.gz

[hg19.sdx](https://s3.amazonaws.com/bystro-source/hg19.sdx)
[hg19.seq](https://s3.amazonaws.com/bystro-source/hg19.seq)
[hg38.sdx](https://s3.amazonaws.com/bystro-source/hg38.sdx)
[hg38.seq](https://s3.amazonaws.com/bystro-source/hg38.seq)
