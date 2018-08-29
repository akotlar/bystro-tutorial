#!/bin/bash

exFolder=./example

outFolder=./out-example

rm -rf $outFolder

mkdir -p $outFolder && cd $_;
mkdir -p ./log;

# * needed instead of ./ to prevent ./filename.annotaiton.tsv.gz
annotation=1000g_100klines_example_old.annotation.tsv.gz;
annotationBaseName=1000g_100klines_example.annotation;
annotationVcfGz=1000g_100klines_old_example.annotation.20klines.vcf.gz;
annotationSampleList=1000g_100klines_example.annotation.sample_list;

fam=1000g_100klines_example.fake.fam
cov=1000g_100klines_example.fake.cov

cp ../example/$annotation ./
cp ../example/$annotationVcfGz ./
cp ../example/$annotationSampleList ./

cp ../example/$fam ./ 
cp ../example/$cov ./

annotationSampleListWithFam=$annotationSampleList.with_fam
perl ../bin/add_fake_fam_to_sample_list.pl --in $annotationSampleList --out $annotationSampleListWithFam

annotationPlink="$annotationBaseName.plink";
annotationPlinkHwe="$annotationPlink.hwe_1e-6";
annotationPlinkHweSkat="$annotationPlinkHwe.skat";

skatFile="skat_$annotationPlinkHwe.R";

plink --vcf $annotationVcfGz --keep-allele-order --const-fid seq --make-bed --out $annotationPlink &> "./log/making_$annotationPlink.log";
mv "$annotationPlink.fam" "$annotationPlink.fam_nosex_affectation";

mv $fam "$annotationPlink.fam"
plink --bfile $annotationPlink --hwe 1e-6 midp --keep $annotationSampleListWithFam --make-bed --out $annotationPlinkHwe &> "./log/making_$annotationPlinkHwe.log";

perl ../bin/skat_set_id.pl --in $annotation --out $annotationPlinkHweSkat &> "./log/making_weights_for_$annotationPlinkHweSkat";

cp ../bin/SKAT_template.R ./;
# escape backslashes https://stackoverflow.com/questions/27787536/how-to-pass-a-variable-containing-slashes-to-sed
mv SKAT_template.R $skatFile;
# perl -pi -w -e "s/cwd_placeholder/.///\//\\/}/" $skatFile;
perl -pi -w -e "s/inName_placeholder/$annotationPlinkHwe/" $skatFile;
perl -pi -w -e "s/rCorr_placeholder/$rCorr/" $skatFile;
perl -pi -w -e "s/cov_placeholder/$cov/" $skatFile;

Rscript $skatFile;
cd ../;

