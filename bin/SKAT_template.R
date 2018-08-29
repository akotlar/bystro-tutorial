library("SKAT")

# This is a template; needs to be modified before use. Everything named "_placeholder" needs to sed replaced with a normal value
args = commandArgs(trailingOnly=TRUE)

if(length(args) == 0) {
  args = c(FALSE,FALSE,FALSE,FALSE);
}

if(args[1]) {
  inName = args[1];
} else {
  inName = 'inName_placeholder';
}

tCwd = 'cwd_placeholder'
if(args[2]) {
  tCwd = args[2];
}

if(tCwd != 'cwd_placeholder') {
  setwd(tCwd)
}

rCorr = 'rCorr_placeholder'
if(args[3]) {
  rCorr = as.numeric(args[3]);
}

file.cov = 'cov_placeholder'
if(args[4]) {
  file.cov = args[4];
}

if(file.cov == 'cov_placeholder') {
  file.cov = NULL;
}

# Just for backward compat
if(!is.numeric(rCorr)) {
  rCorr = 0
}

baseName = paste(inName, "skat",collapse=".",sep=".")

# Using only rare variants below
file.bed <- paste(inName, 'bed', collapse=".",sep=".")
file.bim <- paste(inName, 'bim', collapse=".",sep=".")
file.fam <- paste(inName, 'fam', collapse=".",sep=".")

file.geneSet <- paste(baseName, 'setid_genes.tsv', collapse=".",sep=".")

file.gnomadBetaWeights = Read_SNP_WeightFile(FileName=paste(baseName, 'weights_gnomad_beta.tsv',  collapse=".", sep="."))
file.gnomadBetaMadsenWeights = Read_SNP_WeightFile(FileName=paste(baseName, 'weights_gnomad_beta_madsen.tsv',  collapse=".", sep="."))

file.ssdGenes <- paste(baseName, 'genes.ssd',  collapse=".", sep=".") # This file is about to be created with the below command
file.infoGenes <- paste(baseName, 'genes.info',  collapse=".", sep=".")

Generate_SSD_SetID(file.bed, file.bim, file.fam, file.geneSet, file.ssdGenes, file.infoGenes)
ssd.infoGenes <- Open_SSD(file.ssdGenes, file.infoGenes)

## DEFINE YOUR OWN COVARIATES
if(!is.null(file.cov)) {
  famCov <- Read_Plink_FAM_Cov(file.fam, file.cov, Is.binary=TRUE, cov_header=FALSE)
  
  # Define your own
  x1 <- famCov$COV1
  x2 <- famCov$COV2
  x3 <- famCov$COV3
  x4 <- famCov$COV4
  
  yCov <- famCov$Phenotype
  sex <- famCov$Sex

  obj<-SKAT_Null_Model(yCov ~ x1 + x2 + x3 + x4 + sex, out_type="D") #n.Resampling=1000, type.Resampling='bootstrap'
} else {
  famCov <- Read_Plink_FAM(file.fam, Is.binary=TRUE)
 
  yCov <- famCov$Phenotype
  sex <- famCov$Sex

  obj<-SKAT_Null_Model(yCov ~ sex, out_type="D") #n.Resampling=1000, type.Resampling='bootstrap'
}


warnings()

#### Run basic model with gnomadBeta weights
writeLines(paste("\nTesting",baseName,"gene set with gnomadBetaWeights",collapse=" "))

out <- SKATBinary.SSD.All(ssd.infoGenes, obj, kernel = "linear.weighted", method="optimal.adj", r.corr=rCorr, obj.SNPWeight=file.gnomadBetaWeights)
out.results <- out$results
toWrite = out.results[order(out.results$P.value), c(1, 2, 4)]
write.table(toWrite, file=paste(baseName, 'genes_gnomadBetaWeights.tsv',collapse=".", sep="."), col.names=T, row.names=F, quote=F, sep='\t')
warnings()

min(out.results$P.value)
signif = Resampling_FWER(out,FWER=0.05)
signif
write.table(signif$result, file=paste(baseName, 'genes_gnomadBetaWeights_fwer.tsv',collapse=".", sep="."), col.names=T, row.names=F, quote=F, sep='\t')


#### Run basic model with no weights (should be used only for our < 1e-4 model)
writeLines(paste("\nTesting",baseName,"gene set unweighted",collapse=" "))

out <- SKATBinary.SSD.All(ssd.infoGenes, obj, kernel = "linear.weighted", method="optimal.adj", r.corr=rCorr, weights.beta=c(1,1))
out.results <- out$results
toWrite = out.results[order(out.results$P.value), c(1, 2, 4)]
write.table(toWrite, file=paste(baseName, 'genes_unweighted.tsv',collapse=".", sep="."), col.names=T, row.names=F, quote=F, sep='\t')
warnings()

min(out.results$P.value)
signif = Resampling_FWER(out,FWER=0.05)
signif
write.table(signif$result, file=paste(baseName, 'genes_unweighted_fwer.tsv',collapse=".", sep="."), col.names=T, row.names=F, quote=F, sep='\t')
