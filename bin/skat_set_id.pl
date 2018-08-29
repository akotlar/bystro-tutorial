#! /usr/bin/env perl
use 5.10.0;
use strict;
use warnings;

use List::Util qw/max/;

use DDP;
use Math::Complex;
use PDL::Stats;
use PDL::Stats::Distr;
use List::MoreUtils qw/first_index uniq/;
use Getopt::Long;
use Yaml::XS qw/LoadFile/;
use Scalar::Util qw/looks_like_number/;

my ($ann, $out, $pathwayYaml);

GetOptions (
  'in|i=s' => \$ann,
  'config|c' => \$pathwayYaml,
  'out|o=s' => \$out,
);

if(!$ann) {
  die "--in|i|ann|a --out|o (optional)"
}

my ($setColName, $popName);

if($pathwayYaml) {
  my $config = LoadFile($pathwayYaml);

  $setColName = $config->{gene};
  $popName = $config->{population};

  if(!$popName || $popName eq 'all') {
    $popName = '';
  } else {
    $popName = '_' . $popName;
  }
}

$popName //= '';
$setColName //= 'refSeq.nearest.name2';

my $annFh;
my $outName = $out;
if($ann =~ /\.gz$/) {
  $outName //= substr($ann, 0, length($ann) - 3);

  # Don't think this actually would die, since separate process forked via pipe
  open($annFh, '-|', "gunzip -c $ann") or die "Couldn't open annotated snp file.";
} else {
  $outName //= $ann;
  open($annFh, '<', $ann) or die "Couldn't open annotated snp file.";
}

open(my $outFh, '>', ($outName) . ".setid_genes.tsv") or die "Couldn't open set id out file";

open(my $weightGnomadBetaFh, '>', ($outName) . ".weights_gnomad_beta.tsv") or die "Couldn't open weight out file";
open(my $weightBetaMadsenFh, '>', $outName . ".weights_gnomad_beta_madsen.tsv") or die "Couldn't open weight out file";

my $typeColIdx = 2;

my $firstLine = <$annFh>;
chomp $firstLine;

my $gnomadGenomesAf = 'gnomad.genomes.af' . $popName;
my $gnomadExomesAf = 'gnomad.exomes.af' . $popName;
my $gnomadGenomesAn = 'gnomad.genomes.an' . $popName;
my $gnomadExomsAn = 'gnomad.exomes.an' . $popName;

my @headers = split("\t", $firstLine);

p @headers; 

my %sets = (
  genes => {},
);

my %setFhs = (
  genes => $outFh,
);

my %weights = (
  beta => {},
  madsen => {},
);

my %weightFhs = (
  beta => $weightGnomadBetaFh,
  madsen => $weightBetaMadsenFh,
);

my $setColIdx = first_index{ $_ eq $setColName } @headers;
my $gnomadGenomesIdIdx = first_index { $_ eq 'gnomad.genomes.id' } @headers;
my $gnomadExomesIdIdx = first_index { $_ eq 'gnomad.exomes.id' } @headers;
my $dbSnpNameIdx = first_index { $_ eq 'dbSNP.name' } @headers;

my $gnomadGenomesAfIdx = first_index { $_ eq $gnomadGenomesAf } @headers;
my $gnomadExomesAfIdx = first_index { $_ eq $gnomadExomesAf } @headers;
my $gnomadGenomesNfeAnIdx = first_index { $_ eq $gnomadGenomesAn } @headers;
my $gnomadExomesNfeAnIdx = first_index { $_ eq $gnomadExomsAn } @headers;

my $sampleMafIdx = first_index { $_ eq 'sampleMaf' } @headers;

my $altIdx = first_index { $_ eq 'alt' } @headers;

if($setColIdx == -1) {
  die "set column $setColName not found in $ARGV[0] header";
}

if($gnomadGenomesIdIdx == -1) {
  die "Couldn't find gnomad genomes id col";
}

if($gnomadExomesIdIdx == -1) {
  die "Couldn't find gnomad exomes id col";
}

if($dbSnpNameIdx == -1) {
  die "Couldn't find dbSNP.name col";
}

if($gnomadGenomesAfIdx == -1) {
  die "couldn't find gnomad.genomes.af_nfe $ARGV[0] header";
}

if($gnomadExomesAfIdx == -1) {
  die "couldn't find gnomad.exomes.af_nfe $ARGV[0] header";
}

if($gnomadExomesAfIdx == -1) {
  die "couldn't find gnomad.exomes.af_nfe $ARGV[0] header";
}

my %variantCounts;
my %weightTrackerBeta;
my %weightTrackerBetaMadsen;
my %unique;

while (<$annFh>) {
  chomp;
  my @fields = split '\t', $_;

  my $idVal;
  # my $gnomadInverseWeight = getGnomadInverseWeight(\@fields);
  my $gnomadBetaWeight = getGnomadBetaWeight(\@fields, 1, 25);
  # my $gnomadInterpWeight = getGnomadInterpolatedWeight(\@fields);
  # my $gnomadRawWeight = getGnomadRawWeight(\@fields);
  my $gnomadBetaMadsenWeight = getGnomadBetaWeight(\@fields, 0.5, 0.5);

  $idVal = $fields[0] . "_" . $fields[1] . "_" . $fields[$altIdx];

  my $setValRef = makeUniq($fields[$setColIdx], \@fields);

  for(my $i = 0; $i < @{$setValRef}; $i++) {
    my $setVal = $setValRef->[$i];

    if(!defined $setVal) {
      # chr_pos
      $setVal = $fields[0] . "_" . $fields[1];
    }

    $variantCounts{$setVal}++;
    my $setSnp = "$setVal\t$idVal";

    if(!$unique{$setSnp}) {
      $sets{genes}{$setVal} //= [];
      push @{$sets{genes}{$setVal}}, $idVal;
      $unique{$setSnp} = 1;
    }
  }


  if(!defined $weights{beta}{$idVal}) {
    $weights{beta}{$idVal} = $gnomadBetaWeight;
    $weights{madsen}{$idVal} = $gnomadBetaMadsenWeight;
  } else {
    say STDERR "Previously have seen $idVal\n";
    say STDERR join("\t", @fields);

    if($gnomadBetaWeight > $weights{beta}{$idVal}) {
      say STDERR "Updating $idVal beta weight from $weightTrackerBeta{$idVal} to $gnomadBetaWeight\n";
      $weights{beta}{$idVal} = $gnomadBetaWeight;
    }

    if($gnomadBetaMadsenWeight > $weights{madsen}{$idVal}) {
      say STDERR "Updating $idVal madsen weight from $weightTrackerBetaMadsen{$idVal} to $gnomadBetaMadsenWeight\n";
      $weights{madsen}{$idVal} = $gnomadBetaMadsenWeight;
    }
  }
}
say "\n";

# Remove sets with 2 or fewer variants
for my $setType (keys %sets) {
  for my $setName (keys %{$sets{$setType}}) {
    if($variantCounts{$setName} < 3) {
      say STDERR "Dropping set $setName because found $variantCounts{$setName} variants";
      delete $sets{$setType}{$setName};
    }
  }
}

for my $setType (keys %sets) {
  my $fh = $setFhs{$setType};

  for my $setName (keys %{$sets{$setType}}) {
    for my $idVal (@{$sets{$setType}{$setName}}) {
      say $fh "$setName\t$idVal";

      for my $type (keys %weights) {
        if(!defined $weights{$type}{$idVal}) {
          next;
        }

        my $fh = $weightFhs{$type};

        say $fh "$idVal\t$weights{$type}{$idVal}";
      }
    }
  }
}

sub getGnomadInterpolatedWeight {
  my $fieldsAref = shift;

  if($fieldsAref->[$gnomadExomesAfIdx] ne '!') {
    my ($err, $af) = getGnomadWeight($fieldsAref->[$gnomadExomesNfeAnIdx], $fieldsAref->[$gnomadExomesAfIdx]);

    if($err) {
      die $err;
    }

    return $af;
  }

  if($fieldsAref->[$gnomadGenomesAfIdx] ne '!' ) {
    my ($err, $af) = getGnomadWeight($fieldsAref->[$gnomadGenomesNfeAnIdx], $fieldsAref->[$gnomadGenomesAfIdx]);

    if($err) {
      die $err;
    }

    return $af;
  }

  # say STDERR "Couldn't find gnomad weight, using sampleMaf for; $fieldsAref->[0]\:$fieldsAref->[1] ; genomes: $fieldsAref->[$gnomadGenomesAfIdx] ; exomes: $fieldsAref->[$gnomadExomesAfIdx] \n";

  return $fieldsAref->[$sampleMafIdx];
}

sub getGnomadWeight {
  my ($an, $af) = @_;

  if(!(looks_like_number($an) && looks_like_number($af))) {
    return ("Not a nubmer", undef);
  }

  if($af == 0) {
    return (undef, 1/($an + 1));
  }

  return (undef, $af);
}

sub getGnomadBetaWeight {
  my ($fieldsAref, $a, $b) = @_;

  if(!defined $a) {
    $a = 1;
  }

  if(!defined $b) {
    $b = 25;
  }

  my $af = getGnomadInterpolatedWeight($fieldsAref);

  return pdf_beta($af, $a, $b);
}



# If there are MULTIALLELIC sites, or indels, take only the first allel or first
# position's value
# If there are multiple overlapping transcripts
my %uniq;
sub makeUniq {
  my ($setVal, $fieldAref) = @_;

  my @allSets;

  # NOTE:
  # DGCR14 renamed to ESS2 ; it is one of our target regions
  # Excluding sema6a-as1 because 12 variants, all overlapping SEMA6A, which we want
  # Eclusing UBXN7-AS1 for the same reasons
  # Excluding KTI12 because not wanted and all variants overlap TXNDC12
  # Excluding PIGF for same reason
  # NOT SURE WHY BUT WE HAVE SEMA6A-AS2, which doesn't fully overlap SEMA6A variants; keeping for now
  # Excluding FANCL because not wanted, completely overlaps variants of VRK2, which is wanted
  # Excluding RIMBP3C because all but 10 variants (UTR) overlap RIMBP3B, and 
  # the requested targets are always for RIMBP3B_RIMBP3C, never RIMBP3C alone
  # NOT SURE WHY BUT WE HAVE ATP6V1E2; doesn't overlap anything
  # Excluding NCBP2-AS1 because overlaps NCBP2 and not wanted
  for my $setData (split /[|]/, $setVal) {
    my @setVals = uniq sort { $a cmp $b } split /[;]/, $setData;

    push @allSets, \@setVals;
  }

  my %major;
  my $max = 0;
  my $maxGene;

  for my $sRef (@allSets) {
    for my $s (@$sRef) {
      $major{$s}++;

      if($major{$s} > $max) {
        $max = $major{$s};
        $maxGene = $s;
      }
    }
  }

  # say STDERR "Major is $maxGene";
  # p %major;
  return [$maxGene];
}

say STDERR "\n\nDone.\nSet values found: " . (scalar keys %unique);