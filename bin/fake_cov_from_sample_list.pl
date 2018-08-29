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

my ($sampleList, $out, $pathwayYaml);

GetOptions (
  'in|i=s' => \$sampleList,
  'out|o=s' => \$out,
);

if(!$sampleList) {
  die "--in|i sampleList --out|o (optional)"
}

my $annFh;
my $outName = $out;
if($sampleList =~ /\.gz$/) {
  $outName //= substr($sampleList, 0, length($sampleList) - 3) . ".fake.cov";

  # Don't think this actually would die, since separate process forked via pipe
  open($annFh, '-|', "gunzip -c $sampleList") or die "Couldn't open sampleList file.";
} else {
  $outName //= $sampleList . ".fake.cov";
  open($annFh, '<', $sampleList) or die "Couldn't open annotated snp file.";
}

open(my $outFh, '>', $outName) or die "Couldn't make fam file";

my @lines = <$annFh>;
chomp @lines;

my @range = ( -10000 .. 10000 );

my $i = -1;
for my $id (@lines) {
    $i++;

    my $cov1 = $range[rand(@range)] / 10_000;
    my $cov2 = $range[rand(@range)] / 10_000;
    my $cov3 = $range[rand(@range)] / 10_000;
    my $cov4 = $range[rand(@range)] / 10_000;

    my $fam = $id;

    my $dad = -9;
    my $mom = -9;

    say $outFh "$fam\t$id\t$cov1\t$cov2\t$cov3\t$cov4";
}