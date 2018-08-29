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
  'out|i=s' => \$out,
);

if(!$sampleList) {
  die "--in|i|ann|a --out|o (optional)"
}

my $annFh;
my $outName = $out;
if($sampleList =~ /\.gz$/) {
  $outName //= substr($sampleList, 0, length($sampleList) - 3) . ".fake.sample_list";

  # Don't think this actually would die, since separate process forked via pipe
  open($annFh, '-|', "gunzip -c $sampleList") or die "Couldn't open sampleList file.";
} else {
  $outName //= $sampleList . ".fake.sample_list";
  open($annFh, '<', $sampleList) or die "Couldn't open annotated snp file.";
}

open(my $outFh, '>', $outName) or die "Couldn't make sample_list file";

my @lines = <$annFh>;
chomp @lines;

my $i = -1;
for my $id (@lines) {
    $i++;

    say $outFh "$id\t$id";
}