#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use Data::Dumper qw(Dumper);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);
use Pod::Usage;

=pod

=head1 NAME

outlierSR	aggregate outliers from isoforms by gsp from split reads in the bam

=head1 SYNOPSIS

  --gene, -g	Gene
  --job, -j	Job
  --minread, -mr	Minimum read cutoff to be an outlier
  --bgtype, -bg	Background cohort to use for cutoffs (can be "all" or a file name with cases)
  --aggdir, -ad	Directory where result files are stored
  --aggdir2, -ad2	Extra directory where result files are stored
  --aggdir3, -ad3	Extra directory where result files are stored
  --help, -h    Print this help

=head1 VERSION

    1.0

=cut

GetOptions(
	"gene|g=s" => \(my $gene = "IKZF1"),
	"case|c=s" => \(my $case = ""),
	"minread|mr=s" => \(my $minreads = 10),
        "bgtype|bg=s" => \(my $bgtype = "all"),
	"rmfile|r=s" => \(my $rmfile = ""),
	"aggdir|ad=s" => \(my $aggdir = "out_iso"),
	"aggdir2|ad2=s" => \(my $aggdir2 = ""),
	"aggdir3|ad3=s" => \(my $aggdir3 = ""),
	"mincut|mc=s" => \(my $mincut = 0.01),
	"scalecut|sc=s" => \(my $scalecut = 5),
	"quantile|q=s" => \(my $quantile = 0.9), 
	"help|h" => sub { HelpMessage(0) },
) or HelpMessage(1);
#

my %hashSR = ();
my %vafSR = ();
my %countSR = ();
my @cells;
my %bghash = ();
my %bgvafs = (); # background vafs
my %bgcohort = ();
my %rmcohort = ();

my @aggdirs = ( $aggdir );
if( $aggdir2 ne "" ) { push @aggdirs, $aggdir2; }
if( $aggdir3 ne "" ) { push @aggdirs, $aggdir3; }

if( $bgtype ne "all" ) {
  open( FI, $bgtype ) or die $!;
  while(<FI>) {
    chomp;
    $bgcohort{$_} = 1;
  }
  close FI;
}

if( $rmfile ne "" ) {
  open( FI, $rmfile ) or die $!;
  while(<FI>) {
    chomp;
    $rmcohort{$_} = 1;
  }
  close FI;
}

my @cases = ();
my %cutoffs = ();
my %wt = ();

for( my $i=1; $i<=50; $i++ ) {
  $wt{$gene . "_Exon" . $i . ":Exon" . ($i+1)} = 1;
}

for $aggdir (@aggdirs) {
  opendir( DH, $aggdir . "/" . $gene ) or die $!;
  my @files = readdir DH;
  close DH;
  my @casefiles = grep( /_vafs.txt$/, @files);

  foreach my $casefile (sort @casefiles) {
    open(FI, $aggdir . "/" . $gene . "/" . $casefile ) or die $!;
    $casefile =~ s/_${gene}_vafs.txt$//; #print $casefile . "\n";
    push @cases, $casefile;
    while(<FI>) {
      chomp;
      @cells = split(/\t/, $_);
      $hashSR{$cells[1]}{$casefile} = $cells[3];
      $vafSR{$casefile}{$cells[1]} = $cells[3];
      $countSR{$casefile}{$cells[1]} = $cells[4];
    }
    close FI;
    if( !exists( $rmcohort{$casefile} ) && $bgtype ne "all" && 
	exists( $bgcohort{$casefile} ) && !exists( $bgvafs{$casefile} ) ) {
      foreach my $isokey (keys %{$vafSR{$casefile}}) { 
        $bgvafs{$casefile}{$isokey} = $vafSR{$casefile}{$isokey};
        $bghash{$isokey}{$casefile} = $vafSR{$casefile}{$isokey};
      }
    }
  }
}


my $cutN = int( scalar( @cases ) * (1-$quantile) );
if( $bgtype eq "all" ) { 
  %bghash = %hashSR;
} else {
  $cutN = int( scalar( keys %bgcohort ) * (1-$quantile) );
}

foreach my $isokey ( keys %hashSR ) {
  if( !exists($wt{$isokey}) ) { $cutoffs{$isokey} = 0; }
}

foreach my $isokey ( keys %bghash ) {
  if( !exists($wt{$isokey}) ) {
    my @ss = sort {$bghash{$isokey}{$b} <=> $bghash{$isokey}{$a}} keys %{$bghash{$isokey}};
    my @counts = sort {$b <=> $a} values %{$bghash{$isokey}};
    my $n = scalar(@cases) - scalar(@counts);
    push @counts, (0) x $n;
    $cutoffs{$isokey} = min( 1, $counts[$cutN] * $scalecut ); 
  }
}

my %outliers = ();
foreach my $casefile (sort keys %vafSR) {
  foreach my $isokey ( sort keys %{$vafSR{$casefile}} ) {
    if( exists($cutoffs{$isokey}) && $vafSR{$casefile}{$isokey} >= max($mincut, $cutoffs{$isokey}) &&
	$countSR{$casefile}{$isokey} >= $minreads && $cutoffs{$isokey} < 1 )  {
      $outliers{$casefile}{$isokey} = $vafSR{$casefile}{$isokey};
      @cells = split(/_/, $casefile);
      if( $cells[0] eq $case ) {
        print $casefile . "\t" . $isokey . "\t" . $countSR{$casefile}{$isokey} . "\t" . 
	  $vafSR{$casefile}{$isokey} . "\t" . $cutoffs{$isokey} . "\n";
      }
    }
  }
}

