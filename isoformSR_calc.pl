#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum);
use List::Util qw(max);
use List::Util qw(min);
use List::MoreUtils qw(first_index indexes);
use Data::Dumper qw(Dumper);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);
use Getopt::Long qw(HelpMessage);
use Pod::Usage;

=pod

=head1 NAME

isoformSR_calc	isoforms by gsp from split reads in the bam

=head1 SYNOPSIS

  --bam, -b     Input bam
  --inpath, -i  Input path
  --outpath, -o Output path
  --annpath, -a Annotation files path
  --tmppath, -t Tmp output path
  --gene, -g	Gene
  --gtf, -gtf	GTF of gsp probes for Archer assays
  --bamx, -bx   bam suffix
  --dets, -d	Save dets
  --cutoff, -c	Cutoff for novel breakpoints
  --help, -h    Print this help

=head1 VERSION

    1.0

=cut

GetOptions(
        "bam|b=s" => \(my $bam = ""),
        "inpath|i=s" => \(my $inpath = ""),
        "outpath|o=s" => \(my $outpath = "out_iso"),
        "annpath|a=s" => \(my $annpath = "helper"),
        "tmppath|t=s" => \(my $tmppath = "out_tmp"),
	"gene|g=s" => \(my $genetarget = "IKZF1"),
	"gtf|gtf=s" => \(my $gtf = ""),
	"bamx|bx=s" => \(my $bamsuffix = ".bam" ),
        "dets|d=s" => \(my $dets = "none" ),  # can be "all", "split", or "nonsplit"
	"cutoff|c=s" => \(my $cutoff = 100 ),
	"help|h" => sub { HelpMessage(0) },
) or HelpMessage(1);

sub parse_cigardet {
    my ($text) = @_; 
    my @line_cells = split( /\_/, $text );
    my $cigar = $line_cells[0]; 
    my $editdist = $line_cells[1];
    my $rcflag = "+"; if( ($line_cells[3] & 16)==16 ) { $rcflag = "-"; }
    my $sizeleftclip = 0; my $sizerightclip = 0;
    if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $sizeleftclip = $1; }
    if( $cigar =~ /(\d+)([SH])$/ ) { $sizerightclip = $1; }
    my $lpos = $line_cells[2]; my $rpos = $line_cells[2] - 1; my $sizeseq = 0;
    while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
	$cigar = $3; 
	if( $2 ne "D" ) { $sizeseq += $1; }
	if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
    }
    $cigar = $line_cells[0];
    if( $cigar =~ /^(\d+)([SH])(.*)/ ) { $cigar = $3; }
    if( $cigar =~ /(.*)([MDI])(\d+)([SH])$/ ) { $cigar = $1 . $2; }
    return( $line_cells[0] . "_" . ( $sizeleftclip + 1 ) . "_" . ($sizeseq - $sizerightclip) . "_" . $lpos . "_" . $rpos . "_" . $cigar . "_" . $editdist . "_" . $rcflag );
}


my $proc_splits = 1;
my $proc_nonsplits = 1;
my $proc_from_bam = 1;
my $proc_archive = 0;

my $syscmd;
my $chr; my $start = -1; my $end = -1;
my $key; my $readname; my $cigar; my $seq; my $qa; my $i; my $j; my $k;
my $alnstart; my $alnend; my $rloc; my $lloc; my $editdist; my $rcflag;
my @line_cells; my @sub_cells; my @sub_cells1; 
my @p_cells; my @s_cells; my @a; my @b;
my %bksKeep = ();

my $samcmd = "samtools";
my $gene;
my @genes = ($genetarget);
my %geneChrs = map { $_ => "*" } @genes;
my %geneStrands = map { $_ => "*" } @genes;
my %geneStarts = map { $_ => -1 } @genes;
my %geneEnds = map { $_ => -1 } @genes;
my %geneRegions = map { $_ => "" } @genes;
my %transcripts = map { $_ => "*" } @genes;
my %numExons = map { $_ => "*" } @genes;
my %cstrands = map { $_ => "*" } @genes;

open(FI, $annpath . "/helper_bksKeep.txt" ) or die $!;
while(<FI>) {
    chomp;
    $bksKeep{$_} = 1;
}
close FI;

open(FI, $annpath . "/helper_genes.txt" ) or die $!;
while(<FI>) {
  chomp;
  @line_cells = split(/,/, $_);
  $gene = $line_cells[0];
  if( exists( $transcripts{$gene} ) ) {
    @sub_cells = split( /\./, $line_cells[1] );
    my $tt = $sub_cells[0]; $tt =~ s/\_//;
    $transcripts{$gene} = $tt;
  }
}
close FI;

if( $genetarget eq "KMT2A" ) { $transcripts{"KMT2A"} = "NM005933"; }
if( $genetarget eq "ERG" ) { $transcripts{"ERG"} = "NM001136154"; }

if( not -e "${annpath}/tx_stored/tx_${genetarget}.txt" ) {
  if ( system( "awk -v gene=$genetarget -F '\t' '{ if(\$13 == gene) { print } }' ${annpath}/tx_ncbiRefSeq_hg19_20220316.txt > ${annpath}/tx_stored/tx_${genetarget}.txt" ) ) { die "Failed awk extract from transcript file. Exiting..."; }
  if( system( "awk -v gene=$genetarget -F '\t' '{ if(\$13 == gene) { print } }' ${annpath}/tx_extra.txt >> ${annpath}/tx_stored/tx_${genetarget}.txt" ) ) { die "Failed awk extract from transcript file. Exiting..."; }
  if( system( "awk -v gene=$genetarget -F '\t' '{ if(\$13 == gene) { print } }' ${annpath}/tx_gencode_v28lift_hg19.txt >> ${annpath}/tx_stored/tx_${genetarget}.txt" ) ) { die "Failed awk extract from transcript file. Exiting..."; }
}

my %exonStarts = (); my %exonEnds = ();
my %exonStartsMain = (); my %exonEndsMain = ();
my %startsByExon = (); my %endsByExon = ();

open(FI, "${annpath}/tx_stored/tx_${genetarget}.txt" ) or die $!;
while(<FI>) {
    chomp;
    @line_cells = split(/\t/, $_);
    @sub_cells = split( /\./, $line_cells[1] );
    my $tt = $sub_cells[0]; $tt =~ s/\_//;
    $gene = $line_cells[12];
    if( ( $gene eq "P2RY8" || $gene eq "CRLF2" ) && $line_cells[2] eq "chrY" ) { next; }
    if( $gene eq "NOTCH4" && $line_cells[2] ne "chr6" ) { next; }
    if( exists( $transcripts{$gene} ) && $tt eq $transcripts{$gene} ) {
	$numExons{$gene} = $line_cells[8]; 
	$cstrands{$gene} = $line_cells[3];
    }
    if( exists( $geneStarts{$gene} ) ) {
	if( $geneChrs{$gene} ne "*" && $geneChrs{$gene} ne $line_cells[2] ) { print( "Chromosome change for " . $gene . "\n" ); exit(1); } 
	if( $geneStarts{$gene} < 0 || $line_cells[4] < $geneStarts{$gene} ) { $geneStarts{$gene} = $line_cells[4]; }
	if( $geneEnds{$gene} < 0 || $line_cells[5] > $geneEnds{$gene} ) { $geneEnds{$gene} = $line_cells[5]; }
	$geneChrs{$gene} = $line_cells[2];
	$geneStrands{$gene} = $line_cells[3];

	my $offset = 1; my $sign = 1;
	if( $line_cells[3] eq "-" ) { $offset = -$line_cells[8]; $sign = -1; } 
	@a = split( /,/, $line_cells[9] ); 
	for( my $i=0; $i<scalar(@a); $i++ ) {
	    $startsByExon{$gene}{$tt}{"Exon" . (($i+$offset) * $sign)} = $a[$i]+1;
	    $exonStarts{$gene}{$a[$i]+1}{$tt} = ($i+$offset) * $sign;
	    if( $tt eq $transcripts{$gene} ) {
		$exonStartsMain{$gene}{$a[$i]+1} = ($i+$offset) * $sign;
	    }
	}
	@a = split( /,/, $line_cells[10] ); 
	for( my $i=0; $i<scalar(@a); $i++ ) {
	    $endsByExon{$gene}{$tt}{"Exon" . (($i+$offset) * $sign)} = $a[$i];
	    $exonEnds{$gene}{$a[$i]}{$tt} = ($i+$offset) * $sign;
	    if( $tt eq $transcripts{$gene} ) {
		$exonEndsMain{$gene}{$a[$i]} = ($i+$offset) * $sign;
	    }
	}
    }	
}
close(FI);

foreach $gene ( @genes ) {
    $geneRegions{$gene} = $geneChrs{$gene} . ":" . ( $geneStarts{$gene} - 1e6 ) . "-" . ( $geneEnds{$gene} + 1e6 );
}

my %gsps = ();
my %gspsLoc = ();
my %gspsEnd = ();
my $igtf = -1;
if( $gtf ne "" ) {
  open(FI, $gtf) or die $!;
  while(<FI>) {
    chomp;
    @line_cells = split(/\t/,$_);
    if( $line_cells[0] eq $geneChrs{$genetarget} &&
        $line_cells[3] >= $geneStarts{$genetarget} && 
        $line_cells[4] <= $geneEnds{$genetarget} ) {	  
      @sub_cells = split(/\ /, $line_cells[8]);
      if( $igtf == -1 ) {
	$igtf = first_index { /_GSP2\";/ } @sub_cells;
      }
      $gsps{$sub_cells[$igtf]} = $line_cells[0] . "_" . $line_cells[6] . "_" . $line_cells[3] . "_" . $line_cells[4]; 
      if( $line_cells[6] eq "+" ) {
	$gspsEnd{$sub_cells[$igtf]}{$line_cells[6]} = $line_cells[4];
	for( my $i=-10; $i<=10; $i++ ) { $gspsLoc{ $line_cells[0]."_".$line_cells[6]."_".($line_cells[3]+$i)} = $sub_cells[$igtf]; }
      } elsif( $line_cells[6] eq "-" ) {
	$gspsEnd{$sub_cells[$igtf]}{$line_cells[6]} = $line_cells[3];
	for( my $i=-10; $i<=10; $i++ ) { $gspsLoc{ $line_cells[0]."_".$line_cells[6]."_".($line_cells[4]+$i)} = $sub_cells[$igtf]; }
      }
    }
  }
  close(FI);
}

my @gspkeys = sort keys %gsps;
push @gspkeys, "other";

my @casefiles;
if( $bam eq "" ) {
  opendir( DH, $inpath ) or die $!;
  @casefiles = readdir DH;
  @casefiles = grep( /${bamsuffix}$/, @casefiles);
  close DH;
} else {
  my @cells = split( /\//, $bam );
  @casefiles = ( $cells[scalar(@cells)-1] );
  if( scalar(@cells)==1 ) { $inpath = "./"; } else { $inpath = join( "/", @cells[0..(scalar(@cells)-2)] ); }
}

foreach my $case ( sort @casefiles ) {
    my %splicedata = (); my %countdata = ();
    my %moldets = (); 
    my %basereads = ();
    my %rnacounts = (); my %dnacounts = (); my %tnacounts = (); my %mixcounts = (); my %othercounts = ();
    my %wtcounts = (); my %altcounts1 = (); my %altcounts2 = ();
    my %dnadedups = (); my %rnadedups = (); my %mixdedups = (); my %otherdedups = (); my %tnadedups = ();
    my %wtnodedup = (); my %altnodedup1 = (); my %altnodedup2 = (); my %bkptnodedup = ();
    my %dnanodedup = (); my %inexcounts = (); my %inexcounts2 = ();
    my %dnapicard = (); my %rnapicard = (); my %mixpicard = (); my %otherpicard = (); my %tnapicard = ();
    my %wtpicard = (); my %altpicard1 = (); my %altpicard2 = ();
    my %picardFlag = ();

    my %dnadedupsGsp = (); my %rnadedupsGsp = (); my %mixdedupsGsp = (); my %otherdedupsGsp = (); my %tnadedupsGsp = ();
    my %wtnodedupGsp = (); my %altnodedup1Gsp = (); my %altnodedup2Gsp = (); my %bkptnodedupGsp = ();

    my $runkey = $case; $runkey =~ s/$bamsuffix//;
    print "***** $runkey *****\n";

    @sub_cells = split(/_/, $runkey );
    @line_cells = split(/-/, $sub_cells[0] );

    my $newrunkey = $line_cells[0] . "-" . $line_cells[1];
    if( scalar(@line_cells) > 2 ) { $newrunkey .= "-" . $line_cells[2]; }
    $newrunkey .= "_" . $sub_cells[1];

    foreach $gene ( @genes ) {
	print "************* $gene *****************\n";
	my $fbase = $tmppath . "/" . $runkey . "_" . $gene; print $fbase . "\n";
	if( $proc_from_bam ) {
	    my $tmparchivename = $inpath . "/" . $case;
	    $syscmd = sprintf( "%s view %s %s > %s",
			$samcmd, $tmparchivename, $geneRegions{$gene}, $fbase . ".sam" );
	    print($syscmd . "\n");
	    if( system($syscmd) ) { die "Failed samtools extract gene region. Exiting..."; }
	}

	my %psam = (); my %ssam = (); my %leftclipsam = (); my %rightclipsam = (); my %pseq = ();
	my $readsuffix; my %pdata = (); my %sdata = (); my %bkdata = (); my %primdata = ();
	my %masterMap = ();

	if( $genetarget ne "allgenes" && $gene ne $genetarget ) { next; }

	sleep(10);
	open( FI, $fbase . ".sam" ) or die $!;
	while(<FI>){
	    if( !/^@/ ) {
		@line_cells = split( /\t/, $_ );
		if( ( $line_cells[1] & 64 ) == 64 ) {
		    $readsuffix = "/1";
		} elsif( ( $line_cells[1] & 128 ) == 128 ) {
		    $readsuffix = "/2";
		} else {
		    $readsuffix = "";
		}

		my $mm = -1;
		if( $_ =~ /NM:i:(\d+)/ ) { $mm = $1; }
		if( ( $line_cells[1] & 256 ) != 256 && ( $line_cells[1] & 2048 ) != 2048 &&
		    !exists($bkdata{$line_cells[0] . $readsuffix}{"primary"}) ) {
		    $bkdata{$line_cells[0] . $readsuffix}{"primary"} = $line_cells[5] . "_" . $mm . "_" . $line_cells[3] . "_" . $line_cells[1];
		    $bkdata{$line_cells[0] . $readsuffix}{"rcBit"} = ($line_cells[1] & 16);
		    $picardFlag{$line_cells[0] . $readsuffix} = ($line_cells[1] & 1024);
		    $pseq{$line_cells[0] . $readsuffix} = $line_cells[9];
		} elsif( ( $line_cells[1] & 2048 ) == 2048 ) {
		    $bkdata{$line_cells[0] . $readsuffix}{"secondary"}{$line_cells[5] . "_" . $mm . "_" . $line_cells[3] . "_" . $line_cells[1] } = ($line_cells[1] & 16);
		}
	    }
	}
	close FI;
	    
	foreach( keys %bkdata ) {
	    if( !exists( $bkdata{$_}{"primary"} ) ) {
		delete $bkdata{$_}; 
	    } elsif( !exists( $bkdata{$_}{"secondary"} ) ) { 
		$primdata{$_} = $bkdata{$_};
		delete $bkdata{$_}; 
	    } else { 
	        # also check primary and secondary align to same strand
		foreach my $aln ( keys %{$bkdata{$_}{"secondary"}} ) {
		    if( $bkdata{$_}{"secondary"}{$aln} != $bkdata{$_}{"rcBit"} ) { delete $bkdata{$_}{"secondary"}{$aln}; }
		}
		if( scalar( keys %{$bkdata{$_}{"secondary"}} ) == 0 ) { delete $bkdata{$_}; }
	    }
	}

	my %R2 = ();
	my %locsR2 = ();
	my %gspsR2 = ();
	my %gspByR2 = ();
	foreach( keys %primdata ) {
	  if( /\/2$/ ) {
	    @line_cells = split(/\_/, $primdata{$_}{primary} );
	    if( $primdata{$_}{rcBit} == 0 ) {
	      $R2{$_} = $line_cells[2]; 
	      $locsR2{$line_cells[2]} += 1;
	      if( exists($gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}) ) {
		$gspsR2{$gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}} += 1;
		$gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]};
              } else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	    } elsif( $primdata{$_}{rcBit} == 16 ) {
	      $cigar = $line_cells[0]; my $rpos = $line_cells[2]-1;
	      while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
		$cigar = $3; 
		if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
	      }
	      if( exists($gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}) ) {
		$gspsR2{$gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}} += 1;
		$gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_-_".$rpos};
              } else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	    }
	  }
	}

	foreach( keys %bkdata ) {
	  if( /\/2$/ ) {
	    my $alnkey = $bkdata{$_}{primary};
	    @line_cells = split(/\_/, $alnkey ); $cigar = $line_cells[0];
	    if( $bkdata{$_}{rcBit} == 0 ) {
	      if( $cigar =~ /^(\d+)M/ ) {
	        $R2{$_} = $line_cells[2];
	        $locsR2{$line_cells[2]} += 1;
	        if( exists($gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}) ) {
		  $gspsR2{$gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}} += 1;
		  $gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]};
		} else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	      } else {
		foreach $alnkey (keys %{$bkdata{$_}{secondary}}) {
		  @line_cells = split(/\_/, $alnkey ); $cigar = $line_cells[0];
		  if( $cigar =~ /^(\d+)M/ ) {
		    $R2{$_} = $line_cells[2];
		    $locsR2{$line_cells[2]} += 1;
		    if( exists($gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}) ) {
		      $gspsR2{$gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]}} += 1;
		      $gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_+_".$line_cells[2]};
		    } else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
		  }
		}
		if( !exists($gspByR2{$_}) ) { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	      }
	    } elsif( $bkdata{$_}{rcBit} == 16 ) {
	      if( $cigar =~ /M$/ ) {
	        my $rpos = $line_cells[2]-1;
	        while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
		  $cigar = $3; 
		  if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
	        }
	        if( exists($gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}) ) {
		  $gspsR2{$gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}} += 1;
		  $gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_-_".$rpos};
		} else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	      } else {
		foreach $alnkey (keys %{$bkdata{$_}{secondary}}) {
		  @line_cells = split(/\_/, $alnkey ); $cigar = $line_cells[0];
		  if( $cigar =~ /M$/ ) {
		    my $rpos = $line_cells[2]-1;
		    while( $cigar =~ /(\d+)([SHDIMN])(.*)/ ) {
		      $cigar = $3; 
		      if( $2 eq "M" || $2 eq "N" || $2 eq "D" ) { $rpos += $1; }
		    }
		    if( exists($gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}) ) {
		      $gspsR2{$gspsLoc{$geneChrs{$genetarget}."_-_".$rpos}} += 1;
		      $gspByR2{$_} = $gspsLoc{$geneChrs{$genetarget}."_-_".$rpos};
		    } else { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
		  }
		}
		if( !exists($gspByR2{$_}) ) { $gspsR2{other} += 1; $gspByR2{$_} = 'other'; }
	      }
	    }
	  }
	}

	if( $proc_nonsplits ) {
	    my $primbin;
	    my @starts = sort { $a <=> $b } keys %{$exonStartsMain{$gene}};
	    my @ends = sort { $a <=> $b } keys %{$exonEndsMain{$gene}};
	    if( $dets eq "all" || $dets eq "nonsplit" ) {
	      if( !-d $outpath . "/dets/" . $gene ) { make_path $outpath . "/dets/" . $gene or die "Failed to create dets path"; }
	      open( FO, ">" . $outpath . "/dets/" . $gene . "/" .
	        $newrunkey . "_" . $gene . "_nonsplitdets.txt");
	    } 

	    foreach( sort keys %primdata ) {
		@line_cells = split( /\_/, parse_cigardet( $primdata{$_}{"primary"} ) );
		$alnstart = $line_cells[1];
		$alnend = $line_cells[2];
		$lloc = $line_cells[3];
		$rloc = $line_cells[4];
		$cigar = $line_cells[5];
		$editdist = $line_cells[6];
		$rcflag = $line_cells[7];
		$primbin = $gene . "_" . $lloc . "_" . $rloc . "_" . $cigar . "_" . $editdist . "_" . $rcflag;

		my $rex = -1; my $rin = -1; my $lex = -1; my $lin = -1;
		for( my $e=0; $e < scalar(@ends); $e++ ) {
		    if( $lloc >= $starts[$e] && $lloc <= $ends[$e]) {
			$lex = $e;
			last;
		    } elsif($e < scalar(@ends)-1 && $lloc > $ends[$e] && $lloc < $starts[$e+1]) {
			$lin = $e;
			last;
		    }
		}
		for( my $e=0; $e < scalar(@ends); $e++ ) {
		    if( $rloc >= $starts[$e] && $rloc <= $ends[$e]) {
			$rex = $e;
			last;
		    } elsif($e < scalar(@ends)-1 && $rloc > $ends[$e] && $rloc < $starts[$e+1]) {
			$rin = $e;
			last;
		    }
		}
		my $primkey = "";
		if( $lex >= 0 ) {
		    $primkey .= "Exon" . $exonStartsMain{$gene}{$starts[$lex]};
		} elsif( $lin >= 0 ) {
		    if( $geneStrands{$gene} eq "-" ) {
			$primkey .= "Intron" . ($exonStartsMain{$gene}{$starts[$lin]}-1);
		    } else {
			$primkey .= "Intron" . $exonStartsMain{$gene}{$starts[$lin]};
		    }
		} else {
		    $primkey .= "Other";
		}

		if( $rex >= 0 ) {
		    if( ( $lex >= 0 && $lex != $rex ) || $lin >= 0 ) { $primkey .= "Exon" . $exonStartsMain{$gene}{$starts[$rex]}; }
		} elsif ( $rin >= 0) {
		    if( ( $lin >= 0 && $lin != $rin ) || $lex >= 0 ) {
			if( $geneStrands{$gene} eq "-" ) {
			    $primkey .= "Intron" . ($exonStartsMain{$gene}{$starts[$rin]}-1); 
			} else {
			    $primkey .= "Intron" . $exonStartsMain{$gene}{$starts[$rin]}; 
			}
		    }
		} else {
		    $primkey .= "Other";
		}
		if( $dets eq "all" || $dets eq "nonsplit" ) { print FO $_ . "\t" . $primbin . "\t" . $primkey . "\n"; }

		@sub_cells = split( /\_/, $primbin );
		if( $sub_cells[3] eq "*" ) { next; }
		my $baseread = $_; $baseread =~ s/\/1//g; $baseread =~ s/\/2//g; $basereads{$baseread}=1;
		if( $primkey =~ /Intron/ ) {
		    $moldets{$gene}{$_}{dna}{$primkey . "_" . $sub_cells[1] . "_" . $sub_cells[2] . "_" . $sub_cells[3] . "_" . $sub_cells[5]} = 1;
		} else {
		    my $excount = () = $primkey =~ /Exon/gi;
		    if( $excount > 1 ) {
			$moldets{$gene}{$_}{dna}{$primkey . "_" . $sub_cells[1] . "_" . $sub_cells[2] . "_" . $sub_cells[3] . "_" . $sub_cells[5] } = 1;
		    } else {
			$moldets{$gene}{$_}{tna}{$primkey . "_" . $sub_cells[1] . "_" . $sub_cells[2] . "_" . $sub_cells[3] . "_" . $sub_cells[5]} = 1;
		    }
		}
	    }
	    if( $dets eq "all" || $dets eq "nonsplit" ) { close FO; }
	}
	
	if( $proc_splits ) {	    
	  my $bkbin; my %agggood = (); my %aggother = (); my $tmpkey;
	  my @starts = sort { $a <=> $b } keys %{$exonStarts{$gene}};
	  my @ends = sort { $a <=> $b } keys %{$exonEnds{$gene}};
	  my @startsMain = sort { $a <=> $b } keys %{$exonStartsMain{$gene}};
	  my @endsMain = sort { $a <=> $b } keys %{$exonEndsMain{$gene}};
	  if( $dets eq "all" || $dets eq "split" ) {
	    if( !-d $outpath . "/dets/" . $gene ) { make_path $outpath . "/dets/" . $gene or die "Failed to create dets path"; }
	    open( FO, ">" . $outpath . "/dets/" . $gene . "/" .
		$newrunkey . "_" . $gene . "_splitdets.txt" ) or die $!;
	  }
	  foreach( sort keys %bkdata ) {
	    @line_cells = split( /\_/, parse_cigardet( $bkdata{$_}{"primary"} ) );
	    my %alnstarts = ( $line_cells[1] => 0 );
	    my @alnends = ( $line_cells[2] );
	    my @llocs = ( $line_cells[3] );
	    my @rlocs = ( $line_cells[4] );
	    my @cigars = ( $line_cells[5] );
	    my @editdists = ( $line_cells[6] );
	    my @rcflags = ( $line_cells[7] );
	    $rcflag = $line_cells[7]; # enforced earlier that rcflags of secondary must be same as primary; could change this

	    my @tmpkeys = keys %{$bkdata{$_}{"secondary"}};
	    for( my $i=0; $i < scalar(@tmpkeys); $i++ ) {
		@line_cells = split( /\_/, parse_cigardet( $tmpkeys[$i] ) );
		$alnstarts{ $line_cells[1] } = $i+1;
		push @alnends, $line_cells[2];
		push @llocs, $line_cells[3];
		push @rlocs, $line_cells[4];
		push @cigars, $line_cells[5];
		push @editdists, $line_cells[6];
		push @rcflags, $line_cells[7];
	    }
	    
	    @tmpkeys = sort {$a<=>$b} keys %alnstarts;
	    for( my $i=0; $i < scalar(@tmpkeys)-1; $i++ ) {
		my $j = $alnstarts{$tmpkeys[$i]}; my $k = $alnstarts{$tmpkeys[$i+1]};
		my $overlap = min( $alnends[$j], $alnends[$k] ) - $tmpkeys[$i+1] + 1;
		$bkbin = $gene . "_" . $rlocs[$j] . "_" . $llocs[$k] . "_overlap" . "_" . $overlap . "_" . $cigars[$j] . "_" . $cigars[$k] . "_" . $editdists[$j] . "_" . $editdists[$k] . "_" . $rcflag; 

		my $e=0; my $s=0; my $splicekey; my $resolved=0;
		if( scalar(@endsMain)>0 && scalar(@startsMain)>0 ) {
		  while( $rlocs[$j] > $endsMain[$e] && $e+1 < scalar(@endsMain) ) { $e++; }
		  if( $e > 0 && $e < scalar(@endsMain)-1 && $endsMain[$e]-$rlocs[$j] > $rlocs[$j]-$endsMain[$e-1] ) { $e = $e-1; }
		  while( $llocs[$k] > $startsMain[$s] && $s+1 < scalar(@startsMain) ) { $s++; }
		  if( $s > 0 && $s < scalar(@startsMain)-1 && $startsMain[$s]-$llocs[$k] > $llocs[$k]-$startsMain[$s-1] ) { $s = $s-1; }
		
		  if( $geneStrands{$gene} eq "+" ) {
		    $splicekey = "Exon" . $exonEndsMain{$gene}{$endsMain[$e]} . ":Exon" . $exonStartsMain{$gene}{$startsMain[$s]};
		  } else {
		    $splicekey =  "Exon" . $exonStartsMain{$gene}{$startsMain[$s]} . ":Exon" . $exonEndsMain{$gene}{$endsMain[$e]};
		  }

		  if( $overlap >= -5 && $rlocs[$j]-$endsMain[$e]+$startsMain[$s]-$llocs[$k] == $overlap ) {  # allow negative overlap of -5 (likely seq error)
		    $bkdata{$_}{"breakpoints"}{$bkbin} = $splicekey;
		    $splicedata{$gene}{$splicekey} += 1;  
		    $countdata{$gene} += 1; 
		    $agggood{$splicekey} += 1;
		    $resolved=1;
		  }
	        }
		if( !$resolved ) {
		  $e=0; while( $rlocs[$j] > $ends[$e] && $e+1 < scalar(@ends) ) { $e++; }
		  if( $e > 0 && $e < scalar(@ends)-1 && $ends[$e]-$rlocs[$j] > $rlocs[$j]-$ends[$e-1] ) { $e = $e-1; }
		  $s=0; while( $llocs[$k] > $starts[$s] && $s+1 < scalar(@starts) ) { $s++; }
		  if( $s > 0 && $s < scalar(@starts)-1 && $starts[$s]-$llocs[$k] > $llocs[$k]-$starts[$s-1] ) { $s = $s-1; }

		  my %exonLefts; my %exonRights;
		  if( $geneStrands{$gene} eq "+" ) {
		    %exonLefts = %{$exonEnds{$gene}{$ends[$e]}};
		    %exonRights = %{$exonStarts{$gene}{$starts[$s]}};
		  } else {
		    %exonRights = %{$exonEnds{$gene}{$ends[$e]}};
		    %exonLefts = %{$exonStarts{$gene}{$starts[$s]}};
		  }

		  if( exists( $exonLefts{$transcripts{$gene}} ) && exists( $exonRights{$transcripts{$gene}} ) ) {
		    $splicekey = "Exon" . $exonLefts{$transcripts{$gene}} . ":Exon" . $exonRights{$transcripts{$gene}};
		    if( $overlap >= -5 && $rlocs[$j]-$ends[$e]+$starts[$s]-$llocs[$k] == $overlap ) {  # allow negative overlap of -5 (likely seq error)
			$bkdata{$_}{"breakpoints"}{$bkbin} = $splicekey;
			$splicedata{$gene}{$splicekey} += 1; 
			$countdata{$gene} += 1; 
			$agggood{$splicekey} += 1;
		    } else {
			$tmpkey = $splicekey . "_offset_" . ($rlocs[$j]-$ends[$e]) . "_" . ($starts[$s]-$llocs[$k]) . "_" . (-$overlap) . "_" . $gene . "_" . $rlocs[$j] . "_" . $llocs[$k] . "_overlap" . "_" . $overlap;
			$bkdata{$_}{"breakpoints"}{$bkbin} = $tmpkey;
			$aggother{$tmpkey} += 1;
		    }
		  } else {
		    $splicekey = "";
		    foreach my $transcript ( sort keys %exonLefts ) {
			if( exists( $exonRights{$transcript} ) ) {
			    $splicekey = "Exon" . $exonLefts{$transcript} . ":Exon" . $exonRights{$transcript} . "_" . $transcript;
			}
		    }
		    if( length($splicekey) == 0 ) {
			@a = sort keys %exonLefts;
			@b = sort keys %exonRights;
			$splicekey = "Exon" . $exonLefts{$a[0]} . ":Exon" . $exonRights{$b[0]} . "_" . $a[0] . ":" . $b[0];
		    }

		    if( $overlap >= -5 && $rlocs[$j]-$ends[$e]+$starts[$s]-$llocs[$k] == $overlap ) { # allow negative overlap of -5 (likely seq error)
			$bkdata{$_}{"breakpoints"}{$bkbin} = $splicekey;
			$splicedata{$gene}{$splicekey} += 1;  
			$countdata{$gene} += 1; 
			$agggood{$splicekey} += 1;
		    } else {
			$tmpkey = $splicekey . "_offset_" . ($rlocs[$j]-$ends[$e]) . "_" . ($starts[$s]-$llocs[$k]) . "_" . (-$overlap) . "_" . $gene . "_" . $rlocs[$j] . "_" . $llocs[$k] . "_overlap" . "_" . $overlap;
			$bkdata{$_}{"breakpoints"}{$bkbin} = $tmpkey;
			$aggother{$tmpkey} += 1;
		    }
		  }
		}

		my $bksplice = $bkdata{$_}{"breakpoints"}{$bkbin};
		if( $dets eq "all" || $dets eq "split" ) { print FO $_ . "\t" . $bkbin . "\t" . $bksplice . "\t" . $pseq{$_} . "\n"; }

		@line_cells = split( /\t/, $_ );
		@sub_cells = split( /\_/, $bksplice );
		@sub_cells1 = split( /\_/, $bkbin );
		my $baseread = $_; $baseread =~ s/\/1//g; $baseread =~ s/\/2//g; $basereads{$baseread}=1;
		if( scalar( @sub_cells ) == 1 ) {
		    ( my $exons = $sub_cells[0] ) =~ s/Exon//g;
		    @a = split( /\:/, $exons );
		    if( $a[1]-$a[0] == 1 ) {
			$moldets{$gene}{$_}{wt}{$bksplice . "_" . $sub_cells1[5] . "_" . $sub_cells1[6] . "_" . $sub_cells1[3] . $sub_cells1[4] . "_" . $sub_cells1[9] } = 1;
		    } else {
			$moldets{$gene}{$_}{alt1}{$bksplice . "_" . $sub_cells1[5] . "_" . $sub_cells1[6] . "_" . $sub_cells1[3] . $sub_cells1[4] . "_" . $sub_cells1[9] } = 1;
		    }
		} elsif( scalar( @sub_cells ) == 2 ) {
		    $moldets{$gene}{$_}{alt2}{$bksplice . "_" . $sub_cells1[5] . "_" . $sub_cells1[6] . "_" . $sub_cells1[3] . $sub_cells1[4] . "_" . $sub_cells1[9] } = 1;
		} else {
		    $moldets{$gene}{$_}{other}{$bksplice . "_" . $sub_cells1[5] . "_" . $sub_cells1[6] . "_" . $sub_cells1[3] . $sub_cells1[4] . "_" . $sub_cells1[9] } = 1;
		}
	    }
	  }
	  if( $dets eq "all" || $dets eq "split" ) { close FO; }
	}

	if( !$proc_archive ) { system( "rm " . $fbase . ".sam" ); } 

	my $tmpout1 = ""; my $tmpout2 = ""; my $tmpout3 = "";
	foreach my $baseread ( keys %basereads ) {
	    if( !exists( $moldets{$gene}{$baseread . "/1"} ) || !exists( $moldets{$gene}{$baseread . "/2"} ) ) { next; }
	    my $dna = 0; my $tna = 0; my $wt = 0; my $alt1 = 0; my $alt2 = 0;  my $other = 0;
	    if( exists($moldets{$gene}{$baseread . "/1"}{dna}) || exists($moldets{$gene}{$baseread . "/2"}{dna}) ) { $dna = 1; }
	    if( exists($moldets{$gene}{$baseread . "/1"}{tna}) || exists($moldets{$gene}{$baseread . "/2"}{tna}) ) { $tna = 1; }
	    if( exists($moldets{$gene}{$baseread . "/1"}{wt}) || exists($moldets{$gene}{$baseread . "/2"}{wt}) ) { $wt = 1; }
	    if( exists($moldets{$gene}{$baseread . "/1"}{alt1}) || exists($moldets{$gene}{$baseread . "/2"}{alt1}) ) { $alt1 = 1; }
	    if( exists($moldets{$gene}{$baseread . "/1"}{alt2}) || exists($moldets{$gene}{$baseread . "/2"}{alt2}) ) { $alt2 = 1; }
	    if( exists($moldets{$gene}{$baseread . "/1"}{other}) || exists($moldets{$gene}{$baseread . "/2"}{other}) ) { $other = 1; }

	    my $readkey = "";
	    foreach my $type ( ("dna","tna", "wt", "alt1", "alt2", "other") ) {
		if( exists($moldets{$gene}{$baseread . "/1"}{$type}) ) {
		    @a = sort keys %{$moldets{$gene}{$baseread . "/1"}{$type}};
		    foreach my $key (@a) {
			if(length($readkey)>0) { $readkey .= ";"; } 
			$readkey .= $type. "_" . $key; 
		    }
		}
	    }
	    my $readkey2 = "";
	    foreach my $type ( ("dna","tna", "wt", "alt1", "alt2", "other") ) {
		if( exists($moldets{$gene}{$baseread . "/2"}{$type}) ) {
		    @a = sort keys %{$moldets{$gene}{$baseread . "/2"}{$type}};
		    foreach my $key (@a) {
			if(length($readkey2)>0) { $readkey2 .= ";"; } 
			$readkey2 .= $type . "_" . $key; 
		    }
		}
	    }
	    $readkey .= "|" . $readkey2;

	    if( $dna && !$wt && !$alt1 && !$alt2 && !$other ) {
		$dnacounts{$gene} += 1;
		$dnadedups{$gene}{$readkey} += 1;
		$dnadedupsGsp{$gene}{$readkey}{$gspByR2{$baseread . "/2"}} += 1; 
		if( $picardFlag{$baseread . "/1"} == 0 ) { $dnapicard{$gene}{$readkey} += 1; }
	    } elsif( ( $wt || $alt1 || $alt2 ) ) { 
		$rnacounts{$gene} += 1;
		$rnadedups{$gene}{$readkey} += 1;
		$rnadedupsGsp{$gene}{$readkey}{$gspByR2{$baseread . "/2"}} += 1;
		if( $picardFlag{$baseread . "/1"} == 0 ) { $rnapicard{$gene}{$readkey} += 1; }
	    } elsif( $tna && !$dna && !$wt && !$alt1 && !$alt2 && !$other ) {
		$tnacounts{$gene} += 1;
		$tnadedups{$gene}{$readkey} += 1;
		$tnadedupsGsp{$gene}{$readkey}{$gspByR2{$baseread . "/2"}} += 1; 
		if( $picardFlag{$baseread . "/1"} == 0 ) { $tnapicard{$gene}{$readkey} += 1; }
	    }

	    #a read which is part rna and part other will be counted in both categories
	    if( $other ) {
		$othercounts{$gene} += 1;
		$otherdedups{$gene}{$readkey} += 1;
		$otherdedupsGsp{$gene}{$readkey}{$gspByR2{$baseread . "/2"}} += 1;
		if( $picardFlag{$baseread . "/1"} == 0 ) { $otherpicard{$gene}{$readkey} += 1; }
	    }
	}

	foreach my $readkey (keys %{$dnadedups{$gene}}) {
	    @line_cells = split(/\|/, $readkey);
	    @a = split(/;/, $line_cells[0]);
	    push @a, split(/;/, $line_cells[1]);
	    my %dnaseen = ();
	    foreach my $readdet (@a) {
		@b = split(/\_/, $readdet);
		$dnanodedup{$gene}{$b[1]} += $dnadedups{$gene}{$readkey};
		if( !exists($dnaseen{$b[1]}) ) {
		    $dnaseen{$b[1]}=1;
		    if( $b[1] =~ /^(Intron)(\d+)(Intron)(\d+)$/ ) {
			my $i0 = $2; my $i1 = $4;
			if( $cstrands{$gene} eq "-" ) { $i0 = $4; $i1 = $2; }
			if( $i0 >= $i1 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			for( $i = $i0; $i < $i1; $i++ ) {
			    $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
			    $inexcounts{$gene}{"Exon".($i+1)} += 1; 
			    $inexcounts{$gene}{"Exon".($i+1).":Intron".($i+1)} += 1; 
			    if( $i < $i1-1) { $inexcounts{$gene}{"Intron".($i+1)} += 1; }
			} 
		    } elsif( $b[1] =~ /^(Exon)(\d+)(Exon)(\d+)$/ ) {
			my $i0 = $2; my $i1 = $4;
			if( $cstrands{$gene} eq "-" ) { $i0 = $4; $i1 = $2; }
			if( $i0 >= $i1 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			for( $i = $i0; $i < $i1; $i++ ) {
			    $inexcounts{$gene}{"Exon".$i.":Intron".$i} += 1; 
			    $inexcounts{$gene}{"Intron".$i} += 1; 
			    $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
			    if( $i < $i1-1) { $inexcounts{$gene}{"Exon".($i+1)} += 1; }
			} 
		    } elsif( $b[1] =~ /^(Intron)(\d+)(Exon)(\d+)$/ ) {
			if( $cstrands{$gene} eq "-" ) {
			    if( $4 > $2 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			    for( $i=$4; $i<$2; $i++ ) {
				$inexcounts{$gene}{"Exon".$i.":Intron".$i} += 1; 
				$inexcounts{$gene}{"Intron".$i} += 1; 
				$inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1)} += 1; 
			    }
			    $inexcounts{$gene}{"Exon".$i.":Intron".$i} += 1; 
			} else {
			    if( $2 >= $4 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			    for( $i=$2; $i<$4-1; $i++ ) {
				$inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1).":Intron".($i+1)} += 1; 
				$inexcounts{$gene}{"Intron".($i+1)} += 1; 
			    }
			    $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
			}
		    } elsif( $b[1] =~ /^(Exon)(\d+)(Intron)(\d+)$/ ) {
			if( $cstrands{$gene} eq "-" ) {
			    if( $4 >= $2 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			    for( $i=$4; $i<$2-1; $i++ ) {
				$inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1).":Intron".($i+1)} += 1; 
				$inexcounts{$gene}{"Intron".($i+1)} += 1; 
			    }
			    $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
			} else {
			    if( $2 > $4 ) { print "Error in introns: $readkey ...\n"; exit(1); }
			    for( $i=$2; $i<$4; $i++ ) {
				$inexcounts{$gene}{"Exon".$i.":Intron".$i} += 1; 
				$inexcounts{$gene}{"Intron".$i} += 1; 
				$inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} += 1; 
				$inexcounts{$gene}{"Exon".($i+1)} += 1; 
			    }
			    $inexcounts{$gene}{"Exon".$i.":Intron".$i} += 1; 
			}
		    } else {
			$inexcounts{$gene}{$b[1]} += 1; 
		    }
		}
	    }
	}
	for( $i = 1; $i < $numExons{$gene}; $i++ ) {
	    if( !exists( $inexcounts{$gene}{"Intron".$i} ) ) { $inexcounts{$gene}{"Intron".$i}=0; }
	    if( !exists( $inexcounts{$gene}{"Exon".$i.":Intron".$i} ) ) { $inexcounts{$gene}{"Exon".$i.":Intron".$i}=0; }
	    if( !exists( $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)} ) ) { $inexcounts{$gene}{"Intron".$i.":Exon".($i+1)}=0; }
	}

	foreach my $readkey (keys %{$otherdedups{$gene}}) {
	    @line_cells = split(/\|/, $readkey);
	    @a = split(/;/, $line_cells[0]);
	    push @a, split(/;/, $line_cells[1]);
	    my %bkptsseen = ();
	    foreach my $readdet (@a) {
		@b = split(/\_/, $readdet);		
		if( $b[0] eq "other" ) {
		    my $bkptkey = "";
		    if( $b[6] eq $gene ) { $bkptkey = $b[7]."_".$b[8]."_".$b[9]."_".$b[10]; }
		    if( $b[7] eq $gene ) { $bkptkey = $b[8]."_".$b[9]."_".$b[10]."_".$b[11]; }
		    if( $gene eq "FGFR1" && length($bkptkey) > 0 ) { $bkptkey .= "_near_".$b[1]."_".$b[2]; }
		    if( $bkptkey ne "" && !exists($bkptsseen{$bkptkey}) ) { 
		    	$bkptnodedup{$gene}{$bkptkey} += $otherdedups{$gene}{$readkey}; $bkptsseen{$bkptkey}=1;
			foreach my $gsp (keys %{$otherdedupsGsp{$gene}{$readkey}} ) {
			    $bkptnodedupGsp{$gene}{$bkptkey}{$gsp} += $otherdedupsGsp{$gene}{$readkey}{$gsp};
			}
		    }
		}
	    }
	}

	foreach my $readkey (keys %{$rnadedups{$gene}}) {
	    @line_cells = split(/\|/, $readkey);
	    @a = split(/;/, $line_cells[0]);
	    push @a, split(/;/, $line_cells[1]);
	    my %exonsseenwt = (); my %exonsseenalt1 = (); my %exonsseenalt2 = ();
	    foreach my $readdet (@a) {
		@b = split(/\_/, $readdet);		
		if( $b[0] eq "wt" ) {
		    if(!exists($exonsseenwt{$b[1]}) ) { 
			$wtnodedup{$gene}{$b[1]} += $rnadedups{$gene}{$readkey};
		    	if(exists($rnapicard{$gene}{$readkey})) { $wtpicard{$gene}{$b[1]} += $rnapicard{$gene}{$readkey}; }
		    	$wtcounts{$gene}{$b[1]} += 1; $exonsseenwt{$b[1]}=1;
			foreach my $gsp (keys %{$rnadedupsGsp{$gene}{$readkey}} ) {
			    $wtnodedupGsp{$gene}{$b[1]}{$gsp} += $rnadedupsGsp{$gene}{$readkey}{$gsp};
			}
		    }
		} elsif( $b[0] eq "alt1" ) {
		    if(!exists($exonsseenalt1{$b[1]}) ) { 
			$altnodedup1{$gene}{$b[1]} += $rnadedups{$gene}{$readkey};
		    	if(exists($rnapicard{$gene}{$readkey})) { $altpicard1{$gene}{$b[1]} += $rnapicard{$gene}{$readkey}; }
		    	$altcounts1{$gene}{$b[1]} += 1; $exonsseenalt1{$b[1]}=1;
			foreach my $gsp (keys %{$rnadedupsGsp{$gene}{$readkey}} ) {
			    $altnodedup1Gsp{$gene}{$b[1]}{$gsp} += $rnadedupsGsp{$gene}{$readkey}{$gsp};
			}
		    }
		} elsif( $b[0] eq "alt2" ) {
		    if(!exists($exonsseenalt2{$b[1]."_".$b[2]}) ) {
			$altnodedup2{$gene}{$b[1]."_".$b[2]} += $rnadedups{$gene}{$readkey};
		    	if(exists($rnapicard{$gene}{$readkey})) { $altpicard2{$gene}{$b[1]."_".$b[2]} += $rnapicard{$gene}{$readkey}; }
		    	$altcounts2{$gene}{$b[1]."_".$b[2]} += 1;  $exonsseenalt2{$b[1]."_".$b[2]}=1;
			foreach my $gsp (keys %{$rnadedupsGsp{$gene}{$readkey}} ) {
			    $altnodedup2Gsp{$gene}{$b[1]."_".$b[2]}{$gsp} += $rnadedupsGsp{$gene}{$readkey}{$gsp};
			}
		    }
		}
	    }
	}

	my %countsL = (); my %countsR = (); my %bks = ();; my @txs = ();
	my $bkL; my $bkR; my $isokey; my $countSR; my $exL; my $exR; my $txL; my $txR;
	my $fres = $newrunkey . "_" . $gene . "_counts.txt";
	my $fvaf = $newrunkey . "_" . $gene . "_vafs.txt";
	if( !-d $outpath . "/" . $gene ) { make_path $outpath . "/" . $gene or die "Failed to create output gene path"; }
	open(FO, ">" . $outpath . "/" . $gene . "/" . $fres ) or die $!;
	for( my $i=1; $i<$numExons{$gene}; $i++ ) {
	    $exL = "Exon".$i; $exR = "Exon".($i+1); $isokey = $exL . ":" . $exR;
	    if( $geneStrands{$gene} eq '+' ) {
	      $bkL = $endsByExon{$gene}{$transcripts{$gene}}{$exL};
	      $bkR = $startsByExon{$gene}{$transcripts{$gene}}{$exR};
	    } else {
	      $bkL = $startsByExon{$gene}{$transcripts{$gene}}{$exL};
	      $bkR = $endsByExon{$gene}{$transcripts{$gene}}{$exR};
	    }
	    $bks{$isokey}{L} = $bkL; $bks{$isokey}{R} = $bkR;
	    if( exists($wtnodedup{$gene}{$isokey} ) ) {
	        print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" . $wtnodedup{$gene}{$isokey};
		for( my $j=0; $j<scalar(@gspkeys); $j++ ) {
		  if( exists( $wtnodedupGsp{$gene}{$isokey}{$gspkeys[$j]} ) ) {
		    $countSR = $wtnodedupGsp{$gene}{$isokey}{$gspkeys[$j]};
		    print FO "\t" . $countSR;
		    if( exists( $gspsEnd{$gspkeys[$j]}{'+'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsL{$bkL}{$isokey} += $countSR;
		      } else {
		        $countsR{$bkR}{$isokey} += $countSR;
		      }
		    } elsif( exists( $gspsEnd{$gspkeys[$j]}{'-'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsR{$bkR}{$isokey} += $countSR;
		      } else {
		        $countsL{$bkL}{$isokey} += $countSR;
		      }
		    } elsif( $gspkeys[$j] eq 'other' ) {
		      $countsL{$bkL}{$isokey} += $countSR/2;
		      $countsR{$bkR}{$isokey} += $countSR/2;
		    }
		  } else { print FO "\t0"; }
		}
	        print FO "\n";
	    }
	}
	if( exists( $altnodedup1{$gene} ) ) {
	    foreach $isokey ( sort keys %{$altnodedup1{$gene}} ) {
		@sub_cells = split(/:/, $isokey);
		$exL = $sub_cells[0]; $exR = $sub_cells[1];
		if( $geneStrands{$gene} eq '+' ) {
		  $bkL = $endsByExon{$gene}{$transcripts{$gene}}{$exL};
		  $bkR = $startsByExon{$gene}{$transcripts{$gene}}{$exR};
		} else {
		  $bkL = $startsByExon{$gene}{$transcripts{$gene}}{$exL};
		  $bkR = $endsByExon{$gene}{$transcripts{$gene}}{$exR};
		}
	        $bks{$isokey}{L} = $bkL; $bks{$isokey}{R} = $bkR;
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" . $altnodedup1{$gene}{$isokey}; 
		for( my $j=0; $j<scalar(@gspkeys); $j++ ) {
		  if( exists( $altnodedup1Gsp{$gene}{$isokey}{$gspkeys[$j]} ) ) {
		    $countSR = $altnodedup1Gsp{$gene}{$isokey}{$gspkeys[$j]};
		    print FO "\t" . $countSR;
		    if( exists( $gspsEnd{$gspkeys[$j]}{'+'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsL{$bkL}{$isokey} += $countSR;
		      } else {
		        $countsR{$bkR}{$isokey} += $countSR;
		      }
		    } elsif( exists( $gspsEnd{$gspkeys[$j]}{'-'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsR{$bkR}{$isokey} += $countSR;
		      } else {
		        $countsL{$bkL}{$isokey} += $countSR;
		      }
		    } elsif( $gspkeys[$j] eq 'other' ) {
		      $countsL{$bkL}{$isokey} += $countSR/2;
		      $countsR{$bkR}{$isokey} += $countSR/2;
		    }
		  } else { print FO "\t0"; }
		}
		print FO "\n";
	    }
	}
	if( exists( $altnodedup2{$gene} ) ) {
	    foreach $isokey ( sort keys %{$altnodedup2{$gene}} ) {
		@sub_cells1 = split(/_/, $isokey);
		@sub_cells = split(/:/, $sub_cells1[0]);
		$exL = $sub_cells[0]; $exR = $sub_cells[1];
		@txs = split(/:/, $sub_cells1[1]);
		$txL = $txs[0];
		if( scalar(@txs) == 1 ) { $txR = $txs[0]; } else { $txR = $txs[1]; }
		if( $geneStrands{$gene} eq '+' ) {
		  $bkL = $endsByExon{$gene}{$txL}{$exL};
		  $bkR = $startsByExon{$gene}{$txR}{$exR};
		} else {
		  $bkL = $startsByExon{$gene}{$txL}{$exL};
		  $bkR = $endsByExon{$gene}{$txR}{$exR};
		}
	        $bks{$isokey}{L} = $bkL; $bks{$isokey}{R} = $bkR;
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" . $altnodedup2{$gene}{$isokey}; 
		for( my $j=0; $j<scalar(@gspkeys); $j++ ) {
		  if( exists( $altnodedup2Gsp{$gene}{$isokey}{$gspkeys[$j]} ) ) {
		    $countSR = $altnodedup2Gsp{$gene}{$isokey}{$gspkeys[$j]};
		    print FO "\t" . $countSR;
		    if( exists( $gspsEnd{$gspkeys[$j]}{'+'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsL{$bkL}{$isokey} += $countSR;
		      } else {
		        $countsR{$bkR}{$isokey} += $countSR;
		      }
		    } elsif( exists( $gspsEnd{$gspkeys[$j]}{'-'} ) ) {
		      if( $geneStrands{$gene} eq '+' ) {
		        $countsR{$bkR}{$isokey} += $countSR;
		      } else {
		        $countsL{$bkL}{$isokey} += $countSR;
		      }
		    } elsif( $gspkeys[$j] eq 'other' ) {
		      $countsL{$bkL}{$isokey} += $countSR/2;
		      $countsR{$bkR}{$isokey} += $countSR/2;
		    }
		  } else { print FO "\t0"; }
		}
		print FO "\n";
	    }
	}
	if( exists( $bkptnodedup{$gene} ) ) {
	    foreach $isokey ( sort keys %{$bkptnodedup{$gene}} ) {
		if( $bkptnodedup{$gene}{$isokey} >= $cutoff || exists( $bksKeep{$gene ."_" . $isokey} ) ) {
		  @sub_cells = split(/_/, $isokey);
		  $bkL = $sub_cells[0]; $bkR = $sub_cells[1];
	          $bks{$isokey}{L} = $bkL; $bks{$isokey}{R} = $bkR;
		  print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" . $bkptnodedup{$gene}{$isokey}; 
		  for( my $j=0; $j<scalar(@gspkeys); $j++ ) {
		    if( exists( $bkptnodedupGsp{$gene}{$isokey}{$gspkeys[$j]} ) ) {
		      $countSR = $bkptnodedupGsp{$gene}{$isokey}{$gspkeys[$j]};
		      print FO "\t" . $countSR;
		      if( exists( $gspsEnd{$gspkeys[$j]}{'+'} ) ) {
			$countsL{$bkL}{$isokey} += $countSR;
		      } elsif( exists( $gspsEnd{$gspkeys[$j]}{'-'} ) ) {
			$countsR{$bkR}{$isokey} += $countSR;
		      } elsif( $gspkeys[$j] eq 'other' ) {
			$countsL{$bkL}{$isokey} += $countSR/2;
			$countsR{$bkR}{$isokey} += $countSR/2;
		      }
		    } else { print FO "\t0"; }
		  }
		  print FO "\n";
		}
	    }
	}
	close FO;

	my %vafs = ();
	foreach $isokey ( keys %bks ) {
	  $bkL = $bks{$isokey}{L};
	  $bkR = $bks{$isokey}{R};
	  my $numerator = 0;
	  if( exists( $countsL{$bkL}{$isokey} ) ) { $numerator += $countsL{$bkL}{$isokey}; }
	  if( exists( $countsR{$bkR}{$isokey} ) ) { $numerator += $countsR{$bkR}{$isokey}; }
	  my $denominator = 0;
 
	  if( exists( $countsL{$bkL} ) && scalar( keys %{$countsL{$bkL}} ) > 0 ) { 
	    $denominator += sum( values %{$countsL{$bkL}} );
	    for( my $i=1; $i<=10; $i++ ) {
	      if( exists( $countsL{($bkL+$i)} ) && scalar( keys %{$countsL{($bkL+$i)}} ) > 0 ) { 
		$denominator += sum( values %{$countsL{($bkL+$i)}} );
	      }
	      if( exists( $countsL{($bkL-$i)} ) && scalar( keys %{$countsL{($bkL-$i)}} ) > 0 ) { 
		$denominator += sum( values %{$countsL{($bkL-$i)}} );
	      }
	    }
	  }
	  if( exists( $countsR{$bkR} ) && scalar( keys %{$countsR{$bkR}} ) > 0 ) {
	    $denominator += sum( values %{$countsR{$bkR}} );
	    for( my $i=1; $i<=10; $i++ ) {
	      if( exists( $countsR{($bkR+$i)} ) && scalar( keys %{$countsR{($bkR+$i)}} ) > 0 ) { 
		$denominator += sum( values %{$countsR{($bkR+$i)}} );
	      }
	      if( exists( $countsR{($bkR-$i)} ) && scalar( keys %{$countsR{($bkR-$i)}} ) > 0 ) { 
		$denominator += sum( values %{$countsR{($bkR-$i)}} );
	      }
	    }
	  }
	  if( $denominator > 0 ) { $vafs{$isokey} = $numerator / $denominator; }
	}

	open( FO, ">" . $outpath . "/" . $genetarget . "/" . $fvaf ) or die $!;
	for( my $i=1; $i<$numExons{$gene}; $i++ ) {
	    $isokey = "Exon" . $i . ":Exon".($i+1);
	    if( exists( $vafs{$isokey} ) && $vafs{$isokey} > 0 ) {
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" .
			$bks{$isokey}{L} . "|" . $bks{$isokey}{R} . "\t" .  
			$vafs{$isokey} . "\t" . $wtnodedup{$gene}{$isokey} . "\n";
	    }
	}
	foreach $isokey ( sort keys %{$altnodedup1{$gene}} ) {
	    if( exists( $vafs{$isokey} ) ) {
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" .
			$bks{$isokey}{L} . "|" . $bks{$isokey}{R} . "\t" .  
			$vafs{$isokey} . "\t" . $altnodedup1{$gene}{$isokey} . "\n";
	    }
	}
	foreach $isokey ( sort keys %{$altnodedup2{$gene}} ) {
	    if( exists( $vafs{$isokey} ) ) {
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" .
			$bks{$isokey}{L} . "|" . $bks{$isokey}{R} . "\t" .  
			$vafs{$isokey} . "\t" . $altnodedup2{$gene}{$isokey} . "\n";
	    }
	}
	foreach $isokey ( sort keys %{$bkptnodedup{$gene}} ) {
	    if( exists( $vafs{$isokey} ) ) {
		print FO $newrunkey . "\t" . $gene . "_" . $isokey . "\t" .
			$bks{$isokey}{L} . "|" . $bks{$isokey}{R} . "\t" .  
			$vafs{$isokey} . "\t" . $bkptnodedup{$gene}{$isokey} . "\n";
	    }
	}
	close FO;
    }
}
