# isoformSR
Basic usage: perl isoformSR_calc.pl -b bamfile -g gene <br>
Output is written to iso_out/gene (2 files: vafs and counts) <br>
Requires a bam file with local alignments (uses supplementary alignments with flag bit 2048 to determine split-reads) <br>
Samtools executable should be in the $PATH (otherwise change $samcmd in isoformSR_calc.pl to the location of the samtools command) <br>
