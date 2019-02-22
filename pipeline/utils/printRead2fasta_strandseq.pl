#!/usr/bin/perl -w
use strict;
use List::Util qw/shuffle/;

my $bamFolder = $ARGV[0];
my @files = glob("$bamFolder/*.bam");

#my @files_shuff = shuffle(@files);

#my @files_select = @files_shuff[0..99];

#open OUT1,'>>', "strandS_libs_NA12878_allLibs_noDups.fq" or die "Could not write to output file: $!\n";

open OUT2,'>', "strandS_libs_NA12878_list.txt" or die "Could not write to output file: $!\n";

#print "@files_select\n";
#Sample_NW150212-III_42/
#foreach my $file (@files_select) {
foreach my $file (@files) {
        my $file_id = $file;

        $file_id =~ /\/(.+_\d+)_/;
        $file_id = $1;
        warn("Working on library $file_id\n");

        print OUT2 "$file_id\n";

	open (BAM, "samtools view $file|") or die "Data processing failed";

	while(<BAM>) {
		chomp;
		my ($readName, $flag, $chr, $pos, $seq, $qual) = (split "\t", $_)[0,1,2,3,9,10];
		
		if ($chr eq "*") {
                        $chr = 'unknown';
                }

		#optional filtering criteria
		#next if $chr !~ /chr21|chr22/;
		#next if $flag&1024;
		if ($flag&64) {  #check if first in pair
		
			if ($flag&16) { #read reverse strand
				$seq =~ tr/ACGTacgt/TGCAtgca/;		
				$seq = reverse($seq);
			}

			my $ID = "@".$readName."_".$file_id."_".$flag."_".$chr."_".$pos;		

			print "$ID\n";
			print "$seq\n";
			print "+\n";
			print "$qual\n";
		}
	}

}

