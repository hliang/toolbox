#!/usr/bin/perl
# Hanquan Liang
use strict;
use warnings;

my $usage = "perl count_completeness.pl seqB_len.txt blast_seqA_vs_seqB_m8";

##### ===== head -n 6 Ta_est_rikenJP.id.len.txt ===== #####
#RFL_Contig1	590
#RFL_Contig2	496
#RFL_Contig3	579
#RFL_Contig4	1399
#RFL_Contig5	704
#RFL_Contig6	630
##### ===== head -n 6 sl_vs_est_JP_e0.001m8 ===== #####
#scaffold00006	RFL_Contig409	86.84	114	15	0	162	275	861	748	4e-22	 107
#scaffold00009	RFL_Contig2225	97.14	35	1	0	9685	9719	3722	3688	4e-08	61.9
#scaffold00011	RFL_Contig5767	100.00	6239	0	0	10942	17180	6239	1	0.0	1.233e+04
#scaffold00011	RFL_Contig3311	100.00	6239	0	0	10942	17180	6239	1	0.0	1.233e+04
#scaffold00012	RFL_Contig3542	86.14	166	23	0	6567	6732	544	709	1e-33	 147
#scaffold00012	RFL_Contig3542	83.33	120	20	0	8205	8324	1096	1215	2e-13	79.8

##### ===== perl count_completeness.pl Ta_est_rikenJP.id.len.txt sl_vs_est_JP_e0.001m8 ===== #####

my %ests;
my $maxlen=0;

open SEQLEN, shift @ARGV or die $usage;
while(<SEQLEN>) {
	chomp;
	unless (/^#/) {
		my @line = split /\t/;
		push @{$ests{$line[0]}}, $line[1];
		$maxlen = $line[1] if ($line[1] > $maxlen) ;
	}
}


open BLAST, shift @ARGV or die $usage;
while(<BLAST>) {
	chomp;
	my @line = split /\t/;
	push @{$ests{$line[1]}}, ($line[8], $line[9]);
}

print "#id\tlength\tcov_len\tbrkpnt:\n";
for my $estid (sort keys %ests) {
	my @cov_ends = &expandbin(@{$ests{$estid}}[0..$#{$ests{$estid}}]); # slice
	print join("\t", $estid, ${$ests{$estid}}[0], $cov_ends[0], "brkpnt:", @cov_ends[1..$#cov_ends]);
#	print "$estid\t@{$ests{$estid}}";
	print "\n";
}

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}


sub expandbin {
        my @list = @_;
        my @x2y;
	my $lengthmax = shift @list;
        foreach (1..$lengthmax+2) {
                push @x2y, 0;
        }
#       my $i; my $j;
        while (@list) {
                my $i = shift @list;
                my $j = shift @list;
		($i, $j) = ($j, $i) if ($i>$j); #make sure $i(start_pos) < $j(end_pos)
#		($i, $j) = (sort { $b <=> $a } [$i, $j]);
                foreach ($i-1..$j-1) {
                        $x2y[$_] +=1;
                }
        }
#       my $max = (sort { $b <=> $a } @x2y)[0];
        my @ends; my $n=0; my $in=0; my $covered=0;
        foreach (@x2y) {
		$n++;
                if (($_>0) && ($in==0) ) { #start of match
                        push (@ends, $n);
                        $in++;
			$covered++;
                } elsif ( ($_==0) && ($in>0) ) { #end of match
                        push (@ends, $n-1);
                        $in=0;
                } elsif (($_>0) && ($in>0) && ($n==scalar@x2y) ) { #end of match at the end of sequence, (not necessary because $maxlen+2
                        push (@ends, $n);
			$in=0;
			$covered++;
                } elsif (($_>0) && ($in>0) ) { #match
			$covered++;
		}
        }
        return ($covered, @ends);
}
