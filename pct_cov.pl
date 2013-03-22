#!/usr/bin/perl

use strict;
use warnings;

#PGTG_16_coverage.allsites.txt
#PGTG_gtf_gene_exons.3.txt

# output of GATK:DepthOfCoverage
my $fh_cov = shift @ARGV || "PGTG_16_coverage.allsites.txt";
# custom file, containing start and end pos of each gene
my $fh_genepos = shift @ARGV || "PGTG_gtf_gene_exons.3.txt";


my $min_cov = 1;

my %HoHoA;
my @samples;
print STDERR "reading $fh_cov ...\n";
open COV, $fh_cov;
while(<COV>){
	print STDERR " $.\r";
	chomp;
	my @tmp=split(/\s|:/, $_);
	#header line contains sample IDs
	if($.==1){
		for(my $i=3; $i<=$#tmp; $i++) {
			$tmp[$i] =~ s/Depth_for_//;
			$samples[$i]=$tmp[$i];
		}
	}else {
		for(my $i=4; $i<=$#tmp; $i++){ #only this fields contain the depth for samples
			#sample->contig->pos
			${$HoHoA{$samples[$i-1]}{$tmp[0]}}[$tmp[1]] +=0;
			${$HoHoA{$samples[$i-1]}{$tmp[0]}}[$tmp[1]]++ if $tmp[$i]>$min_cov;
		}
	}
}
close COV;
print STDERR "\ndone\n";

print "gene_id\t", join("\t", @samples[3..$#samples]), "\n";

print STDERR "\nreading $fh_genepos ...\n";
open POS, $fh_genepos;
while(<POS>) {
	if(s/^#//){
		print STDERR " $.\r";
		chomp;
		my @tmp=split(/ /, $_);
		
		my $thisContig=$tmp[0];
		my $thisGene=$tmp[1];
		my $length = $tmp[-1]-$tmp[3]+1;
		my %HoH_cov; # sample -> gene -> percentage
		print "$thisGene";
		foreach my $thisSample (@samples[3..$#samples]) {
			for my $i ($tmp[3]..$tmp[-1]){
				#$HoH_cov{$thisSample}{$thisGene}++ if (exists(${$HoHoA{$thisSample}{$thisContig}}[$i]) && ${$HoHoA{$thisSample}{$thisContig}}[$i]>0);
				$HoH_cov{$thisSample}{$thisGene}+=0 if (exists(${$HoHoA{$thisSample}{$thisContig}}[$i]) ); #&& ${$HoHoA{$thisSample}{$thisContig}}[$i]>0);
				$HoH_cov{$thisSample}{$thisGene}++ if (exists(${$HoHoA{$thisSample}{$thisContig}}[$i]) && ${$HoHoA{$thisSample}{$thisContig}}[$i]>0);
			}
			my $thisPct = sprintf "%.3f", $HoH_cov{$thisSample}{$thisGene}/$length;
			#print "\t$HoH_cov{$thisSample}{$thisGene}/$length=$thisPct";# if (exists(@{$HoHoA{$thisSample}{$thisContig}}) );
			print "\t$HoH_cov{$thisSample}{$thisGene}/$length=$thisPct" if $#{$HoHoA{$samples[4]}{$thisContig}}>0 ; #the "if condition" is useful for splitting pos file into smaller files for parallel processing.
		}
		print "\n";
	}

}
close POS;
print STDERR "\ndone\n";

