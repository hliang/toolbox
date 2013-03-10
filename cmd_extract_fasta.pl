#!/usr/bin/perl
# Hanquan Liang
# Fri Oct 22 2010
# usage: ./cmd_extract_fasta.pl idsfile input.fasta > output.fasta

use warnings;
use strict;

my $idsfile = $ARGV[0];
my $seqfile = $ARGV[1];
my %ids = ();

# build a hash for read names
open IDSFILE, $idsfile or die $!;
while (<IDSFILE>) {
	chomp;
	$ids{$_} += 1;
}
close IDSFILE;

# read by FASTA record
local $/ = ">";

open FASTA, $seqfile or die $!;
while (<FASTA>) {
	unless ($_ eq ">" ) {
		chomp;
		my $seq = $_;
# get the readname of each fasta record.
		my ($id) = $seq =~ /^>*(\S+)/;
		if (exists($ids{$id})) {
#		print "$ids{$id}\n";
# remove fasta header
#		$seq =~ s/^>*.+\n//;
#remove endlines, concatenate all sequences in one line
#		$seq =~ s/\n//g;
			print ">","$seq";
		} 
	}
}
close FASTA;

