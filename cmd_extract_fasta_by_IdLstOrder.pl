#!/usr/bin/perl
# Hanquan Liang
# Fri Oct 22 2010
# usage: ./cmd_extract_fasta.pl idsfile input.fasta > output.fasta

use warnings;
use strict;

my $usage = "./cmd_extract_fasta_...pl idsfile input.fasta > output.fasta";
my $idsfile = $ARGV[0];
my $seqfile = $ARGV[1];

# save "line end"
my $line_end = $/;
# read by FASTA record
local $/ = ">";

my %seqhash;
open FASTA, $seqfile or die $usage;
while (<FASTA>) {
	unless ($_ eq ">") {
	chomp;
	my $seq = $_;
# get the readname of each fasta record.
	my ($id) = $seq =~ /^>*(\S+)/;
	$seqhash{$id} = ">".$seq ; #unless ($id eq "");
#	print $seqhash{$id};
	}
}
close FASTA;

# restore "line end" to normal "\n"
local $/ = $line_end;
#$/ = "\n";
open IDSFILE, $idsfile or die $usage;
while (<IDSFILE>) {
	chomp;
	if (exists($seqhash{$_}) ) {
		print "$seqhash{$_}";
#		print length($seqhash{$_}),"\n";
	} else {
		die "sequence $_ does not exists in input seq file.\n";
	}
}
close IDSFILE;
