#!/usr/bin/perl
# Hanquan Liang
# Tue Jan 11 2011
# Extract the infomation about gene locations
# perl hit_gene_loc.pl rice_all_msu.gff3 Bradi_1.0.gff2 ZmB73_5a.59_WGS.gff.gz.tmp Sb_plantGDB.gff rice_cds_vs_sl_e-10m8b1 brachypodium_cds_vs_sl_e-10m8b1 maize_cds_vs_sl_e-10m8b1 sorghum_cds_vs_sl_e-10m8b1 > geneloc4.txt

use strict;
use warnings;

my $usage = "perl gff2geneloc.pl blast.slim.out rice.gff brachy.gff sorghum.gff maize.gff";

# header
# print "#species\tchr#\tgene_id\tstart\tend\tmiddle\thit\n";


#########################
####### GFF file ########
#########################

my %loc_rice;
open GFFRICE, shift @ARGV or die $usage;
while (my $line = <GFFRICE>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @gene_loc = split /\t/, $line;
#		print "@gene_loc\n";
#		print "$line[2]\t";
#		if (($line[2] eq "gene") && ($line[8] =~ /Alias\=LOC_Os\dg\d+$/) ) {

#		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /Alias\=(LOC_Os\d+g\d+)$/) ) {
#		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /ID\=(.+);Name/) ) {
		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /Alias\=(\w+)$/) ) {
			my $middle = ($gene_loc[3]+$gene_loc[4])/2;
# specie, chr#, geneid, start, end, middle
#			push @{$loc_rice{$1}}, "rice", $gene_loc[0],$1,$gene_loc[3],$gene_loc[4],$middle;
			push @{$loc_rice{$1}}, "rice", $gene_loc[0],$1,$gene_loc[3],$gene_loc[4],$middle,0;
#			print "rice\t$gene_loc[0]\t$1\t$gene_loc[3]\t$gene_loc[4]\t$middle\n";
		}
	}
}
close GFFRICE;

my %loc_brachy;
open GFFBRACHY, shift @ARGV or die $usage;
while (my $line = <GFFBRACHY>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @gene_loc = split /\t/, $line;
		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /(Bradi\d+\w\d+)$/) ) {
			my $middle = ($gene_loc[3]+$gene_loc[4])/2;
			push @{$loc_brachy{$1}}, "brachy", $gene_loc[0],$1,$gene_loc[3],$gene_loc[4],$middle,0;
#			print "brachy\t$gene_loc[0]\t$1\t$gene_loc[3]\t$gene_loc[4]\t$middle\n";
		}
	}
}
close GFFBRACHY;

my %loc_maize;
open GFFMAIZE, shift @ARGV or die $usage;
while (my $line = <GFFMAIZE>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @gene_loc = split /\t/, $line;
		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /ID\=(.+);Name/) ) {
			my $middle = ($gene_loc[3]+$gene_loc[4])/2;
			push @{$loc_maize{$1}}, "maize", $gene_loc[0],$1,$gene_loc[3],$gene_loc[4],$middle,0;
#			print "maize\t$gene_loc[0]\t$1\t$gene_loc[3]\t$gene_loc[4]\t$middle\n";
		}
	}
}
close GFFMAIZE;

my %loc_sorghum;
open GFFSORGHUM, shift @ARGV or die $usage;
while (my $line = <GFFSORGHUM>) {
	chomp $line;
	unless ($line =~ /^#/) {
		my @gene_loc = split /\t/, $line;
###### HEADACHE: sorghum gene id ########
#		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /locus_tag\=(.+)\.\d+$/) ) {
		if (($gene_loc[2] eq "gene") && ($gene_loc[8] =~ /locus_tag\=(.+)\.\d+$/) ) {
			my $middle = ($gene_loc[3]+$gene_loc[4])/2;
			unless ( exists($loc_sorghum{$1}) ) {
				push @{$loc_sorghum{$1}}, "sorghum", $gene_loc[0],$1,$gene_loc[3],$gene_loc[4],$middle,0;
			}
#			print "sorghum\t$gene_loc[0]\t$1\t$gene_loc[3]\t$gene_loc[4]\t$middle\n";
		}
	}
}
close GFFSORGHUM;



#########################
###### BLAST file #######
#########################

open BLASTRICE, shift @ARGV or die $usage;
while (<BLASTRICE>) {
	chomp;
	my @hit_info = split /\./;
	my $hitgene = $hit_info[0];
#	my $hitgene = ($_ =~ s/\.\d+.+$//g);
	if (exists($loc_rice{$hitgene})) {
#		push @{$loc_rice{$hit_info[1]}}, 1;
		${$loc_rice{$hitgene}}[-1] += 1;
	} else {
		print "WARNING: $hitgene does not exists in rice gff file\n";
	}	
}
close BLASTRICE;

open BLASTBRACHY, shift @ARGV or die $usage;
while (<BLASTBRACHY>) {
	chomp;
	my @hit_info = split /\./;
	my $hitgene = $hit_info[0];
#	my $hitgene = $_ =~ s/\.\d+.+$//g;
	if (exists($loc_brachy{$hitgene})) {
		${$loc_brachy{$hitgene}}[-1] += 1;
	} else {
		print "WARNING: $hitgene does not exists in brachy gff file\n";
	}	
}
close BLASTBRACHY;


open BLASTMAIZE, shift @ARGV or die $usage;
while (<BLASTMAIZE>) {
	chomp;
	my @hit_info = split /\t/;
	my $hitgene = $hit_info[0];
	$hitgene =~ s/_T\d+//;
	$hitgene =~ s/_FGT/_FG/;
	if (exists($loc_maize{$hitgene})) {
		${$loc_maize{$hitgene}}[-1] += 1;
	} else {
		print "WARNING: $hitgene does not exists in maize gff file\n";
	}	
}
close BLASTMAIZE;


open BLASTSORGHUM, shift @ARGV or die $usage;
while (<BLASTSORGHUM>) {
	chomp;
	my @hit_info = split /\./;
	my $hitgene = $hit_info[0];
# because gff file contais only chromosome gene info,
# location info for scaffold genes (e.g. Sb0099s002010, Sb####s###### ) are unknown from the gff file we have
	if ( exists($loc_sorghum{$hitgene}) ) {
		${$loc_sorghum{$hitgene}}[-1] += 1;
	} elsif ($hitgene =~ /Sb\d+g\d+/) {
		print "WARNING: $hitgene does not exists in sorghum gff file\n";
	}	
}
close BLASTSORGHUM;


#########################
##### output report #####
#########################
print "#species\tchr\tgene_id\tstart\tend\tmiddle\thit\n";

for my $gene (sort keys %loc_rice) {
#	print "$gene\t";
	print join ("\t", @{$loc_rice{$gene}} );
	print "\n";
}

for my $gene (sort keys %loc_brachy) {
	print join ("\t", @{$loc_brachy{$gene}} );
	print "\n";
}

for my $gene (sort keys %loc_maize) {
	print join ("\t", @{$loc_maize{$gene}} );
	print "\n";
}

for my $gene (sort keys %loc_sorghum) {
	print join ("\t", @{$loc_sorghum{$gene}} );
	print "\n";
}
	



