#!/usr/bin/perl

#A script that identifies indels in a protein alignment, and removes corresponding DNA sequence
#Requires an aligned protein sequence, and corresponding unaligned DNA sequence files in fasta format
#This script can be easily be changed to print out the indel-removed protein sequences

use strict;
use warnings;
use Getopt::Long;

my $prot_seq;
my $prot_header;
my $c;
my $dna_seq;
my $dna_header;
my %dna;
my %dna_indel;
my %indel;
my @indel;
my $prot_alignment;
my $dna;
my $out;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'prot_alignment:s' => \$prot_alignment,
		'dna:s' => \$dna,
		'outfile:s' => \$out,
		'help' => \$help
	);
}
if ($help) {
	print "All infiles (--prot_alignment and --dna) must be in fasta format)\n";
	print "Usage: --prot_alignment=\"name of protein alignment fasta file\" --dna=\"fasta file of unaligned, correspononding DNA sequences\" --outfile=\"file name of outfile to write to – will be in fasta format\"\n";
	exit;
}

open ALIGNMENT, "<$prot_alignment" or die;
open DNA, "<$dna" or die;
open OUT, ">$out" or die;

#Stores indel locations of protein alignment

while (<ALIGNMENT>) {
	chomp;
	if ($_ =~ /^>(.*)/) {
		if ($prot_seq) {
			foreach my $aa (split //, $prot_seq) {
				$c++;
				if ($aa eq "-") {
					push @{$indel{$prot_header}}, $c;
				}
			}
		}
		$prot_header = $1;
		$prot_seq = '';
	} 
	else {
		$prot_seq .= $_;
	}
}

$c = 0;

foreach my $aa (split //, $prot_seq) {
	$c++;
	if ($aa eq "-") {
		push @{$indel{$prot_header}}, $c;
	}
}

#Tallies the number of indels, and prints to screen...large numbers probably means the alignment is bad

my $sum;

foreach (keys %indel) {
	my $size = @{$indel{$_}};
	$sum += $size;
}

warn "There is a total of $sum indels in the protein alignment. If this scares you, go back and check the alignment. If not, keep calm and carry on...\n";

#Stores DNA sequences

while (<DNA>) {
	chomp;
	if ($_ =~ /^>(.*)/) {
		if ($dna_seq) {
			$dna{$dna_header} = $dna_seq;
		}
		$dna_header = $1;
		$dna_seq = '';
	}
	else {
		$dna_seq .= $_;
	}
}

$dna{$dna_header} = $dna_seq;

#Inserts sequence specific indels in DNA sequences, stores locations of indels, and removes stop codons

foreach my $header (keys %dna) {
	$c = 0;
	if (exists $indel{$header}) {
		foreach my $i (@{$indel{$header}}) {
			substr $dna{$header}, $i*3-3, 0, '---';
		}
	}
	my $codon = '';
	foreach my $nt (split //, $dna{$header}) {
		$codon .= $nt;
		if (length($codon) == 3) {
			$c++;
			if ($codon eq "---") {
				push (@indel, $c);
			}
			elsif ($codon =~ /TAG|TAA|TGA/) {
				next;
			}
			else {
				$dna_inel{$header}{$c} = $codon;
			}
			$codon = '';
		}
	}
}

#Deletes (reciprocally) the indels in each DNA sequence. If all goes well you should be left with an in-frame DNA alignment

while (my ($header, $positions) = each %dna_indel) {
	while (my ($pos, $n) = each %$positions) {
		foreach my $i (@indel) {
			if ($i eq $pos) {
				delete ($dna_indel{$header}{$pos});
			}
			else {
				next;
			}
		}
	}
}

#Prints in-frame DNA alignment to outfile

foreach my $header (keys %dna_indel) {
	print OUT ">$header\n";
	foreach my $pos (sort {$a <=> $b} keys %{$dna_indel{$header}}) {
		print OUT "$dna_indel{$header}{$pos}";
	}
	print OUT "\n";
}
