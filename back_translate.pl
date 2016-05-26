#!/usr/bin/perl

#Forces a nucleotide sequence alignment onto a protein alignment
#Need protein alignment in fasta format, and DNA sequence in fasta format
#Sequence names in protein and DNA fasta must be the same

use strict;
use warnings;
use Getopt::Long;

my $prot_seq;
my $prot_header;
my $c;
my $dna_seq;
my $dna_header;
my %dna;
my %indel;

my $prot;
my $dna;
my $out;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'prot:s' => \$prot,
		'dna:s' => \$dna,
		'outfile:s' => \$out,
		'help' => \$help
	);
}
if ($help) {
	print "All infiles (--prot and --dna) must be in fasta format)\n";
	print "Usage: --prot=\"name of protein alignment fasta file\" --dna=\"fasta file of unaligned, correspononding DNA sequences\" --outfile=\"file name of outfile to write to – will be in fasta format\"\n";
	exit;
}

open ALIGNMENT, "<$prot" or die;
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
	$c = 0;
}

foreach my $aa (split //, $prot_seq) {
	$c++;
	if ($aa eq "-") {
		push @{$indel{$prot_header}}, $c;
	}
}

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

#Inserts sequence specific indels in DNA sequences

foreach my $header (sort {$a cmp $b} keys %dna) {
	if (exists $indel{$header}) {
		foreach my $i (@{$indel{$header}}) {
			if ($i*3-3 < length($dna{$header})) {
				substr $dna{$header}, $i*3-3, 0, '---';
			}
			else {
				$dna{$header} .= '---';
			}
		}
		print OUT ">", $header, "\n", $dna{$header}, "\n";
	}
}
