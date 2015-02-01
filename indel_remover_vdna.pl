#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
no warnings 'uninitialized';

my $dna_alignment;
my $len;
my $out;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'dna_alignment:s' => \$dna_alignment,
		'length:s' => \$len,
		'outfile:s' => \$out,
		'help' => \$help
	);
}

if ($help) {
	print "DNA alignment infile must be in fasta format\n";
	print "Usage: --dna_alignment=\"name of DNA alignment fasta file\" --length=\"length (bp) requirement of final concatenated sequence - 300 bp is reasonable for most downstream applications\" --outfile=\"file name of outfile to write to – will be in fasta format\"\n";
	exit;
}

open FASTA, "<$dna_alignment" or die;
open OUT, ">$out" or die;

my %codontable = (
'AAA' => 'K', 'AAC' => 'N', 'AAG' => 'K', 'AAT' => 'N',
'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',
'AGA' => 'R', 'AGC' => 'S', 'AGG' => 'R', 'AGT' => 'S',
'ATA' => 'I', 'ATC' => 'I', 'ATG' => 'M', 'ATT' => 'I',
'CAA' => 'Q', 'CAC' => 'H', 'CAG' => 'Q', 'CAT' => 'H',
'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',
'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R',
'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L',
'GAA' => 'E', 'GAC' => 'D', 'GAG' => 'E', 'GAT' => 'D',
'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A',
'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',
'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',
'TAA' => '*', 'TAC' => 'Y', 'TAG' => '*', 'TAT' => 'Y',
'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S',
'TGA' => '*', 'TGC' => 'C', 'TGG' => 'W', 'TGT' => 'C',
'TTA' => 'L', 'TTC' => 'F', 'TTG' => 'L', 'TTT' => 'F'
);

my $seq;
my %sp1;
my %sp2;
my $c;
my %fasta;
my @header;

while (<FASTA>) {
	chomp;
	if ($_ =~ /^>(.*)/) {
		push (@header, $_);
		if ($seq) {
			foreach my $n (split //, $seq) {
				$c++;
				$sp1{$c} = $n;
			}
		}
		$seq = '';
	} 
	else {
		$seq .= $_;
	}
}

$c = 0;

foreach my $n (split //, $seq) {
	$c++;
	$sp2{$c} = $n;
}

my $header1 = $header[0];
my $header2 = $header[1];

$c = 1;

my %h1;
my %h2;

#Puts DNA bases into fragments using indels as delimiters

foreach (sort {$a<=>$b} keys %sp1) {
	if (exists $sp2{$_} && $sp2{$_} ne "-" && $sp1{$_} ne "-") {
		$h1{$c} .= $sp1{$_};
		$h2{$c} .= $sp2{$_};
	}
	elsif (exists $sp2{$_} && $sp2{$_} eq "-" || $sp1{$_} eq "-") {
		$c++;
	}	
}

#Removes stop codon positions from end of last fragment in both species

foreach (sort {$a<=>$b} (keys %h1)[-1]) {
	if (exists $h2{$_} && ($h1{$_} =~ /TAA$|TGA$|TAG$/ || $h2{$_} =~ /TAA$|TGA$|TAG$/)) {
		$h1{$_} = substr($h1{$_}, 0, -3);
		$h2{$_} = substr($h2{$_}, 0, -3);
	}
	else {
		next;
	}
}

for (keys %sp1) {
	delete $sp1{$_};
}
for (keys %sp2) {
	delete $sp2{$_};
}

my @ff;

#Puts each DNA fragment into each coding frame, and tests for stop codons in species 1

foreach (sort {$a<=>$b} keys %h1) {
	my @seq;
	my $codon;
	my $prot;
	my $dna;
	
	@seq = split('', $h1{$_});
	for my $i (0 .. $#seq) {
		$codon .= $seq[$i];
		if (($i + 1) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp1{$_."_1"} = $dna;
		push (@ff, $_."_1");
	}

	@seq = ();
	$codon = '';
	$prot = '';
	$dna = '';

	for my $i (0 .. $#seq) {
		$codon .= $seq[$i];
		if (($i + 2) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp1{$_."_2"} = $dna;
		push (@ff, $_."_2");
	}

	@seq = ();
	$codon = '';
	$prot = '';
	$dna = '';

	for my $i (0 .. $#seq) {	
		$codon .= $seq[$i];
		if (($i + 3) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp1{$_."_3"} = $dna;
		push (@ff, $_."_3");
	}
}

#Puts each DNA fragment into each coding frame, and tests for stop codons in species 2

foreach (sort {$a<=>$b} keys %h2) {
	my @seq;
	my $codon;
	my $prot;
	my $dna;
	
	@seq = split('', $h2{$_});
	for my $i (0 .. $#seq) {
		$codon .= $seq[$i];
		if (($i + 1) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp2{$_."_1"} = $dna;
		push (@ff, $_."_1");
	}

	@seq = ();
	$codon = '';
	$prot = '';
	$dna = '';

	for my $i (0 .. $#seq) {
		$codon .= $seq[$i];
		if (($i + 2) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp2{$_."_2"} = $dna;
		push (@ff, $_."_2");
	}

	@seq = ();
	$codon = '';
	$prot = '';
	$dna = '';

	for my $i (0 .. $#seq) {	
		$codon .= $seq[$i];
		if (($i + 3) % 3 == 0) {
			$prot .= $codontable{$codon};
			$dna .= $codon;
			$codon = '';
		}
	}
	if ($prot !~ /\*/ && $dna =~ /[A|T|C|G]|[a|t|c|g]/) {
		$sp2{$_."_3"} = $dna;
		push (@ff, $_."_3");
	}
}

my @uniq_ff;

@uniq_ff = uniq(@ff);

@uniq_ff = sort @uniq_ff;

my $seq1;
my $seq2;

#Combines fragments that are found in both species and in the same frame

foreach my $i (@uniq_ff) {
	if (exists $sp1{$i} && $sp2{$i}) {
		$seq1 .= $sp1{$i};
		$seq2 .= $sp2{$i};
	}
	else {
		next;
	}
}

if (length($seq1) <= $len) {
	warn "Sequences $header1 and $header2 have length $len, and are too short!\n";
}
else {
	print OUT $header1, "\n", $seq1, "\n";
	print OUT $header2, "\n", $seq2, "\n";
}

sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
