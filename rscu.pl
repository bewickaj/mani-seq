#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $fasta;
my $complete;
my $out;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'fasta:s' => \$fasta,
		'complete:s' => \$complete,
		'outfile:s' => \$out,
		'help' => \$help
	);
}

if ($help) {
	print "Usage: --fasta=\"name of fasta file\" --complete=\"yes: the majority of your sequences start with ATG and end with TAA, TGA, or TAG; no: most of your sequences are partially complete, but still in-frame\" --outfile=\"file name of outfile to write to\"\n";
	exit;
}

open FASTA, "<$fasta" or die;
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

my %syncount = (
'F' => 2, 'L' => 6, 'I' => 3, 'M' => 1,
'V' => 4, 'S' => 6, 'P' => 4, 'T' => 4,
'A' => 4, 'Y' => 2, 'H' => 2, 'Q' => 2,
'N' => 2, 'K' => 2, 'D' => 2, 'E' => 2,
'C' => 2, 'W' => 1, 'R' => 6, 'G' => 4,
'*' => 3
);

my $seq;
my $header;
my %fasta;

while (<FASTA>) {
	chomp;
	if ($_ =~ /^>(.*)/) {
		if ($seq) {
			$header =~ s/ .*//g;
			$fasta{$header} = $seq;
		}
		$header = $1;
		$seq = '';
	} 
	else {
		$seq .= $_;
	}
}
$header =~ s/ .*//g;
$fasta{$header} = $seq;

print OUT "locus\taa\tcodon\tsyn\tobs\ttotal\tRSCU\n";

foreach (keys %fasta) {
	my %hoh;
	if ($complete eq "yes") {
		if ($fasta{$_} =~ /^ATG/ && $fasta{$_} =~ /[TAA|TAG|TGA]$/) {
			my $codon;
			my @seq = split('', $fasta{$_});
			for my $nt (@seq) {
				$codon .= $nt;
				if (length($codon) == 3) {
					if ($codon =~ /[W|S|M|K|R|Y|B|D|H|V|N]/) {
						next;
					}
					else {
						$hoh{$codontable{$codon}}{$codon}++;
						$codon = '';
					}
				}
			}
			foreach my $aa (keys %hoh) {
				my $sum;		
				foreach my $codon (keys %{$hoh{$aa}}) {
					$sum += $hoh{$aa}{$codon};
				}
				foreach my $codon (keys %{$hoh{$aa}}) {
					my $rscu = $hoh{$aa}{$codon} / ((1 / $syncount{$aa}) * ($sum));
					print OUT "$_\t$aa\t$codon\t$syncount{$aa}\t$hoh{$aa}{$codon}\t$sum\t$rscu\n";
				}
			}
		}
	}
	elsif ($complete eq "no") {
		my $codon;
		my @seq = split('', $fasta{$_});
		for my $nt (@seq) {
			$codon .= $nt;
			if (length($codon) == 3) {
				if ($codon =~ /[W|S|M|K|R|Y|B|D|H|V|N]/) {
					next;
				}
				else {
					$hoh{$codontable{$codon}}{$codon}++;
					$codon = '';
				}
			}
		}
		foreach my $aa (keys %hoh) {
			my $sum;		
			foreach my $codon (keys %{$hoh{$aa}}) {
				$sum += $hoh{$aa}{$codon};
			}
			foreach my $codon (keys %{$hoh{$aa}}) {
				my $rscu = $hoh{$aa}{$codon} / ((1 / $syncount{$aa}) * ($sum));
				print OUT "$_\t$aa\t$codon\t$syncount{$aa}\t$hoh{$aa}{$codon}\t$sum\t$rscu\n";
			}
		}
	}		
}
