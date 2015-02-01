#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

my $muscle;
my $rmindel;
my $length;
my $help;

if (@ARGV > 0) {
	GetOptions (
		'muscle:s' => \$muscle,
		'length:s' => \$length,
		'rmindel:s' => \$rmindel,
		'help' => \$help
	);
}

if ($help) {
	print "A wrapper to glob up MUSCLE alignment files, and run the indel_remover_vdna.pl script.\n";
	print "The outfiles from the indel_remover_vdna.pl script are then globbed-up and PAML's yn00 program is executed.\n";
	print "Usage: ./indel_remover_vdna_wrapper.pl --muscle=\"name of MUSCLE DNA alignment fasta file EXTENSION to glob for indel_remover_vdna.pl script\" --length=\"length (bp) requirement of final concatenated sequence for indel_remover_vdna.pl script\" --rmindel=\"name of indel-removed DNA alignment fasta file EXTENSION to glob for PAML analysis\"\n";
	exit;
}

while (defined(my $file = glob "\*$muscle")) {
	open my $fh, "<", $file;
	foreach (defined( my $line = <$fh> )) {		
		my $ext_len = length($muscle);
		my $core = substr ($file, 0, -($ext_len));

		system("perl indel_remover_vdna.pl --dna_alignment=$file --length=$length --outfile=$core$rmindel");
	}
	close $fh;
}

system("find . -size 0 -delete");

while (defined(my $file = glob "\*$rmindel")) {
	open my $fh, "<", $file;
	foreach (defined( my $line = <$fh> )) {
		my $ext_len = length($muscle);
		my $core = substr ($file, 0, -($ext_len));
		
		open CTL, ">$core\.ctl", or die;

		print CTL "seqfile = $file\n";
		print CTL "outfile = $core\_paml.out\n";
		print CTL "verbose = 0\n";
		print CTL "icode = 0\n";
		print CTL "weighting = 0\n";
		print CTL "commonf3x4 = 0\n";
	
		system("/usr/local/paml/yn00 $core\.ctl");
		system("mv 2YN.dN $core\_2YN.dN.out");
		system("mv 2YN.dS $core\_2YN.dS.out");
		system("mv 2YN.t $core\_2YN.t.out");
		system("mv rst $core\_rst.out");
		system("mv rst1 $core\_rst1.out");
		system("mv rub $core\_rub.out");
	}
	close $fh;
}
