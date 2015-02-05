##Mani-seq

A series of scripts to manipulate DNA and protein sequences.

---

####indel_remover_vprot.pl: Removes indels in a pairwise protein sequence alignment

A script that removes indels in DNA sequences that were identified in the corresponding pairwise protein alignment. The output will be two, in-frame DNA sequences, which is particularly useful for downstream analyzes. For example, calculating pairwise non/synonymous substitution rates (*dN*, and *dS*, respectively), and adaptive evolution (*omega* = *dN/dS*). This script requires two files in **fasta** format: (i) Pairwise **protein alignment**, and (ii) Corresponding **unaligned DNA** sequences. I have provided example infiles from a reciprocally best, BLAST test I did between *A. thaliana* and *O. sativa*. I performed the protein alignment using [MUSCLE](http://www.drive5.com/muscle/index.htm). I opted to go from a protein alignment, and essentially back translate to DNA, because often I am comparing species that are distantly diverged from one another, and many silent changes could occur causing messy DNA alignments. Also, it is easier to remove indels in protein alignments and keep the DNA sequence in-frame.

A word of caution: Alignment programs don't have a "biological intuition", so it is a good idea to manually check your alignments from time to time. I have added a warning message in the script that prints to screen the cumulative sum of indels. Remember the indels are in the protein sequence, so the number of nucleotides missing is the indel number multiplied by three.

Run as so:

````bash
./indel_remover_vprot.pl --prot_alignment= --dna= --outfile=
````

---

####indel_remover_vdna.pl: Removes indels in a pairwise DNA sequence alignment

A similar script to indel_remover_vprot.pl, however only a DNA pairwise DNA sequence alignment (in fasta format) is required. It is not a pretty script, and could probably use a couple subroutines, but it works. This script has more steps because of framing issues, and these are worth mentioning:

* Indels are identified and removed reciprocally from both sequences
* While removing, the sequences are split into fragments based on the location of these indels
* Each fragment, for each species, is shifted into three different frames, and a logical statement is used to retain sequences that do not have stop codons
* If fragments in the same frame are found in both species they are kept, and concatenated into a single sequence

I have added a --length flag, so you can specify a length (bp) cutoff for your concatenated sequences. How short do you go? The general consensus is no less than 300bp (100aa) for downstream applications like phylogenetic tree construction and substitution rate estimates.

Run as so:

````bash
./indel_remover_vdna.pl --dna_alignment= --length= --outfile=
````

---

###indel_remover_vdna_wrapper.pl: Wrapper to run indel_remover_vdna.pl and PAML

A wrapper to glob up MUSCLE alignment files, and run the indel_remover_vdna.pl script. The outfiles from the indel_remover_vdna.pl script are then globbed-up and PAML's yn00 program for pairwise sequence alignment is executed. You need to provide the name of MUSCLE DNA alignment fasta file **EXTENSION** to glob for indel_remover_vdna.pl script, the length (bp) requirement of final concatenated sequence for indel_remover_vdna.pl script, and the name of indel-removed DNA alignment fasta file **EXTENSION** to glob for PAML analysis. All MUSCLE alignments must be in the same folder.

Run as so:

````bash
./indel_remover_vdna_wrapper.pl --muscle= --length= --rmindel=
````

---

###rscu.pl: Estimates Relative Synonymous Codon Usage (RSCU)

![equation](http://www.sciweavers.org/tex2img.php?eq=RSCU_%7Bi%7D%20%3D%20%20%20%5Cfrac%7BX_%7Bi%7D%7D%7B%5Cfrac%7B1%7D%7Bn%7D%20.%20%5Csum_i%5En%20%20X_%7Bi%7D%7D%20%20%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)