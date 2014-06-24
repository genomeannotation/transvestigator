transvestigator
===============

Validates transcriptome and prepares it for submission to the NCBI!!!  Does a lot of things, but we haven't finished documenting it yet.  Contact us if you want help using and get your transcriptome submitted to TSA!

## Basic Usage

Type

    python3 transvestigator.py

in a directory containing files called "transcriptome.fasta" and "transcriptome.gff".

Optionally, you may include a file called "transcriptome.anno" -- this is a tab-separated file of annotations in the format specified [here](http://genomeannotation.github.io/GAG/#functional-annotations).

If you include a file named "transcriptome.blacklist", the transcripts listed in that file will be excluded from your output.

After the program runs, you will find files called "transcriptome.new.fsa" and "transcriptome.new.tbl". These are your input files for [tbl2asn](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2).


## What It Does:

When you run transvestigator, it reads the transcriptome into memory. It then fixes feature lengths, creates starts and stops, removes multiple CDS features, adjusts CDS phase, then writes a .fasta and .tbl file. If the directory contained a "transcriptome.anno" file, the .tbl file will contain functional annotations. If the directory contained a "transcriptome.blacklist" file, those transcripts will be excluded.

#### Fix Feature Lengths

If the indices given for a gene, mRNA or CDS extend beyond the actual length of the sequence which contains them, the end index of the feature is adjusted to fall within the sequence boundaries.


#### Create Starts And Stops

The .tbl file indicates the presence of start and stop codons, but the .gff file associated with a transcriptome may not explicitly indicate the presence or absence of these features. So transvestigator inspects the first/last three bases of each CDS to determine whether they are a start/stop, then updates the transcriptome accordingly. This step is the key to eliminated the dreaded "PartialProblem" errors further downstream.


#### Remove Multiple CDS Features

If a transcript contains multiple CDS features, the longest one is chosen and the others are discarded. Future work may involve more complicated algorithms for making this decision.


#### Adjust CDS Phase

A strange error arises in NCBI TSA submissions when a CDS begins at the second or third base of a sequence. The fix we apply changes the start index of the feature to 1 and adjusts its phase in order to compensate. This adds a "codon_start" annotation to the .tbl file and, more importantly, makes errors go away :)

