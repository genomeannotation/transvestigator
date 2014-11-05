#!/usr/bin/env python
import os
import sys
from src.gff import read_gff, write_gff
from src.annotation import read_annotations, annotate_genes
from src.sequence import read_fasta
from src.transcript import Rsem
from src.transcript_builder import build_transcript_dictionary
from src.transcript_fixer import fix_transcript, fix_phase
from src.blast import read_blast
from src.rsem import read_rsem

fastapath = "transcriptome.fasta"  # Required
gffpath = "transcriptome.gff"  # Required
blastpath = "transcriptome.blastout" #Optional
rsempath = "transcriptome.rsem" # Optional
annopath = "transcriptome.anno"  # Optional
blacklistpath = "transcriptome.blacklist"  # Optional
tblpath = "transcriptome.new.tbl"
outgffpath = "transcriptome.new.gff"
outfastapath = "transcriptome.new.fsa"
outrsempath = "rsem.out"

def verify_path(path):
    if not os.path.isfile(path):
        return False
    return True

def read_transcript_blacklist(io_buffer):
    blacklist = []
    for line in io_buffer:
        transcript = line.strip("\n\t ")
        if transcript:
            blacklist.append(transcript)
    return blacklist

def main():
    path = ""
    rsems = None
    annos = None
    transcript_blacklist = None

    # Check for optional user-provided path (default is '.')
    if len(sys.argv) > 1:
        path = sys.argv[1] + "/"

    # Read required inputs; exit on failure
    print("Reading files ... ")
    with open(path + fastapath, "r") as fastafile:
        seqs = read_fasta(fastafile)
    if not seqs:
        sys.stderr.write("Error reading fasta; exiting now\n")
        exit()
    with open(path + gffpath, "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        sys.stderr.write("Error reading gff; exiting now\n")
        exit()

    # Read optional inputs; exit on failure
    if verify_path(path + blastpath):
        with open(path + blastpath) as blastfile:
            blast_hits = read_blast(blastfile)
        if not blast_hits:
            sys.stderr.write("Error reading blast output file; exiting now\n")
            exit()
    if verify_path(path + rsempath):
        with open(path + rsempath) as rsemfile:
            rsems = read_rsem(rsemfile)
        if not rsems:
            sys.stderr.write("Error reading rsem; exiting now\n")
            exit()
    if verify_path(path + annopath):
        with open(path + annopath, "r") as annofile:
            annos = read_annotations(annofile)
        if not annos:
            sys.stderr.write("Error reading annotations; exiting now\n")
            exit() 
    if verify_path(blacklistpath):
        with open(path + blacklistpath, "r") as transcript_blacklist_file:
            transcript_blacklist = read_transcript_blacklist(transcript_blacklist_file)
        if not transcript_blacklist:
            sys.stderr.write("Error reading blacklist; exiting now\n")
            exit()
    print("done.\n\n")

    # Perform transcript annotation if annotations were provided
    if annos:
        print("Annotating genes ... ")
        annotate_genes(gff.gene, annos)
        print("done.\n\n")

    # Convert all that data into Transcript objects
    print("Mapping gff data to transcripts ... ")
    transcript_dict = build_transcript_dictionary(seqs, gff.gene)
    # transcript_dict maps transcript/sequence ids to Transcript objects
    print("done.\n\n")

    # Do magical fixes on transcripts
    print("Cleaning up transcripts ... ")
    for transcript in transcript_dict.values():
        transcript.remove_contig_from_gene_id()
        transcript.fix_feature_lengths()
        transcript.create_starts_and_stops()
        fix_transcript(transcript) # removes multiple CDS features
        transcript.make_positive()
        fix_phase(transcript)
    print("done.\n\n")

    # Write RSEM info if provided
    if rsems:
        print("Adding RSEM info to transcripts...")
        for rsem in rsems:
            if rsem.transcript_id not in transcript_dict:
                sys.stderr.write("RSEM transcript "+rsem.transcript_id+" does not exist. Skipping...\n")
                continue
            transcript = transcript_dict[rsem.transcript_id]
            transcript.rsem = Rsem(rsem.tpm, rsem.fpkm, rsem.isopct)
        print("done.\n\n")

        print("Writing RSEM info...")
        with open(path + outrsempath, "w") as outrsemfile:
            outrsemfile.write("transcript_id\tnumber_of_CDSs\tcontains_complete_CDS\tTPM\tFPKM\tIsoPct\n")
            for transcript_id, transcript in transcript_dict.items():
                if (transcript_blacklist and\
                        transcript.sequence.header in transcript_blacklist) or\
                        not transcript.passes_filtering():
                    continue
                cds_count = 0
                contains_complete_cds = False
                for gene in transcript.genes:
                    for mrna in gene.mrna:
                        cds_count += len(mrna.cds)
                        if hasattr(mrna, "start_codon") and hasattr(mrna, "stop_codon"):
                            contains_complete_cds = True
                outrsemfile.write("\t".join([transcript_id, str(cds_count), str(contains_complete_cds), str(transcript.rsem.tpm), str(transcript.rsem.fpkm), str(transcript.rsem.isopct)])+"\n")
        print("done.\n\n")

    # Write .tbl file
    print("Writing .tbl file ... ")
    with open(path + tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            if (transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist) or\
                    not transcript.passes_filtering():
                continue
            tblfile.write(transcript.to_tbl())
    print("done.\n\n")

    # Write .gff file
    print("Writing new .gff file...")
    with open(path + outgffpath, "w") as outgfffile:
        for transcript in transcript_dict.values():
            if (transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist) or\
                    not transcript.passes_filtering():
                continue
            write_gff(outgfffile, transcript.genes[0])
    print("done.\n\n")

    # Write .fsa file
    print("Writing .fsa file ... ")
    with open(path + outfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if (transcript_blacklist and\
                    transcript.sequence.header in transcript_blacklist) or\
                    not transcript.passes_filtering():
                continue
            outfastafile.write(transcript.sequence.to_fasta())
    print("done.\n\n")


###########################################

if __name__ == '__main__':
    main()
