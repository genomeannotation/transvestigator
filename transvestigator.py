#!/usr/bin/env python3
# coding=utf-8
import os
import sys
import argparse
from src.gff import read_gff, write_gff
from src.annotation import read_annotations, annotate_genes
from src.sequence import Sequence, read_fasta
from src.sequtil import get_subsequence, translate
from src.transcript import Rsem
from src.transcript_builder import build_transcript_dictionary
from src.blast import read_blast
from src.rsem import read_rsem

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
    # Handle command line args
    parser = argparse.ArgumentParser(
    epilog="""
    Docs at http://genomeannotation.github.io/transvestigator
    Bugs and feature requests at https://github.com/genomeannotation/transvestigator/issues
    """,
    formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-f', '--fasta', required=True)
    parser.add_argument('-g', '--gff', required=True)
    parser.add_argument('-a', '--anno')
    parser.add_argument('-r', '--rsem')
    parser.add_argument('-bl', '--blacklist')
    parser.add_argument('-o', '--out')
    parser.add_argument('--blast')
    args = parser.parse_args()

    # Make output directory
    out = "transvestigator_out"
    if args.out:
        out = args.out
    os.system('mkdir ' + out)


    tblpath = out + "/transcriptome.new.tbl"
    statspath = out + "/transcriptome.new.stats"
    outgffpath = out + "/transcriptome.new.gff"
    outfastapath = out + "/transcriptome.new.fsa"
    outcdsfastapath = out + "/transcriptome.new.cds.fasta"
    outcdspeppath = out + "/transcriptome.new.cds.pep"
    outrsempath = out + "/rsem.out"
    rsemerrpath = out + "/rsem.err"

    rsems = None
    annos = None
    transcript_blacklist = None

    # Read required inputs; exit on failure
    print("Reading files ... ")
    with open(args.fasta, "r") as fastafile:
        seqs = read_fasta(fastafile)
    if not seqs:
        sys.stderr.write("Error reading fasta; exiting now\n")
        exit()
    with open(args.gff, "r") as gfffile:
        gff = read_gff(gfffile)
    if not gff:
        sys.stderr.write("Error reading gff; exiting now\n")
        exit()

    # Read optional inputs; exit on failure
    if args.blast and verify_path(args.blast):
        with open(args.blast, 'r') as blastfile:
            blast_hits = read_blast(blastfile)
        if not blast_hits:
            sys.stderr.write("Error reading blast output file; exiting now\n")
            exit()
    if args.rsem and verify_path(args.rsem):
        with open(args.rsem, 'r') as rsemfile:
            rsems = read_rsem(rsemfile)
        if not rsems:
            sys.stderr.write("Error reading rsem; exiting now\n")
            exit()
    if args.anno and verify_path(args.anno):
        with open(args.anno, "r") as annofile:
            annos = read_annotations(annofile)
        if not annos:
            sys.stderr.write("Error reading annotations; exiting now\n")
            exit()
    if args.blacklist and verify_path(args.blacklist):
        with open(args.blacklist, "r") as transcript_blacklist_file:
            transcript_blacklist = read_transcript_blacklist(transcript_blacklist_file)
        if not transcript_blacklist:
            sys.stderr.write("Error reading blacklist; exiting now\n")
            exit()
    print("done.\n\n")

    # Perform transcript annotation if annotations were provided
    if annos:
        print("Annotating genes ... ")
        annotate_genes(gff["gene"], annos)
        print("done.\n\n")

    # Convert all that data into Transcript objects
    print("Mapping gff data to transcripts ... ")
    transcript_dict = build_transcript_dictionary(seqs, gff["gene"])
    # transcript_dict maps transcript/sequence ids to Transcript objects
    print("done.\n\n")

    # Do magical fixes on transcripts
    print("Cleaning up transcripts ... ")
    for transcript in transcript_dict.values():
        #transcript.remove_contig_from_gene_id()
        transcript.fix_feature_lengths()
        transcript.create_starts_and_stops()
        transcript.fix_multiple_genes()
        transcript.make_positive()
        transcript.match_cds_and_exon_end()
        transcript.fix_phase()
        transcript.fix_long_gene_names()
    print("done.\n\n")

    # Write RSEM info if provided
    if rsems:
        print("Adding RSEM info to transcripts...")
        with open(rsemerrpath, 'w') as rsem_err:
            for rsem in rsems:
                if rsem.transcript_id not in transcript_dict:
                    rsem_err.write("RSEM transcript "+rsem.transcript_id+" does not exist. Skipping...\n")
                    continue
                transcript = transcript_dict[rsem.transcript_id]
                transcript.rsem = Rsem(rsem.tpm, rsem.fpkm, rsem.isopct)
        print("done.\n\n")

        print("Writing RSEM info...")
        with open(outrsempath, "w") as outrsemfile:
            outrsemfile.write("transcript_id\tnumber_of_CDSs\tcontains_complete_CDS\tTPM\tFPKM\tIsoPct\n")
            for transcript_id, transcript in transcript_dict.items():
                if (
                        (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                        or not transcript.passes_filtering()
                ):
                    continue
                cds_count = 0
                contains_complete_cds = False
                for gene in transcript.genes:
                    for mrna in gene["mrna"]:
                        cds_count += len(mrna["cds"])
                        if "start_codon" in mrna and "stop_codon" in mrna:
                            contains_complete_cds = True
                outrsemfile.write("\t".join([transcript_id, str(cds_count), str(contains_complete_cds), str(transcript.rsem.tpm), str(transcript.rsem.fpkm), str(transcript.rsem.isopct)])+"\n")
        print("done.\n\n")

    # Write .tbl file
    print("Writing .tbl file ... ")
    with open(tblpath, "w") as tblfile:
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            tblfile.write(transcript.to_tbl())
    print("done.\n\n")

    # Write .stats file
    print("Writing .stats file ...")
    with open(statspath, "w") as statsfile:
        statsfile.write("transcript_id\tcomplete\tpfam_domain\tgo_annotation"
                "\tgene_name\n")
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            statsfile.write(transcript.stats())
    print("done.\n\n")

    # Write .gff file
    print("Writing new .gff file...")
    with open(outgffpath, "w") as outgfffile:
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            write_gff(outgfffile, transcript.genes[0])
    print("done.\n\n")

    # Write .fsa file
    print("Writing .fsa file ... ")
    with open(outfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            outfastafile.write(transcript.sequence.to_fasta())
    print("done.\n\n")
    print("Writing .cds.fasta file ... ")
    with open(outcdsfastapath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            mrna = transcript.get_gene().get_mrna()
            cds = mrna.get_cds()
            header = mrna.attributes['ID']
            bases = get_subsequence(transcript.sequence.bases, cds.start, cds.end)
            outfastafile.write(Sequence(header, bases).to_fasta())
    print("done.\n\n")
    print("Writing .cds.pep file ... ")
    with open(outcdspeppath, "w") as outfastafile:
        for transcript in transcript_dict.values():
            if (
                    (transcript_blacklist and transcript.sequence.header in transcript_blacklist)
                    or not transcript.passes_filtering()
            ):
                continue
            mrna = transcript.get_gene().get_mrna()
            cds = mrna.get_cds()
            header = mrna.attributes['ID']
            bases = translate(get_subsequence(transcript.sequence.bases, cds.start + cds.phase, cds.end), cds.strand)
            outfastafile.write(Sequence(header, bases).to_fasta())
    print("done.\n\n")


###########################################

if __name__ == '__main__':
    main()
