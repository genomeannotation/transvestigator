from src.gff_feature import GFFFeature
from src.sequtil import get_subsequence, has_start_codon, has_stop_codon, reverse_complement

class Mrna(GFFFeature):
    def __init__(self, seqid=None, source=None, start=None, end=None, score=None, strand=None, phase=0, attributes=None, children=None):
        GFFFeature.__init__(self, seqid, source, "mRNA", start, end, score, strand, phase, attributes, children)

    def from_gff_feature(feature):
        if feature.type == "mRNA":
            return Mrna(feature.seqid, feature.source, feature.start, feature.end, feature.score, feature.strand, feature.phase, feature.attributes, feature.children)
        return None

    def get_cds(self):
        return self["cds"][0]

    def get_exon(self):
        return self["exon"][0]

    def make_positive(self, seq_len):
        if self.strand == '+':
            return
        self.start, self.end = seq_len-self.end+1, seq_len-self.start+1
        self.strand = "+"
        for cds in self["cds"]:
            cds.start, cds.end = seq_len-cds.end+1, seq_len-cds.start+1
            cds.strand = "+"
        for exon in self["exon"]:
            exon.start, exon.end = seq_len-exon.end+1, seq_len-exon.start+1
            exon.strand = "+"

    def match_cds_and_exon_end(self):
        #Check exon/CDS. If no stop codon, make their ends equal.
        #This is a blind attempt to avoid PartialProblem errors from the NCBI.
        if "stop_codon" in self:
            return
        else:
            if self.get_cds().end != self.get_exon().end:
                self.get_cds().end = self.get_exon().end

    def create_starts_and_stops(self, bases):
        cds = self.get_cds()
        subseq = get_subsequence(bases, cds.start, cds.end)
        if cds.strand == '-':
            subseq = reverse_complement(subseq)
        if has_start_codon(subseq):
            seqid = cds.seqid
            source = cds.source
            type = "start_codon"
            codon_start = cds.start
            codon_end = cds.start + 2
            score = None
            strand = cds.strand
            phase = cds.phase
            mrna_id = self.attributes["ID"]
            attributes = {"ID": mrna_id+":start", "Parent": mrna_id}
            start_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                    strand, phase, attributes)
            self.add_child(start_codon)
        if has_stop_codon(subseq):
            seqid = cds.seqid
            source = cds.source
            type = "stop_codon"
            codon_start = cds.end - 2
            codon_end = cds.start
            score = None
            strand = cds.strand
            phase = cds.phase
            mrna_id = self.attributes["ID"]
            attributes = {"ID": mrna_id+":stop", "Parent": mrna_id}
            stop_codon = GFFFeature(seqid, source, type, codon_start, codon_end, score,
                                    strand, phase, attributes)
            self.add_child(stop_codon)
