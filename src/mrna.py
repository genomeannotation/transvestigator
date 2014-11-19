from src.gff_feature import GFFFeature

class Mrna(GFFFeature):
    def __init__(self, seqid=None, source=None, start=None, end=None, score=None, strand=None, phase=0, attributes=None):
        GFFFeature.__init__(self, seqid, source, "mRNA", start, end, score, strand, phase, attributes)

    def from_gff_feature(feature):
        if feature.type == "mRNA":
            return Mrna(feature.seqid, feature.source, feature.start, feature.end, feature.score, feature.strand, feature.phase, feature.attributes)
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
