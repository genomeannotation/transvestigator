from src.gff_feature import GFFFeature

class Gene(GFFFeature):
    def __init__(self, seqid=None, source=None, start=None, end=None, score=None, strand=None, phase=0, attributes=None):
        GFFFeature.__init__(self, seqid, source, "gene", start, end, score, strand, phase, attributes)

    def from_gff_feature(feature):
        if feature.type == "gene":
            return Gene(feature.seqid, feature.source, feature.start, feature.end, feature.score, feature.strand, feature.phase, feature.attributes)
        return None
