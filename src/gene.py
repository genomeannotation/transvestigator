from src.gff_feature import GFFFeature

class Gene(GFFFeature):
    def __init__(self, seqid=None, source=None, start=None, end=None, score=None, strand=None, phase=0, attributes=None, children=None):
        GFFFeature.__init__(self, seqid, source, "gene", start, end, score, strand, phase, attributes, children)

    def from_gff_feature(feature):
        if feature.type == "gene":
            return Gene(feature.seqid, feature.source, feature.start, feature.end, feature.score, feature.strand, feature.phase, feature.attributes, feature.children)
        return None

    def get_mrna(self):
        return self["mrna"][0]

    def get_cds_length(self):
        return self.get_mrna().get_cds().length()

    def is_complete(self):
        has_start, has_stop = False, False
        if "start_codon" in self.get_mrna():
            has_start = True
        if "stop_codon" in self.get_mrna():
            has_stop = True
        if has_start and has_stop:
            return True
        else:
            return False

    def to_tbl(self):
        # Check for starts and stops
        has_start = False
        has_stop = False
        if "start_codon" in self.get_mrna():
            has_start = True
        if "stop_codon" in self.get_mrna():
            has_stop = True
        # Create tbl entry
        tbl = ""
        if not has_start:
            tbl += "<"
        tbl += str(self.start)+"\t"
        if not has_stop:
            tbl += ">"
        tbl += str(self.end)+"\tgene\n"
        # Gene name if it has one
        if "Name" in self.attributes:
            tbl += "\t\t\tgene\t"+self.attributes["Name"]+"\n"
        # Locus tag
        locus_tag = self.attributes["ID"]
        if "|" in locus_tag:
            # Get rid of ugly locus tags that look like c14595_g1_i1|g.5835;
            # the NCBI *hates* them
            fields = locus_tag.split("|")
            if len(fields) > 1:
                locus_tag = fields[1]
        tbl += "\t\t\tlocus_tag\t"+locus_tag+"\n"
        if not has_start:
            tbl += "<"
        tbl += str(self.get_mrna().get_cds().start)+"\t"
        if not has_stop:
            tbl += ">"
        tbl += str(self.get_mrna().get_cds().end)+"\tCDS\n"
        # Codon start
        if self.get_mrna().get_cds().phase != 0:
            tbl += "\t\t\tcodon_start\t"+str(self.get_mrna().get_cds().phase+1)+"\n"
        # Protein id
        tbl += "\t\t\tprotein_id\t"+self.get_mrna().attributes["ID"]+"\n"
        # Dbxref if it has any
        if "Dbxref" in self.get_mrna().attributes:
            for dbxref in self.get_mrna().attributes["Dbxref"].split(","):
                tbl += "\t\t\tdb_xref\t"+dbxref+"\n"
        # product if it has any
        if "product" in self.get_mrna().attributes:
            tbl += "\t\t\tproduct\t"+self.get_mrna().attributes["product"]+"\n"
        else: # no product, write hypothetical protein
            tbl += "\t\t\tproduct\thypothetical protein\n"
        # Ontology_term if it has any
        if "Ontology_term" in self.get_mrna().attributes:
            for term in self.get_mrna().attributes["Ontology_term"].split(","):
                tbl += "\t\t\tOntology_term\t"+term+"\n"
        return tbl

    def remove_contig_from_gene_id(self):
        id_split = self.attributes['ID'].split('|')
        if len(id_split) == 2:
            self.attributes['ID'] = id_split[1]

    def fix_phase(self, bases):
        #Changes start indices and phase values for CDSs starting at 2 or 3.

        #Adjusts start index for partial gene, mRNA and CDS to 1 and adds
        #appropriate phase value for the CDS; we theorize that this
        #is necessary to eliminate errors from the NCBI TSA submission.
        
        gene_start = self.start
        mrna_start = self.get_mrna().start
        cds_start = self.get_mrna().get_cds().start

        # Adjust phase if our feature start on base 2 or 3
        if not "start_codon" in self.get_mrna():
            if gene_start == 2:
                self.start = 1
                self.get_mrna().start = 1
                self.get_mrna().get_cds().start = 1
                self.get_mrna().get_cds().phase = 1
            elif gene_start == 3:
                self.start = 1
                self.get_mrna().start = 1
                self.get_mrna().get_cds().start = 1
                self.get_mrna().get_cds().phase = 2
            if self.get_mrna().get_cds().start == 2:
                self.get_mrna().get_cds().start = 1
                self.get_mrna().get_cds().phase = 1
            elif self.get_mrna().get_cds().start == 3:
                self.get_mrna().get_cds().start = 1
                self.get_mrna().get_cds().phase = 2
        # Adjust end if partial
        if not "stop_codon" in self.get_mrna():
            self.end = len(bases)
            self.get_mrna().end = len(bases)
            self.get_mrna().get_cds().end = len(bases)

    def make_positive(self, seq_len):
        if self.strand == "+":
            return
        self.start, self.end = seq_len-self.end+1, seq_len-self.start+1
        self.strand = "+"
        for mrna in self["mrna"]:
            mrna.make_positive(seq_len)

    def match_cds_and_exon_end(self):
        for mrna in self["mrna"]:
            mrna.match_cds_and_exon_end()

    def create_starts_and_stops(self, bases):
        for mrna in self["mrna"]:
            mrna.create_starts_and_stops(bases)
