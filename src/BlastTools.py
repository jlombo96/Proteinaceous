



def parse_fmt7(blastfile:str=None):
    """
    Parse a blastfmt 7 file into a list of dictionaries
    """
    column_name_map = { "query_id":"qseqid",
                        "subject_id":"sseqid",
                        "q._start":"qstart",
                        "q._end":"qend",
                        "s._start":"sstart",
                        "s._end":"send",
                        "evalue":"evalue",
                        "bit_score":"bitscore",
                        "alignment_length":"alen",
                        "mismatches":"mismatches",
                        "gap_opens":"gap_opens",
                        "%_query_coverage_per_hsp":"qcovhsp",
                        "%_query_coverage_per_subject":"qcovs",
                        "%_subject_coverage_per_query":"scovq",
                        "query_length":"qlen",
                        "subject_length":"slen",
                        "%_positives":"posident",
                        "%_identity":"pident"}
    if not blastfile:
        print("ERROR: No blastfile provided")
        return None
    if not isinstance(blastfile,str):
        print("ERROR: blastfile must be a string")
        return None
    blast_hit_list = []
    with open(blastfile,'r',encoding='utf-8') as blaststream:
        for line in blaststream:
            line = line.strip()
            if line.startswith("#"):
                if "BLAST" in line:
                    continue
                elif "Query" in line:
                    continue
                elif "Database" in line:
                    continue
                elif "hits" in line:
                    continue
                elif "Fields" in line:
                    field_lines = _,field_text = line.split(":")
                    field_names = field_text.split(",")
                    new_names = []
                    for name in field_names:
                        name = name.strip().replace(" ","_")
                        new_names.append(name)
                    nfields = len(new_names)
                else:
                    print(line)
            else:
                fields = line.split("\t")
                if len(fields) != nfields:
                    print(f"WARNING {fields}")
                blast_hsp_dict = {}
                for key,value in zip(new_names,fields):
                    if key in column_name_map:
                        blast_hsp_dict[column_name_map[key]] = value
                    else:
                        blast_hsp_dict[key] = value
                blast_hit_list.append(blast_hsp_dict)
    blaststream.close()
    return blast_hit_list



def create_blast_obj(infile:str=None,outfmt='7'):
    """
    """
    blast_hit_list = []
    if outfmt == '7':
        blast_hit_list = parse_fmt7(infile)
    else:
        raise ValueError(f"ERROR: {outfmt} not a valid outfmt")

class BlastResults:
    """
    """
    def __init__(self,alignments:list=None):
        self._qmap = {}
        self._rmap = {}
        self.results = []
        if alignments is None:
            self._alignments = []
        else:
            self._alignments = alignments
        

class BlastAlignment:
    """
    """
    def __init__(self,hsp_list:list=None):
        if hsp_list is None:
            self._hsp_list = []
        else:
            self._hsp_list = hsp_list
    def sort(self,key:str):
        try:
            new_order = self._hsp_list.sort(key=lambda x: getattr(x,key))
        except AttributeError:
            print(f"ERROR: {key} not a valid key")
        self._hsp_list = new_order
    def append(self,hsp):
        self._hsp_list.append(hsp)
    def extend(self,hsp_list):
        self._hsp_list.extend(hsp_list)
    def __len__(self):
        return len(self.hsp_list)
    def flip(self):
        new_list = []
        for hsp in self._hsp_list:
            new_hsp = BlastHsp(hsp)
            qstart,qend,sstart,send = new_hsp.qstart,new_hsp.qend,new_hsp.sstart,new_hsp.send
            qseqid = new_hsp.qseqid
            sseqid = new_hsp.sseqid
            new_hsp.qstart = sstart
            new_hsp.qend = send
            new_hsp.sstart = qstart
            new_hsp.send = qend
            new_hsp.qseqid = sseqid
            new_hsp.sseqid = qseqid
            new_list.append(new_hsp)
        self._hsp_list = new_list

class BlastHsp:
    """
    """
    def __init__(self,hsp_dict:dict):
        required_keys = {"evalue":float,"bitscore":float,"qstart":int,"qend":int,"sstart":int,"send":int,"qseqid":str,"sseqid":str}
        for key,type_to_assign in required_keys.items():
            if key not in hsp_dict:
                print(f"ERROR: {key} not in hsp_dict")
                return None
            value = hsp_dict[key]
            if type_to_assign == float:
                setattr(self,key,float(value))
            elif type_to_assign == int:
                setattr(self,key,int(value))
            elif type_to_assign == str:
                setattr(self,key,str(value))
            else:
                raise ValueError(f"ERROR: {type_to_assign} not a valid type")
        for k,v in hsp_dict.items():
            if k not in required_keys:
                # Test if the variable type is numeric
                if v.isnumeric():
                    setattr(self,k,float(v))
                else:
                    setattr(self,k,str(v))
    def __str__(self):
        return f"{self.qseqid} {self.sseqid} {self.evalue} {self.bitscore} {self.qstart} {self.qend} {self.sstart} {self.send}"
    def __repr__(self):
        return f"{self.qseqid} {self.sseqid} {self.evalue} {self.bitscore} {self.qstart} {self.qend} {self.sstart} {self.send}"
    def __len__(self):
        return self.qend - self.qstart + 1
    def __getitem__(self,key):
        return getattr(self,key)


def main():
    blast_test_file = "/mnt/c/Users/jonlo/Research_Burton/jons_software/Proteinaceous/test/proteinsearch_GCF_012933605.1_protein.blastp"

    blast_hit_list = parse_fmt7(blast_test_file)
    blast_map = {}
    for hit in blast_hit_list:
        print(hit)
        hsp = BlastHsp(hit)
        blast_key = (hsp.qseqid,hsp.sseqid)
        if blast_key not in blast_map:
            blast_map[blast_key] = BlastAlignment()
        blast_map[blast_key].append(hsp)
    for key,value in blast_map.items():
        print(key)
        for hsp in value._hsp_list:
            print(hsp)
        print("\n\n")
if __name__ == "__main__":
    main()