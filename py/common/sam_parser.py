import sys
import itertools
from typing import Generator


def CIGAR_to_List(cigar):
    delims = ["M", "I", "D", "N", "S", "H", "P", "=", "X"]
    cigar_list = list()
    num_list = list()
    cur_num = ""
    for s in cigar:
        if s in delims:
            cigar_list.append(s)
            num_list.append(int(cur_num))
            cur_num = ""
        else:
            cur_num += s
    return [cigar_list, num_list]

def UpdateAlignmentLength(align_len, cigar_char, cigar_len, seq_len):
    if cigar_char == "M" or cigar_char == "D" or cigar_char == "N":
        align_len += cigar_len
    elif cigar_char == "=":
        align_len = seq_len
    return align_len

############################# SAM Config class #############################

class SAM_Config:
    # prefix strings
    sq_prefix = "@SQ"
    hd_prefix = "@HD"
    rg_prefix = "@RG"
    pg_prefix = "@PG"
    oc_prefix = "@CO"

    # sq fields
    sq_tname_index = 1
    sq_tname_prefix = "SN:"

    # alignment fields
    num_mand_fields = 11

    # alignment indices
    query_index = 0
    flag_index = 1
    target_index = 2
    pos_index = 3
    mapq_index = 4
    cigar_index = 5
    mate_target_index = 6
    mate_pos_index = 7
    tlen_index = 8
    seq_index = 9
    qual_index = 10

############################# SAM Entry class #############################

class SAM_entry:
    query_name = ""         # string
    flag = 0                # int
    target_name = ""        # string
    pos = 0                 # int
    mapping_qiality = 0     # int
    cigar = ""              # string
    mate_target_name = ""   # string
    mate_pos = 0            # int
    tlen = 0                # int
    seq = ""                # string
    qual = ""               # string

    alen = 0
    sam_config = SAM_Config()

    def ComputeAlignmentLength(self):
        lists = CIGAR_to_List(self.cigar)
        char_list = lists[0]
        lens_list = lists[1]
        for i in range(0, len(char_list)):
            self.alen = UpdateAlignmentLength(self.alen, char_list[i], lens_list[i], len(self.seq))

    def __init__(self, alignment_string):
        splits = alignment_string.split()

        assert len(splits) >= self.sam_config.num_mand_fields, "ERROR: Mandatory fields of alignment were not specified"

        self.query_name = splits[self.sam_config.query_index]
        self.flag = int(splits[self.sam_config.flag_index])
        self.target_name = splits[self.sam_config.target_index]
        self.pos = int(splits[self.sam_config.pos_index])
        self.mapping_quality = int(splits[self.sam_config.mapq_index])
        self.cigar = splits[self.sam_config.cigar_index]
        self.mate_target_name = splits[self.sam_config.mate_target_index]
        self.mate_pos = int(splits[self.sam_config.mate_pos_index])
        self.tlen = int(splits[self.sam_config.tlen_index])
        self.seq = splits[self.sam_config.seq_index]
        self.qual = splits[self.sam_config.qual_index]

        self.ComputeAlignmentLength()

        #if self.cigar != "101M":
        #    self.Print()

############################# SAM Parser class #############################

class SAMEntryInfo:
    def __init__(self, tid, tname, pos, alen, seq, flag, name, qual, cigar, tlen):
        self.tid = tid
        self.tname = tname
        self.pos = pos
        self.alen = alen
        self.seq = seq
        self.is_unmapped = (((flag >> 2) & 1) == 1)
        self.query_name = name
        self.proper_alignment = (((flag >> 1) & 1) == 1)
        self.secondary = (((flag >> 8) & 1) == 1)
        self.rc = (((flag >> 4) & 1) == 1)
        self.flag = flag
        self.qual = qual
        self.cigar = cigar
        self.tlen = tlen


    def Print(self):
        sys.stdout.write(self.query_name + " " + str(self.tid) + " " + str(self.pos) + " " + str(self.alen) + " " + str(self.is_unmapped) + " " + str(self.proper_alignment) + " " + str(self.flag) + " " + str(self.secondary) + "\n")


class Samfile:
    headers = list()        # is not used
    queries = list()        # is not used
    targets = list()        # lines corresponding to references. Can be parsed
    programs = list()       # is not used
    comments = list()       # something strange
    entries = list()        # list of SAM_entry objects
    target_map = dict()     # map "target" -> index

    # auxiliaries
    sam_config = SAM_Config()

    def IsLineReferenceDescr(self, line):
        return line.startswith(self.sam_config.sq_prefix)

    def IsLineHeaderDescr(self, line):
        return line.startswith(self.sam_config.hd_prefix)

    def IsLineReadGroupDescr(self, line):
        return line.startswith(self.sam_config.rg_prefix)

    def IsLineProgramDescr(self, line):
        return line.startswith(self.sam_config.pg_prefix)

    def IsLineComment(self, line):
        return line.startswith(self.sam_config.oc_prefix)

    def GetSAMEntry(self, line):
        sam_entry = SAM_entry(line)
        return sam_entry

    def UpdateTargetFields(self, line):
        # add into array
        self.targets.append(line)

        # add into map
        splits = line.split()
        target_name = splits[self.sam_config.sq_tname_index]
        target_name = target_name[len(self.sam_config.sq_tname_prefix):]
        self.target_map[target_name] = len(self.targets) - 1


    def InitFields(self):
        self.targets = list()
        self.headers = list()
        self.queries = list()
        self.programs = list()
        self.comments = list()
        self.target_map = dict()

    def PrintStats(self):
        sys.stdout.write("# targets:\t" + str(len(self.targets)) + "\n")
        sys.stdout.write("# headers:\t" + str(len(self.headers)) + "\n")
        sys.stdout.write("# queries:\t" + str(len(self.queries)) + "\n")
        sys.stdout.write("# programs:\t" + str(len(self.programs)) + "\n")
        sys.stdout.write("# comments:\t" + str(len(self.comments)) + "\n")
        sys.stdout.write("# entries:\t" + str(len(self.entries)) + "\n")

    def __init__(self, handler):
        self.InitFields()
        self.handler = handler
        self.target_map["*"] = -1
        for line in self.handler:
            if not self.processLine(line.strip()):
                self.handler = itertools.chain([line], self.handler)
                break

    def processLine(self, line):
        line = line.strip()
        # line is reference sequence dictionary
        if self.IsLineReferenceDescr(line):
            self.UpdateTargetFields(line)
            return True
        # line is header
        elif self.IsLineHeaderDescr(line):
            self.headers.append(line)
            return True
        # line is read group
        elif self.IsLineReadGroupDescr(line):
            self.queries.append(line)
            return True
        # line is program
        elif self.IsLineProgramDescr(line):
            self.programs.append(line)
            return True
        # line is comment
        elif self.IsLineComment(line):
            self.comments.append(line)
            return True
        return False

    def __iter__(self):
        # type: () -> Generator[SAMEntryInfo]
        for line in self.handler:
            line = line.strip()
            # line is reference sequence dictionary
            if not self.processLine(line):
                entry = self.GetSAMEntry(line)
                if entry.target_name in self.target_map:
                    tid = self.target_map[entry.target_name]
                else:
                    tid = None
                yield SAMEntryInfo(tid, entry.target_name, entry.pos, entry.alen, entry.seq, entry.flag, entry.query_name, entry.qual, entry.cigar, entry.tlen)


    def NumEntries(self):
        return len(self.entries)

    # iterators
    # def __iter__(self):
    #     return SamIter(self)

    def gettid(self, tname):
        return self.target_map[tname]

def chain_iter(iterators):
    for it in iterators:
        for element in it:
            yield element

class SamChain:
    def __init__(self, sam_files):
        self.sam_files = sam_files

    def __iter__(self):
        return chain_iter(self.sam_files)

    def gettid(self, tname):
        for sam in self.sam_files:
            if sam.gettid(tname) != None:
                return sam.gettid(tname)
        return None
