import random
from string import ascii_lowercase, ascii_uppercase

from typing import List, Tuple

from alignment.align_tools import Aligner
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import ContigStorage, Contig
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot


class TestDataset:
    def __init__(self, genome = "", letter_size = 550, error_rate = 0.05, mutation_rate = 0.005, seed = 0):
        random.seed(seed)
        self.reads = [] # type: List[NamedSequence]
        self.disjointigs = [] # type: List[NamedSequence]
        self.contigs = [] # type: List[NamedSequence]
        self.letter_size = letter_size
        self.error_rate = error_rate
        self.mutation_rate = mutation_rate
        self.alphabet = ContigStorage()
        self.matches = dict()
        for c1, c2 in zip(ascii_lowercase, ascii_uppercase):
            seq = self.generate(self.letter_size)
            self.alphabet.add(Contig(seq, c1))
            seq, matches = self.mutate(seq, self.mutation_rate)
            self.alphabet.add(Contig(seq, c2))
            self.matches[c1] = matches
            self.matches[c2] = [(b, a) for a, b in matches]
        self.genome = Contig(self.translate(genome), genome)

    def translate(self, seq):
        return "".join(map(lambda c: self.alphabet[c].seq, seq))

    def addRead(self, read_seq):
        name = "R" + str(len(self.reads)) + "_" + read_seq
        self.reads.append(NamedSequence(self.mutate(self.translate(read_seq), self.error_rate)[0], name))
        return name

    def addDisjointig(self, disjointig_seq):
        # type: (str) -> str
        self.disjointigs.append(NamedSequence(self.mutate(self.translate(disjointig_seq), self.mutation_rate)[0], "D" + str(len(self.disjointigs)) + "_" + disjointig_seq))
        return self.disjointigs[-1].id

    def addContig(self, contig_seq):
        # type: (str) -> str
        name = "C" + str(len(self.contigs)) + "_" + contig_seq
        self.contigs.append(NamedSequence(self.translate(contig_seq), name))
        return name

    def generateReads(self, length = 5, cov = 15, circular = False):
        genome = self.genome.id
        if circular:
            genome = genome + genome[0:length - 1]
        for i in range(0, len(genome) - length + 1):
            for j in range((cov + length - 1) / length):
                self.addRead(genome[i:i + length])

    def generate(self, letter_size):
        # type: (int) -> str
        return "".join([random.choice(["A", "C", "G", "T"]) for i in range(letter_size)])

    def genAll(self, aligner):
        # type: (Aligner) -> Tuple[NewLineStorage, LineDotPlot, ReadCollection]
        disjointigs = DisjointigCollection()
        for dis in self.disjointigs:
            disjointigs.addNew(dis.seq, dis.id)
        from disjointig_resolve.line_storage import NewLineStorage
        lines = NewLineStorage(disjointigs, aligner)
        lines.name_printer = lambda line: line.id + "_" + self.translateBack(line, aligner)
        for line in self.contigs:
            new_line = lines.addNew(line.seq, line.id)
            new_line.initial.add(AlignmentPiece.Identical(new_line.asSegment().asContig().asSegment(), new_line.asSegment()))
        dp = LineDotPlot(lines, aligner)
        dp.construct(aligner)
        lines.alignDisjointigs()
        reads = ReadCollection()
        for read in self.reads:
            reads.addNewRead(read)
        disjointigs.addAlignments(aligner.localAlign(reads, disjointigs))
        return lines, dp, reads

    def mutate(self, seq, rate):
        # type: (str, float) -> Tuple[str, List[Tuple[int, int]]]
        res = [seq[0]]
        matches = []
        matches.append((0,0))
        cur = 1
        for i, c in enumerate(seq):
            if i == 0 or i == len(seq) - 1:
                continue
            if random.random() < rate:
                vars = ["A", "C", "G", "T"]
                vars.remove(c)
                res.append(random.choice([random.choice(vars), "", c + c]))
                cur += len(res[-1])
            else:
                res.append(c)
                matches.append((cur, i))
                cur += 1
        res.append(seq[-1])
        matches.append((len(seq) - 1, cur))
        return "".join(res), matches

    def saveStructure(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.genome.id)
        handler.writeInt(len(self.reads))
        for read in self.reads:
            handler.writeToken(read.id.split("_")[-1])
        handler.writeInt(len(self.disjointigs))
        for disjointig in self.disjointigs:
            handler.writeToken(disjointig.id.split("_")[-1])
        handler.writeInt(len(self.contigs))
        for contig in self.contigs:
            handler.writeToken(contig.id.split("_")[-1])

    @staticmethod
    def loadStructure(handler):
        # type: (TokenReader) -> TestDataset
        random.seed(0)
        res = TestDataset(handler.readToken())
        for i in range(handler.readInt()):
            res.addRead(handler.readToken())
        for i in range(handler.readInt()):
            res.addDisjointig(handler.readToken())
        for i in range(handler.readInt()):
            res.addContig(handler.readToken())
        return res

    def translateBack(self, contig, aligner):
        # type: (Contig, Aligner) -> str
        res = []
        for al in sorted(aligner.overlapAlign([contig], self.alphabet), key = lambda al: al.seg_from.left):
            if len(res) > 0 and al.seg_from.interSize(res[-1].seg_from) > self.letter_size / 2:
                if al.percentIdentity() > res[-1].percentIdentity():
                    res[-1] = al
            else:
                res.append(al)
        return "".join([al.seg_to.contig.id for al in res])