from typing import Dict, Optional, BinaryIO, Callable, Iterator, Generator, Iterable, Union
from common import basic, SeqIO, params
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, UniqueList, Contig, ContigStorage
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from disjointig_resolve.smart_storage import AlignmentStorage


class Disjointig(Contig):
    def __init__(self, seq, id, rc=None):
        # type: (str, str, Optional[Disjointig]) -> None
        self.seq = seq
        self.id = id
        if rc is None:
            self.read_alignments = AlignmentStorage() # type: AlignmentStorage
            rc = Disjointig(basic.RC(seq), basic.Reverse(id), self) # type: Disjointig
            self.rc = rc
        else:
            self.rc = rc
            self.read_alignments = self.rc.read_alignments.rc  # type: AlignmentStorage
        Contig.__init__(self, seq, id, rc)
        self.rc = rc # type:Disjointig

    def addAlignments(self, als):
        # type: (Iterable[AlignmentPiece]) -> None
        self.read_alignments.addAll(als)

    def addAlignment(self, al):
        # type: (AlignmentPiece) -> None
        self.read_alignments.add(al)

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.read_alignments.getAlignmentsTo(seg)

    def allInter(self, seg, min_inter = 1):
        # type: (Segment, int) -> Generator[AlignmentPiece]
        return self.read_alignments.allInter(seg, min_inter)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)
        self.read_alignments.save(handler)

    def loadDisjointig(self, handler, reads):
        # type: (TokenReader, ReadCollection) -> None
        self.id = handler.readToken()
        self.rc.id = basic.Reverse(self.id)
        seq = handler.readToken()
        self.read_alignments.load(handler, reads, self)


class DisjointigCollection(ContigStorage):
    def __init__(self):
        ContigStorage.__init__(self, [], True)
        self.items = self.items # type: Dict[str, Disjointig]
        self.cnt = 1

    def __iter__(self):
        # type: () -> Iterator[Disjointig]
        return self.items.values().__iter__()

    def __getitem__(self, item):
        # type: (str) -> Disjointig
        return self.items[item]

    def addNew(self, seq, name = None):
        # type: (str, Optional[str]) -> Disjointig
        if name is None:
            name = "D" + str(self.cnt)
            self.cnt += 1
        new_disjointig = Disjointig(seq, name)
        self.add(new_disjointig)
        return new_disjointig

    def loadFromFasta(self, handler, save_names = False, int_ids = False, filter = lambda rec: True):
        # type: (BinaryIO, bool, bool, Callable[[NamedSequence], bool]) -> DisjointigCollection
        recs = list(SeqIO.parse_fasta(handler))
        if save_names:
            for rec in recs:
                assert rec.id not in self.items.keys() and basic.Reverse(rec.id) not in self.items.keys()
        for rec in recs:
            if not filter(rec):
                continue
            if save_names:
                number = basic.parseNegativeNumber(rec.id)
                if number is not None:
                    self.cnt = max(self.cnt, int(abs(number)) + 1)
                if int_ids:
                    self.addNew(rec.seq, str(number))
                else:
                    self.addNew(rec.seq)
            else:
                self.addNew(rec.seq)
        return self

    def addAlignments(self, als):
        # type: (Union[Generator[AlignmentPiece], Iterable[AlignmentPiece]]) -> None
        cnt = 0
        all = 0
        for al in als:
            dt = al.seg_to.contig # type: Disjointig
            all += 1
            if al.seg_from.left < 500 and al.rc.seg_from.left < 500 and len(list(al.splitRef())) == 1 and len(al) > params.k:
                cnt += 1
                dt.addAlignment(al)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        handler.writeTokens(map(lambda line: line.id, UniqueList(self.items.values())))
        for disjointig in UniqueList(self.items.values()):
            handler.writeTokenLine(disjointig.id)
            handler.writeTokenLine(disjointig.seq)
        for disjointig in UniqueList(self.items.values()):
            disjointig.save(handler)

    def load(self, handler, reads):
        # type: (TokenReader, ReadCollection) -> None
        self.cnt = int(handler.readToken())
        keys = list(handler.readTokens())
        for key in keys:
            id = handler.readToken()
            assert id == key
            seq = handler.readToken()
            disjointig = self.addNew(seq, key)
            assert key == disjointig.id, key + " " + disjointig.id
        for key in keys:
            self.items[key].loadDisjointig(handler, reads)



