import sys

from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence

sys.path.append("py")
from common import SeqIO, basic, params
from typing import Generator, Dict, Optional, Union, Iterable, Any, BinaryIO

class Contig(NamedSequence):
    def __init__(self, seq, id, rc = None):
        # type: (str, str, Optional[Contig]) -> None
        NamedSequence.__init__(self, seq, id)
        if rc is None:
            rc = Contig(basic.RC(seq), basic.Reverse(id), self)
        self.rc = rc

    def asSegment(self):
        # type: () -> Segment
        return Segment(self, 0, len(self))

    def segment(self, left, right):
        return Segment(self,left, right)

    def suffix(self, pos):
        # type: (int) -> Segment
        if pos < 0:
            pos = self.__len__() + pos
        if pos < 0:
            pos = 0
        if pos > len(self):
            pos = len(self)
        return Segment(self, pos, self.__len__())

    def prefix(self, len):
        # type: (int) -> Segment
        len = min(len, self.__len__())
        return Segment(self, 0, len)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.id)
        handler.writeTokenLine(self.seq)

    @staticmethod
    def loadContig(handler):
        # type: (TokenReader) -> Contig
        id = handler.readToken()
        seq = handler.readToken()
        return Contig(seq, id)

    def print_fasta(self, handler):
        # type: (BinaryIO) -> None
        SeqIO.write(self, handler, "fasta")

class ContigStorage:
    def __init__(self, contigs_list = None, add_rc = True):
        # type: (Optional[Iterable[Contig]], bool) -> None
        self.items = dict() # type: Dict[str, Contig]
        self.add_rc = add_rc
        if contigs_list is not None:
            for item in contigs_list:
                self.add(item)

    def remove(self, item):
        # type: (Contig) -> None
        assert item.id in self.items
        del self.items[item.id]
        if self.add_rc:
            del self.items[item.rc.id]

    def add(self, item):
        # type: (Contig) -> Contig
        self.items[item.id] = item
        if self.add_rc:
            self.items[item.rc.id] = item.rc
        return item

    def addAll(self, items):
        # type: (Iterable[Contig]) -> ContigStorage
        for item in items:
            self.add(item)
        return self

    def __getitem__(self, item):
        # type: (str) -> Optional[Contig]
        if item in self.items:
            return self.items[item]
        elif not self.add_rc and basic.Reverse(item) in self.items:
            return self.items[basic.Reverse(item)].rc
        return None

    def __iter__(self):
        return self.items.values().__iter__()

    def __len__(self):
        return len(self.items)

    def unique(self):
        if self.add_rc:
            for item in self.items.values():
                if basic.isCanonocal(item.id):
                    yield item
        else:
            for item in UniqueList(self).__iter__():
                yield item

    def uniqueSorted(self):
        lids = sorted([l.id for l in self.unique()])
        for lid in lids:
            yield self.items[lid]

    def __contains__(self, item):
        # type: (Contig) -> bool
        return item.id in self.items

    def containsKey(self, key):
        # type: (str) -> bool
        return key in self.items

    def loadFromFasta(self, handler, num_names=True):
        # type: (BinaryIO, bool) -> ContigCollection
        for rec in SeqIO.parse_fasta(handler):
            if num_names:
                self.add(Contig(rec.seq, str(basic.parseNegativeNumberAndMod(rec.id))))
            else:
                self.add(Contig(rec.seq, rec.id))
        return self

    def loadFromFile(self, fname, num_names=True):
        # type: (str, bool) -> ContigCollection
        for rec in SeqIO.parse_by_name(fname):
            if num_names:
                self.add(Contig(rec.seq, str(basic.parseNegativeNumberAndMod(rec.id))))
            else:
                self.add(Contig(rec.seq, rec.id))
        return self

    def writeToFasta(self, handler):
        for contig in self.unique():
            SeqIO.write(contig, handler, "fasta")


class ContigCollection(ContigStorage):
    def __init__(self, contigs_list=None):
        # type: (Optional[Iterable[Contig]]) -> ContigCollection
        ContigStorage.__init__(self, contigs_list, False)

    def filter(self, condition):
        # type: (callable(Contig)) -> ContigCollection
        res = ContigCollection()
        for contig in self.items.values():
            if condition(contig):
                res.add(contig)
        return res

    def print_names(self, handler):
        # type: (file) -> None
        for contig in self:
            handler.write(str(contig.id) + "\n")

    def print_fasta(self, handler):
        # type: (BinaryIO) -> None
        for contig in self.items.values():
            contig.print_fasta(handler)

    def RC(self):
        return ContigCollection(map(lambda contig: contig.rc, self.items.values()))

    def __len__(self):
        # type: () -> int
        return len(self.items)

    def __contains__(self, item):
        return item.id in self.items

    def containsKey(self, key):
        return key in self.items

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine("ContigCollection")
        handler.writeTokens(self.items.keys())
        handler.writeIntLine(len(list(UniqueList(self.items.values()))))
        for contig in UniqueList(self.items.values()):
            contig.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        message = handler.readToken()
        assert message == "ContigCollection", message
        keys = set(handler.readTokens())
        n = handler.readInt()
        for i in range(n):
            contig = Contig.loadContig(handler)
            self.add(contig)
            if contig.rc.id in keys:
                self.add(contig.rc)


class TmpInfo:
    def __init__(self, l):
        self.misc = l


class Segment:
    def __init__(self, contig, left=None, right=None):
        # type: (Union[Contig, NamedSequence], Optional[int], Optional[int]) -> Segment
        if left is None:
            left = 0
        if right is None:
            right = len(contig)
        assert 0 <= left <= right <= len(contig), str([0, left, right, len(contig), contig.id])
        self.contig = contig
        self.left = left
        self.right = right

    def dist(self, other):
        # type: (Segment) -> int
        assert self.contig == other.contig
        if self.inter(other):
            return 0
        if self.right <= other.left:
            return other.left - self.right
        if self.left >= other.right:
            return self.left - other.right
        return 0

    def extendedDist(self, other):
        # type: (Segment) -> int
        assert self.contig == other.contig
        if self.inter(other):
            return -self.interSize(other)
        if self.right <= other.left:
            return other.left - self.right
        if self.left >= other.right:
            return self.left - other.right
        return 0

    def cap(self, other):
        # type: (Segment) -> Segment
        assert self.inter(other)
        return Segment(self.contig, max(self.left, other.left), min(self.right, other.right))

    def cup(self, other):
        # type: (Segment) -> Segment
        assert self.interSize(other) >= 0
        return Segment(self.contig, min(self.left, other.left), max(self.right, other.right))

    def RC(self):
        # type: () -> Segment
        l = len(self.contig)
        return Segment(self.contig.rc, l - self.right, l - self.left)

    def asNamedSequence(self):
        # type: () -> NamedSequence
        return NamedSequence(self.Seq(), self.contig.id + "[" + str(self.left) + "," + str(self.right) + "]")

    def precedes(self, other, delta=0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right <= other.left + delta

    def __le__(self, other):
        # type: (Segment) -> bool
        return self.contig == other.contig and self.left <= other.left and self.right <= other.right

    def connects(self, other, delta=0):
        # type: (Segment, int) -> bool
        return self.contig == other.contig and self.right - delta <= other.left <= self.right + delta

    # Return true if at least one character in common
    def inter(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and not (self.right <= other.left or self.left >= other.right)

    # Return -1 if no intersection, 0 if segments touch, and overlap length if segments overlap
    def interSize(self, other):
        # type: (Segment) -> int
        if self.contig.id != other.contig.id or self.right < other.left or self.left > other.right:
            return -1
        return min(self.right, other.right) - max(self.left, other.left)

    def contains(self, other):
        # type: (Segment) -> bool
        return self.contig.id == other.contig.id and self.left <= other.left and self.right >= other.right

    def asContig(self):
        # type: () -> Contig
        return Contig(self.Seq(),
                      "(" + self.contig.id + ")[" + str(self.left) + "," + str(self.right) + "]")

    def merge(self, other):
        return Segment(self.contig, min(self.left, other.left), max(self.right, other.right))

    def shift(self, val):
        return Segment(self.contig, self.left + val, self.right + val)

    def Seq(self):
        # type: () -> str
        return self.contig.seq[self.left:self.right]

    def __str__(self):
        # type: () -> str
        if (self.right < len(self.contig) * 0.6 or (
                self.right < 0.9 * len(self.contig) and len(self.contig) - self.right > 500)):
            return str(self.contig.id) + "[" + str(self.left) + ":" + str(self.right) + "]"
        else:
            return str(self.contig.id) + "[" + str(self.left) + ":" + str(len(self.contig)) + "-" + str(
                len(self.contig) - self.right) + "]"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.right - self.left

    def __eq__(self, other):
        return self.contig == other.contig and self.left == other.left and self.right == other.right

    def __ne__(self, other):
        return self.contig != other.contig or self.left != other.left or self.right != other.right

    def __hash__(self):
        return hash((self.contig.id, self.left, self.right))

    def changeContig(self, contig):
        # type: (Segment) -> Segment
        return Segment(contig, self.left, self.right)

    def expand(self, range):
        # type: (int) -> Segment
        return Segment(self.contig, max(self.left - range, 0), min(self.right + range, len(self.contig)))

    def expandToSize(self, sz):
        # type: (int) -> Segment
        if len(self) >= sz:
            return self
        d = min(sz, len(self.contig)) - len(self)
        right = min(max((d + 1) / 2, d - self.left), len(self.contig) - self.right)
        left = min(max((d + 1) / 2, d - (len(self.contig) - self.right)), self.left)
        return Segment(self.contig, self.left - left, self.right + right)

    def shrink(self, range):
        # type: (int) -> Segment
        return Segment(self.contig, max(self.left + range, 0), min(self.right - range, len(self.contig)))

    def expandLeft(self, range):
        # type: (int) -> Segment
        return Segment(self.contig, max(self.left - range, 0), self.right)

    def expandRight(self, range):
        # type: (int) -> Segment
        return Segment(self.contig, self.left, min(self.right + range, len(self.contig)))

    def copy(self):
        # type: () -> Segment
        return Segment(self.contig, self.left, self.right)

    def contigAsSegment(self, seg):
        # type: (Segment) -> Segment
        return Segment(seg.contig, self.left + seg.left, self.right + seg.left)

    def suffix(self, pos = None, length = None):
        assert (pos is not None) != (length is not None)
        if length is not None:
            pos = self.right - length
        assert self.left <= pos < self.right
        return Segment(self.contig, pos, self.right)

    def prefix(self, pos = None, length = None):
        assert (pos is not None) != (length is not None)
        if length is not None:
            pos = self.left + length
        assert self.left < pos <= self.right
        return Segment(self.contig, self.left, pos)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.contig.id)
        handler.writeToken(str(self.left))
        handler.writeToken(str(self.right))

    @staticmethod
    def load(handler, contig_or_collection):
        # type: (TokenReader, Any) -> Segment
        id = handler.readToken()
        if isinstance(contig_or_collection, NamedSequence) and contig_or_collection.id == id:
            contig = contig_or_collection
        else:
            contig = contig_or_collection[id]
        assert contig is not None and contig.id == id, id
        return Segment(contig, int(handler.readToken()), int(handler.readToken()))


class Consensus:
    def __init__(self, seq, cov):
        # type: (str, list[int]) -> Consensus
        self.seq = seq
        self.cov = cov

    def suffix(self, pos):
        if pos < 0:
            pos = self.__len__() + pos
        return Consensus(self.seq[pos:], self.cov[pos:])

    def printQuality(self, handler, cov_threshold=params.reliable_coverage):
        # type: (file, int) -> None
        for c, a in zip(self.seq, self.cov):
            if a < cov_threshold:
                handler.write(c.lower())
            else:
                handler.write(c.upper())
        handler.write("\n")

    def printCoverage(self, handler, step):
        # type: (file, int) -> None
        for i, val in list(enumerate(self.cov))[::step]:
            handler.write(str((i, val)) + "\n")

    def cut(self, cov_threshold=params.reliable_coverage, length=None):
        l = len(self.seq)
        while l > 0 and self.cov[l - 1] < cov_threshold:
            l -= 1
        if length is not None:
            l = min(l, length)
        return Consensus(self.seq[:l], self.cov[:l])

    def __len__(self):
        return len(self.seq)

    def RC(self):
        return Consensus(basic.RC(self.seq), self.cov[::-1])

    def merge(self, other, pos=None):
        # type: (Consensus, Optional[int]) -> Consensus
        if pos is None:
            pos = len(self)
        return Consensus(self.seq[:pos] + other.seq, self.cov[:pos] + other.cov)


def UniqueList(sequences):
    # type: (Iterable[NamedSequence]) -> Generator[Any]
    visited = set()
    for seq in sequences:
        if seq.id in visited or basic.Reverse(seq.id) in visited:
            pass
        else:
            yield seq
            visited.add(seq.id)
