import itertools
import sys

from typing import Optional, Iterator, List, Any, Iterable, Generator, Tuple, Callable

from common import params
from common.alignment_storage import AlignmentPiece, AlignedRead
from disjointig_resolve.correction import Correction
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import Segment, Contig

import accurate_line


class LineListener:
    def __init__(self, rc):
        # type: (LineListener) -> None
        self.rc = rc

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (accurate_line.NewLine, Contig, str) -> None
        pass

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (accurate_line.NewLine, Contig, int) -> None
        pass

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        pass

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (accurate_line.NewLine, str, Optional[List[AlignmentPiece]]) -> None
        pass

    def fireAfterCutRight(self, line, pos):
        # type: (accurate_line.NewLine, int) -> None
        pass

    def fireAfterCorrect(self, line, alignments):
        # type: (accurate_line.NewLine, Correction) -> None
        pass

    def fireMerge(self, new_line, line_left, line_right):
        # type: (accurate_line.NewLine, accurate_line.NewLine, accurate_line.NewLine) -> None
        pass

    def fireSplit(self, old_line, line_left, line_right):
        # type: (accurate_line.NewLine, accurate_line.NewLine, accurate_line.NewLine) -> None
        pass


class SmartStorage(LineListener):
    def __init__(self, rc=None):
        # type: (Optional[SmartStorage]) -> None
        self.items = None  # List
        if rc is None:
            rc = SmartStorage(self)
            self.items = []
        LineListener.__init__(self, rc)
        self.sorted = False
        self.key = lambda item: item

    def isCanonical(self):
        return self.items is not None

    def isSorted(self):
        if self.isCanonical():
            return self.sorted
        else:
            return self.rc.isSorted()

    def add(self, seg):
        assert False
        pass

    def addAll(self, items):
        # type: (Iterable) -> SmartStorage
        for item in items:
            self.add(item)
        return self

    def sort(self):
        # type: () -> None
        if self.isCanonical():
            if not self.sorted:
                self.items = sorted(self.items, key=self.key)
                self.sorted = True
        else:
            self.rc.sort()

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeIntLine(len(self))
        for item in self:
            item.save(handler)

    def __len__(self):
        if self.isCanonical():
            return len(self.items)
        else:
            return len(self.rc)

    def clean(self):
        if not self.isCanonical():
            self.rc.clean()
        else:
            self.items = []


# Collection of segments where no segment contains another.
# The reason is that sortings are the same for strait and rc
class SegmentStorage(SmartStorage):
    def __init__(self, rc=None):
        # type: (Optional[SegmentStorage]) -> None
        items = None  # type: List[Segment]
        if rc is None:
            rc = SegmentStorage(self)
            items = []  # type: List[Segment]
        SmartStorage.__init__(self, rc)
        self.items = items
        self.rc = rc # type: SegmentStorage
        self.key = lambda seg: (seg.left, seg.right)

    def __getitem__(self, item):
        # type: (int) -> Segment
        if self.isCanonical():
            self.sort()
            return self.items[item]
        else:
            return self.rc.__getitem__(-1 - item).RC()

    def __iter__(self):
        # type: () -> Iterator[Segment]
        self.sort()
        if self.isCanonical():
            for item in self.items:
                yield item
        else:
            for item in self.rc.items[::-1]:
                yield item.RC()

    def __str__(self):
        if self.isCanonical():
            sign = "+"
        else:
            sign = "-"
        return "ReadStorage" + sign + ":" + str(list(self))

    def __repr__(self):
        return self.__str__()

    def add(self, seg):
        if self.isCanonical():
            self.items.append(seg)
            if self.sorted and (len(self) == 1 or self.key(self.items[-2]) < self.key(self.items[-1])):
                pass
            else:
                self.sorted = False
        else:
            self.rc.add(seg.RC())

    def remove(self, seg):
        # type: (Segment) -> None
        if self.isCanonical():
            l = len(self)
            isin = seg in self.items
            self.items.remove(seg)
            assert not isin or len(self) == l - 1
        else:
            self.rc.remove(seg.RC())

    def addAll(self, segs):
        # type: (Iterable[Segment]) -> SegmentStorage
        for seg in segs:
            self.add(seg)
        return self

    def isIn(self, seg):
        # type: (Segment) -> bool
        if self.isCanonical():
            for seg1 in self.items:
                if seg1.contains(seg):
                    return True
            return False
        else:
            return self.rc.isIn(seg.RC())

    def inter(self, seg, min_size = 0):
        # type: (Segment, int) -> bool
        if self.isCanonical():
            for seg1 in self.items:
                if seg1.interSize(seg) >= min_size:
                    return True
            return False
        else:
            return self.rc.inter(seg.RC(), min_size)

    def interSize(self, seg):
        # type: (Segment) -> bool
        res = 0
        if self.isCanonical():
            for seg1 in self.items:
                res = max(res, seg1.interSize(seg))
            return False
        else:
            return self.rc.interSize(seg.RC())

    # Merge all segments that have intersection of at least inter_size and remove all subsegments
    # After this operation segment ordering is the same when we sort by left bound as if we sort by right bounds
    def mergeSegments(self, inter_size = 0):
        # type: (int) -> None
        if self.isCanonical():
            self.sort()
            if len(self.items) == 0:
                return
            res = [self.items[0]]  # type: List[Segment]
            for seg in self.items[1:]:  # type: Segment
                if seg.left <= res[-1].right - inter_size or seg.right <= res[-1].right:
                    res[-1] = Segment(res[-1].contig, res[-1].left, max(res[-1].right, seg.right))
                else:
                    res.append(seg)
            res = [seg for seg in res if len(seg) >= inter_size]
            self.items = res
        else:
            self.rc.mergeSegments(inter_size)

    def reduce(self, seg):
        # type: (Segment) -> SegmentStorage
        if self.isCanonical():
            res = SegmentStorage()
            for seg1 in self.items:
                if seg.inter(seg1):
                    res.add(seg.cap(seg1))
            return res
        else:
            return self.rc.reduce(seg.RC()).rc

    def filter(self, f):
        # type: (Callable[[Segment], bool]) -> SegmentStorage
        return SegmentStorage().addAll([seg for seg in self if f(seg)])

    def makeCanonical(self):
        if self.isCanonical():
            return
        self.items = [seg.RC() for seg in self.rc.items[::-1]]
        self.rc.items = None
        self.sorted = self.rc.sorted

    def map(self, al):
        # type: (AlignmentPiece) -> SegmentStorage
        matching = al.matchingSequence()
        res = SegmentStorage()
        for seg in self:
             res.add(matching.mapSegDown(al.seg_to.contig, seg))
        return res

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (accurate_line.NewLine, Contig, str) -> None
        self.makeCanonical()


    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (accurate_line.NewLine, Contig, int) -> None
        self.makeCanonical()
        self.sort()
        self.items = [seg for seg in self.items if seg.left < pos] # type: List[Segment]
        if len(self.items) > 0 and self.items[-1].right > pos:
            assert self.items[-1].contig == line
            self.items[-1] = self.items[-1].cap(line.segment(line.left(), pos))

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, correction):
        # type: (Correction) -> None
        self.makeCanonical()
        self.sort()
        self.items = correction.mapSegmentsUp(self.items)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (accurate_line.NewLine, str, Optional[List[AlignmentPiece]]) -> None
        pass

    def fireAfterCutRight(self, line, pos):
        # type: (accurate_line.NewLine, int) -> None
        pass

    def fireAfterCorrect(self, line, alignments):
        # type: (accurate_line.NewLine, Correction) -> None
        self.makeCanonical()
        self.items = [seg.contigAsSegment(line.segment(line.left(), line.right())) for seg in self.items]

    def merge(self, other, inter_size = 0):
        # type: (SegmentStorage, int) -> SegmentStorage
        new_storge = SegmentStorage()
        new_storge.addAll(self)
        new_storge.addAll(other)
        new_storge.mergeSegments(inter_size)
        return new_storge

    def subStorage(self, seg, inter_size = 0):
        # type: (Segment, int) -> SegmentStorage
        if self.isCanonical():
            res = SegmentStorage()
            for seg1 in self:
                if seg1.interSize(seg) > inter_size:
                    res.add(seg1)
            return res
        else:
            return self.rc.subStorage(seg.RC(), inter_size).rc

    # here we assume that segments inside each colection can not have intersection at least min_inter
    def cap(self, other = None, min_inter = 0, seg = None):
        # type: (SegmentStorage, int, Segment) -> SegmentStorage
        if other is None:
            other = SegmentStorage().addAll([seg])
        if not self.isCanonical():
            return self.rc.cap(other.rc, min_inter).rc
        cur = 0
        res = SegmentStorage()
        for seg in other:
            while cur < len(self) and self.items[cur].right <= seg.right:
                if seg.interSize(self.items[cur]) >= min_inter:
                    res.add(seg.cap(self.items[cur]))
                cur += 1
            if cur < len(self) and seg.interSize(self.items[cur]) >= min_inter:
                res.add(seg.cap(self.items[cur]))
        return res

    # This method returns Segment storage that only contain intersections of segments seg1, seg2 when seg1 <= seg2
    def orderedCap(self, other, min_inter = 0):
        # type: (SegmentStorage, int) -> SegmentStorage
        res = SegmentStorage()
        cur = 0
        for seg in list(self):
            while cur < len(other) and other[cur].right < seg.right:
                cur += 1
            if cur < len(other) and seg <= other[cur] and seg.interSize(other[cur]) >= min_inter:
                res.add(seg.cap(other[cur]))
        return res

    def expand(self, radius):
        # type: (int) -> SegmentStorage
        if self.isCanonical():
            return SegmentStorage().addAll([seg.expand(radius) for seg in self.items])
        else:
            return self.rc.expand(radius).rc

    def filterBySize(self, min = 0, max = 10000000000):
        # type: (int, int) -> SegmentStorage
        return self.filter(lambda seg: min <= len(seg) <= max)

    # Maximal collection of segments of size at least inter_size + 1 that have overlap with givel collection of length t most inter_size
    def reverse(self, contig = None, inter_size = 0):
        # type: (Optional[Contig], int) -> SegmentStorage
        assert contig is not None or len(self) != 0
        if contig is None:
            contig = self.items[0].contig
        if len(self) == 0:
            return SegmentStorage().addAll([contig.asSegment()])
        if self.isCanonical():
            self.sort()
            res = SegmentStorage()
            contig = self.items[0].contig
            last = 0
            for seg in self:
                if len(seg) < inter_size:
                    continue
                if seg.left > last:
                    res.add(contig.segment(last, seg.left + inter_size))
                last = max(seg.right - inter_size, last)
            if last + inter_size < len(contig):
                res.add(contig.asSegment().suffix(pos=last))
            return res
        else:
            return self.rc.reverse(contig.rc, inter_size).rc

    def contigAsSegment(self, seg):
        # type: (Segment) -> SegmentStorage
        if self.isCanonical():
            res = SegmentStorage()
            res.addAll(map(lambda seg1: seg1.contigAsSegment(seg), self))
            return res
        else:
            return self.rc.contigAsSegment(seg.RC()).rc

    def load(self, handler, contig):
        # type: (TokenReader, NamedSequence) -> None
        n = handler.readInt()
        for i in range(n):
            self.add(Segment.load(handler, contig))

    #returns the leftmost segment that intersects with seg by at least min_inter (or tockes seg in case min_inter = 0)
    def find(self, seg, min_inter = 0):
        # type: (Segment, int) -> Optional[Segment]
        self.sort()
        for candidate in self:
            if candidate.interSize(seg) >= min_inter:
                return candidate
        return None

    def allInter(self, seg):
        # type: (Segment) -> Generator[Segment]
        for seg1 in self:
            if seg1.inter(seg):
                yield seg1

    def expandTo(self, k):
        if self.isCanonical():
            return SegmentStorage().addAll([seg.expandToSize(k) for seg in self.items])
        else:
            return self.rc.expand(k).rc


class AlignmentStorage(SmartStorage):
    def __init__(self, rc=None):
        # type: (Optional[AlignmentStorage]) -> None
        items = None  # type: List[AlignmentPiece]
        if rc is None:
            rc = AlignmentStorage(self)
            items = []
        SmartStorage.__init__(self, rc)
        self.items = items # type: List[AlignmentPiece]
        self.rc = rc # type: AlignmentStorage
        self.key = lambda al: (al.seg_to.left, al.seg_to.right, al.seg_from.contig.id)
        self.max_len = params.max_read_length

    def checkLine(self, line):
        for al in self:
            if al.seg_to.contig != line:
                return False
        return True

    def __getitem__(self, item):
        # type: (int) -> AlignmentPiece
        if self.isCanonical():
            return self.items[item]
        else:
            return self.rc.__getitem__(-1 - item).rc

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        self.sort()
        if self.isCanonical():
            for item in self.items:
                yield item
        else:
            for item in self.rc.items[::-1]:
                yield item.rc

    def makeCanonical(self):
        if self.isCanonical():
            return
        self.items = [al.rc for al in self.rc.items[::-1]]
        self.rc.items = None
        self.sorted = self.rc.sorted

    def add(self, al):
        if self.isCanonical():
            self.items.append(al)
            if self.sorted and (len(self) == 1 or self.key(self.items[-2]) < self.key(self.items[-1])):
                pass
            else:
                self.sorted = False
            self.max_len = max(self.max_len, len(al.seg_to))
            self.rc.max_len = self.max_len
        else:
            self.rc.add(al.rc)

    def addAndMergeRight(self, al):
        # type: (AlignmentPiece) -> None
        if self.isCanonical():
            for i, al1 in enumerate(self.items): # type: int, AlignmentPiece
                if al.seg_from.inter(al1.seg_from) and al.seg_to.inter(al1.seg_to) and al1.seg_from.right <= al.seg_from.right:
                    tmp = AlignmentPiece.MergeOverlappingAlignments([al1, al])
                    if tmp is not None:
                        self.items[i] = tmp
                        return
            self.add(al)

        else:
            self.rc.addAndMergeLeft(al.rc)

    def addAndMergeLeft(self, al):
        # type: (AlignmentPiece) -> None
        if self.isCanonical():
            for i, al1 in enumerate(self.items): # type: int, AlignmentPiece
                if al.seg_from.inter(al1.seg_from) and al.seg_to.inter(al1.seg_to) and al1.seg_from.left >= al.seg_from.left:
                    tmp = AlignmentPiece.MergeOverlappingAlignments([al, al1])
                    if tmp is not None:
                        self.items[i] = tmp
                        return
            self.add(al)
        else:
            self.rc.addAndMergeRight(al.rc)

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (accurate_line.NewLine, Contig, str) -> None
        self.makeCanonical()
        self.items = [al.targetAsSegment(Segment(new_seq, 0, len(line))) for al in self.items]

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (accurate_line.NewLine, str, Optional[List[AlignmentPiece]]) -> None
        self.makeCanonical()
        self.items = [al.targetAsSegment(line.asSegment()) for al in self.items]

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (accurate_line.NewLine, Contig, int) -> None
        self.makeCanonical()
        new_items = []
        for al in self.items: # type: AlignmentPiece
            if al.seg_to.right <= pos:
                new_items.append(al.changeTargetContig(new_seq))
            elif al.seg_to.left <= pos - params.k:
                new_items.append(al.reduce(target=Segment(line, line.left(), pos)).changeTargetContig(new_seq))
        self.items = new_items # type: List[Segment]

    def fireAfterCutRight(self, line, pos):
        # type: (accurate_line.NewLine, int) -> None
        self.makeCanonical()
        self.items = [al.changeTargetContig(line) for al in self.items]

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, correction):
        # type: (Correction) -> None
        self.makeCanonical()
        self.items = correction.composeQueryDifferences(self.items) # type: List[AlignmentPiece]

    def fireAfterCorrect(self, line, alignments):
        # type: (accurate_line.NewLine, Correction) -> None
        self.makeCanonical()
        assert len(self) == 0 or len(self.items[0].seg_to.contig) == len(line), str([self.items[0].seg_to.contig, line])
        self.items = [al.targetAsSegment(line.asSegment()) for al in self.items] # type: List[AlignmentPiece]

    # Optimize? We are only interested in some of the last alignments.
    def allInter(self, seg, min_inter = 1):
        # type: (Segment, int) -> Generator[AlignmentPiece]
        if self.isCanonical():
            self.sort()
            l = -1 # items[l].left < seg.right
            r = len(self.items) # items[r].left >= seg.right
            while l + 1 < r:
                m = (l + r) / 2
                if self.items[m].seg_to.left >= seg.right:
                    r = m
                else:
                    l = m
            # res = []
            while l >= 0 and self.items[l].seg_to.left >= seg.left - self.max_len:
                if self.items[l].seg_to.interSize(seg) >= min_inter:
                    # res.append(self.items[r])
                    yield self.items[l]
                l -= 1
            # for al in self: # type: AlignmentPiece
            #     if al.seg_to.interSize(seg) >= min_inter:
            #         yield al
        else:
            for al in self.rc.allInter(seg.RC(), min_inter):
                yield al.rc

    def sort(self):
        # type: () -> None
        if self.isCanonical() and self.isSorted():
            return
        self.makeCanonical()
        self.items = sorted(self.items, key=self.key)


    def removeInter(self, seg):
        # type: (Segment) -> None
        if not self.isCanonical():
            self.rc.removeInter(seg.RC())
        else:
            self.items = filter(lambda al: not al.seg_to.inter(seg), self.items)

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        for al in self: # type: AlignmentPiece
            if al.seg_to.contains(seg):
                yield al

    def load(self, handler, collection_from, collection_to):
        # type: (TokenReader, Any, Any) -> None
        n = handler.readInt()
        for i in range(n):
            self.add(AlignmentPiece.load(handler, collection_from, collection_to))

    def reverse(self):
        res = AlignmentStorage()
        for al in self:
            res.add(al.reverse())
        return res

    # This works in square time in worst case bu should work fast if alignments are to left and right sides of the contig
    def merge(self, other):
        # type: (AlignmentStorage) -> AlignmentStorage
        left_items = [(al, -1) for al in self]
        right_items = [(al, 1) for al in other]
        new_items = sorted(left_items + right_items, key = lambda (al, side): (al.seg_to.contig.id, al.seg_from.contig.id, al.seg_from.left))
        new_storge = AlignmentStorage()
        for (c_to, c_from), it in itertools.groupby(new_items, lambda al: (al[0].seg_to.contig, al[0].seg_from.contig)):
            al_sides = list(it)
            als_left = [al for al, side in al_sides if side == -1] # type: List[AlignmentPiece]
            als_left = sorted(als_left, key = lambda al: al.seg_from.right)
            als_right = [al for al, side in al_sides if side == 1] # type: List[AlignmentPiece]
            als_right = sorted(als_right, key = lambda al: al.seg_from.left)
            merged = []
            for al in als_left:
                for j in range(len(als_right)): # type: int
                    if als_right[j] is None:
                        continue
                    if als_right[j].seg_from.left >= al.seg_from.right:
                        break
                    if als_right[j] is not None and al.canMergeTo(als_right[j]):
                        tmp = AlignmentPiece.MergeOverlappingAlignments([al, als_right[j]])
                        if tmp is not None:
                            al = tmp
                            als_right[j] = None
                            break
                merged.append(al)
            new_storge.addAll(merged)
            new_storge.addAll([al for al in als_right if al is not None])
        return new_storge

    def subStorage(self, query = None, target = None, inter_size = 0):
        # type: (Optional[Segment], Optional[Segment], int) -> AlignmentStorage
        if self.isCanonical():
            items = self.items
            if query is not None:
                items = filter(lambda al: al.seg_from.interSize(query) >= inter_size, items)
            if target is not None:
                items = filter(lambda al: al.seg_to.interSize(target) >= inter_size, items)
            self.items[0].reduce(query=query, target=target)
            items = map(lambda al: al.reduce(query=query, target=target), items)
            res = AlignmentStorage()
            res.addAll(items)
            return res
        else:
            return self.rc.subStorage(query, target, inter_size).rc

    def targetAsSegment(self, seg):
        # type: (Segment) -> AlignmentStorage
        if self.isCanonical():
            res = AlignmentStorage()
            res.addAll(map(lambda al: al.targetAsSegment(seg), self))
            return res
        else:
            return self.rc.targetAsSegment(seg.RC()).rc

    def queryAsSegment(self, seg):
        # type: (Segment) -> AlignmentStorage
        if self.isCanonical():
            res = AlignmentStorage()
            res.addAll(map(lambda al: al.queryAsSegment(seg), self))
            return res
        else:
            return self.rc.queryAsSegment(seg.RC())

    def filterByCoverage(self, mi = 0, ma = 1000000, k = 1):
        # type: (int, int, int) -> SegmentStorage
        if len(self) == 0:
            return SegmentStorage()
        segs = list(self.calculateCoverage(k))
        last = None
        contig = self[0].seg_to.contig
        res = SegmentStorage()
        for seg, cov in segs:
            if mi <= cov < ma:
                if last is None:
                    last = seg
                else:
                    last = last.cup(seg)
            elif last is not None and (cov < mi or cov >= ma):
                res.add(last)
                last = None
        if last is not None:
            res.add(last)
        return res

    def calculateCoverage(self, k = 1):
        # type: (int) -> Generator[Tuple[Segment, int]]
        positions = []
        for al in self:
            if len(al.seg_to) >= k:
                positions.append((al.seg_to.left, 1))
                positions.append((al.seg_to.right - k + 1, -1))
        positions = sorted(positions)
        positions = [(pos, sum(delta for pos, delta in iter)) for pos, iter in
                     itertools.groupby(positions, key=lambda pos: pos[0])]
        contig = self[0].seg_to.contig
        cur = 0
        last = 0
        for pos, delta in positions:
            if delta != 0 and last < pos:
                yield contig.segment(last, pos + k - 1), cur
                last = pos
            cur += delta
        if len(contig) - positions[-1][0] >= k:
            yield contig.asSegment().suffix(pos=positions[-1][0]), 0

    def calculateWindowedCoverage(self, sz):
        line = self.__iter__().next().seg_to.contig
        cseg = line.segment(0, sz)
        val= 0
        for seg, c in self.calculateCoverage():
            if seg.inter(cseg):
                val += c * seg.interSize(cseg)
            if seg.right > cseg.left:
                yield cseg, float(val) / sz
                if cseg.right + sz <= len(line):
                    cseg = line.segment(cseg.right, cseg.right + sz)
                    val = c * seg.interSize(cseg)
                else:
                    return
                yield cseg, float(val) / sz


    def __str__(self):
        if self.isCanonical():
            sign = "+"
        else:
            sign = "-"
        return "AlignmentStorage" + sign  + str(list(self))

    def filter(self, condition):
        # type: (Callable[[AlignmentPiece], bool]) -> AlignmentStorage
        return AlignmentStorage().addAll(filter(condition, self))

class ReadAlignmentStorage(AlignmentStorage):
    def __init__(self, rc=None):
        # type: (Optional[ReadAlignmentStorage]) -> None
        items = None  # type: List[AlignmentPiece]
        if rc is None:
            rc = ReadAlignmentStorage(self)
            items = []
        AlignmentStorage.__init__(self, rc)
        self.items = items # type: List[AlignmentPiece]
        self.rc = rc # type: ReadAlignmentStorage

    def add(self, al):
        if self.isCanonical():
            AlignmentStorage.add(self, al)
            read = al.seg_from.contig #type: AlignedRead
            read.addAlignment(al)
        else:
            self.rc.add(al.rc)

    def addAndMergeRight(self, al):
        # type: (AlignmentPiece) -> None
        if self.isCanonical():
            for i, al1 in enumerate(self.items): # type: int, AlignmentPiece
                if al.seg_from.inter(al1.seg_from) and al.seg_to.inter(al1.seg_to) and al1.seg_from.left <= al.seg_from.left:
                    tmp = AlignmentPiece.MergeOverlappingAlignments([al1, al])
                    if tmp is not None:
                        read = al.seg_from.contig #type: AlignedRead
                        read.replaceAlignment(self.items[i], al)
                        self.items[i] = tmp
                        return
            self.add(al)

        else:
            self.rc.addAndMergeLeft(al.rc)

    def addAndMergeLeft(self, al):
        # type: (AlignmentPiece) -> None
        if self.isCanonical():
            for i, al1 in enumerate(self.items): # type: int, AlignmentPiece
                if al.seg_from.inter(al1.seg_from) and al.seg_to.inter(al1.seg_to) and al1.seg_from.left >= al.seg_from.left:
                    tmp = AlignmentPiece.MergeOverlappingAlignments([al, al1])
                    if tmp is not None:
                        read = al.seg_from.contig #type: AlignedRead
                        read.replaceAlignment(self.items[i], al)
                        self.items[i] = tmp
                        return
            self.add(al)
        else:
            self.rc.addAndMergeRight(al.rc)

    def replaceAlignments(self, als):
        self.makeCanonical()
        for al1, al2 in zip(self.items, als):
            read = al1.seg_from.contig # type: AlignedRead
            read.replaceAlignment(al1, al2)
        self.items = filter(lambda al: al is not None, als)

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (accurate_line.NewLine, Contig, str) -> None
        self.makeCanonical()
        self.replaceAlignments([al.targetAsSegment(Segment(new_seq, 0, len(line))) for al in self])

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (accurate_line.NewLine, str, Optional[List[AlignmentPiece]]) -> None
        self.makeCanonical()
        self.replaceAlignments([al.targetAsSegment(line.asSegment()) for al in self.items])

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (accurate_line.NewLine, Contig, int) -> None
        self.makeCanonical()
        new_items = []
        for al in self.items: # type: AlignmentPiece
            if al.seg_to.right <= pos:
                new_items.append(al.changeTargetContig(new_seq))
            elif al.seg_to.left <= pos - params.k:
                new_items.append(al.reduce(target=Segment(line, line.left(), pos)).changeTargetContig(new_seq))
            else:
                new_items.append(None)
        self.replaceAlignments(new_items)

    def fireAfterCutRight(self, line, pos):
        # type: (accurate_line.NewLine, int) -> None
        self.makeCanonical()
        self.replaceAlignments([al.changeTargetContig(line) for al in self.items])

    # alignments from new sequence to new sequence
    def fireBeforeCorrect(self, correction):
        # type: (Correction) -> None
        self.makeCanonical()
        self.replaceAlignments(correction.composeQueryDifferences(self.items))

    def fireAfterCorrect(self, line, alignments):
        # type: (accurate_line.NewLine, Correction) -> None
        self.makeCanonical()
        assert len(self) == 0 or len(self.items[0].seg_to.contig) == len(line), str([self.items[0].seg_to.contig, line])
        self.replaceAlignments([al.targetAsSegment(line.asSegment()) for al in self.items])

    def clean(self):
        if not self.isCanonical():
            self.rc.clean()
        else:
            for al in self.items:
                read = al.seg_from.contig # type: AlignedRead
                read.removeAlignment(al)
            self.items = []



