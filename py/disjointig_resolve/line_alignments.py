from typing import Optional, Generator, Any, List

from common.alignment_storage import AlignmentPiece
from common.save_load import TokenWriter, TokenReader
from common.sequences import Contig, Segment
from disjointig_resolve.correction import Correction
from disjointig_resolve.smart_storage import LineListener, AlignmentStorage


# Unlike all other storages AutoAlignmentStorage stores only half of the alignments. The other half are the reversals of the stored half.
# Also identical alignment can not be stored here but is returned as a part of iteration (see __iter__)
class AutoAlignmentStorage(LineListener):

    def __init__(self, line, rc = None):
        # type: (Contig, Optional[AutoAlignmentStorage]) -> None
        self.line = line
        if rc is None:
            self.content = AlignmentStorage()
            rc = AutoAlignmentStorage(line.rc, self)
            rc.content = self.content.rc
            rc.state = -1
        LineListener.__init__(self, rc)
        self.rc = rc # type: AutoAlignmentStorage
        self.state = 1 # from precedes to

    def makeCanonical(self, al):
        if (self.state == 1) == (al.seg_from.left < al.seg_to.left):
            return al
        else:
            return al.reverse()

    def isCanonical(self, al):
        return (self.state == 1) == (al.seg_from.left < al.seg_to.left)

    def add(self, al):
        # type: (AlignmentPiece) -> None
        if al.isIdentical():
            return
        self.content.add(self.makeCanonical(al))

    def addAll(self, als):
        for al in als:
            self.add(al)
        return self

    def addAndMergeRight(self, al):
        if al.isIdentical():
            return
        if self.isCanonical(al):
            self.content.addAndMergeRight(al)
        else:
            self.content.addAndMergeRight(al.reverse())

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        for al in self.content:
            yield al
        for al in self.content:
            yield al.reverse()
        yield AlignmentPiece.Identical(self.line.asSegment())

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        for al in self:
            if al.seg_to.contains(seg):
                yield al

    def allInter(self, seg):
        for al in self:
            if al.seg_to.inter(seg):
                yield al

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.reverse()
        self.content.fireBeforeExtendRight(line, new_seq, seq)
    # alignments from new sequence to new sequence

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.reverse()
        self.content.fireBeforeCutRight(line, new_seq, pos)

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.reverse()
        self.content.fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.reverse()
        self.content.fireAfterExtendRight(line, seq)

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.reverse()
        self.content.fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line, alignments):
        # type: (Any, Correction) -> None
        self.content.fireAfterCorrect(line, alignments)
        self.reverse()
        self.content.fireAfterCorrect(line, alignments)

    def reverse(self):
        self.state = -self.state
        self.rc.state = -self.rc.state
        self.content = self.content.reverse()
        self.rc.content = self.content.rc

    def merge(self, other):
        # type: (AutoAlignmentStorage) -> AutoAlignmentStorage
        if self.state != other.state:
            self.reverse()
        res = AutoAlignmentStorage(self.line)
        res.state = self.state
        res.content.addAll(self.content.merge(other.content))
        return res

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.content.load(handler, self.line, self.line)

    def setState(self, state):
        assert state in [-1, 1]
        if self.state != state:
            self.reverse()


class RCAlignmentStorage(LineListener):

    def __init__(self, line, rc = None):
        # type: (Contig, Optional[RCAlignmentStorage]) -> None
        self.line = line
        if rc is None:
            self.content = AlignmentStorage()
            rc = RCAlignmentStorage(line.rc, self)
            rc.content = self.content.rc
        LineListener.__init__(self, rc)
        self.rc = rc # type: AutoAlignmentStorage

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        return self.content.__iter__()

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.content.getAlignmentsTo(seg)

    def allInter(self, seg):
        return self.content.allInter(seg)

    def add(self, alignment):
        self.content.add(alignment)
        self.content.add(alignment.reverse().rc)

    def addAndMergeRight(self, al):
        # type: (AlignmentPiece) -> None
        self.content.addAndMergeRight(al)
        self.content.addAndMergeLeft(al.reverse().rc)

    def addAll(self, als):
        for al in als:
            self.add(al)
        return self

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.reverse()
        self.content.fireBeforeExtendRight(line, new_seq, seq)

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.reverse()
        self.content.fireBeforeCutRight(line, new_seq, pos)
    # alignments from new sequence to new sequence

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.reverse()
        self.content.fireBeforeCorrect(alignments)

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.reverse()
        self.content.fireAfterExtendRight(line, seq)

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.reverse()
        self.content.fireAfterCutRight(line, pos)

    def fireAfterCorrect(self, line, alignments):
        # type: (Any, Correction) -> None
        self.content.fireAfterCorrect(line, alignments)
        self.reverse()
        self.content.fireAfterCorrect(line, alignments)
    # This is CRAAAZY!!! But correct.

    def reverse(self):
        self.rc.content = self.content.reverse()
        self.content = self.rc.content.rc

    def merge(self, other):
        # type: (RCAlignmentStorage) -> RCAlignmentStorage
        res = RCAlignmentStorage(self.line)
        res.content.addAll(self.content.merge(other.content))
        return res

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler):
        # type: (TokenReader) -> None
        self.content.load(handler, self.line.rc, self.line)


class TwoLineAlignmentStorage(LineListener):

    def __init__(self, line_from, line_to, rc = None, reverse = None):
        # type: (Contig, Contig, Optional[TwoLineAlignmentStorage], Optional[TwoLineAlignmentStorage]) -> None
        assert line_from.id != line_to.id and line_from.rc.id != line_to.id
        self.line_from = line_from
        self.line_to = line_to
        self.reverse = reverse
        if rc is None:
            self.content = AlignmentStorage()
            self.rc = TwoLineAlignmentStorage(line_from.rc, line_to.rc, self, None)
        else:
            self.rc = rc
            self.content = rc.content.rc # type: AlignmentStorage
        LineListener.__init__(self, self.rc)
        self.rc = self.rc # type: TwoLineAlignmentStorage
        if reverse is None and rc is None:
            reverse = TwoLineAlignmentStorage(line_to, line_from, None, self)
            self.reverse = reverse
            self.reverse.reverse = self
            self.rc.reverse = self.reverse.rc
            self.rc.reverse.reverse = self.rc

    def add(self, al):
        # type: (AlignmentPiece) -> None
        assert al.seg_from.contig == self.line_from
        assert al.seg_to.contig == self.line_to
        self.content.add(al)
        reverse = al.reverse()
        self.reverse.content.add(reverse)

    def addAll(self, als):
        for al in als:
            self.add(al)
        return self

    def __iter__(self):
        # type: () -> Generator[AlignmentPiece]
        return self.content.__iter__()

    def getAlignmentsTo(self, seg):
        # type: (Segment) -> Generator[AlignmentPiece]
        return self.content.getAlignmentsTo(seg)

    def allInter(self, seg):
        return self.content.allInter(seg)

    def normalizeReverse(self):
        self.reverse.content = self.content.reverse()
        self.reverse.rc.content = self.reverse.content.rc

    def fireBeforeExtendRight(self, line, new_seq, seq):
        # type: (Any, Contig, str) -> None
        self.content.fireBeforeExtendRight(line, new_seq, seq)
        self.normalizeReverse()

    def fireBeforeCutRight(self, line, new_seq, pos):
        # type: (Any, Contig, int) -> None
        self.content.fireBeforeCutRight(line, new_seq, pos)
        self.normalizeReverse()

    # alignments from new sequence to new sequence

    def fireBeforeCorrect(self, alignments):
        # type: (Correction) -> None
        self.content.fireBeforeCorrect(alignments)
        self.normalizeReverse()

    def fireAfterExtendRight(self, line, seq, relevant_als = None):
        # type: (Any, str, Optional[List[AlignmentPiece]]) -> None
        self.content.fireAfterExtendRight(line, seq)
        self.normalizeReverse()

    def fireAfterCutRight(self, line, pos):
        # type: (Any, int) -> None
        self.content.fireAfterCutRight(line, pos)
        self.normalizeReverse()

    def fireAfterCorrect(self, line, alignments):
        # type: (Any, Correction) -> None
        self.content.fireAfterCorrect(line, alignments)
        self.normalizeReverse()

    def addAndMergeRight(self, al):
        self.content.addAndMergeRight(al)
        self.normalizeReverse()

    def merge(self, other):
        # type: (TwoLineAlignmentStorage) -> TwoLineAlignmentStorage
        res = TwoLineAlignmentStorage(self.line_from, self.line_to)
        res.content.addAll(self.content.merge(other.content))
        res.normalizeReverse()

    def save(self, handler):
        # type: (TokenWriter) -> None
        self.content.save(handler)

    def load(self, handler, lines):
        # type: (TokenReader, Any) -> None
        self.content.load(handler, lines, lines)
        self.normalizeReverse()