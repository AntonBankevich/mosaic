import sys

from typing import Dict, List, Iterator, Optional, Iterable, BinaryIO, Tuple, Generator

from alignment.align_tools import Aligner
from common import SeqIO, params, basic
from common.alignment_storage import AlignmentPiece, ReadCollection
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from common.sequences import ContigStorage, UniqueList, Contig, ContigCollection, Segment
from disjointig_resolve.accurate_line import NewLine, ExtensionHandler, Knot
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.line_alignments import TwoLineAlignmentStorage
from disjointig_resolve.smart_storage import AlignmentStorage


class NewLineStorage(ContigStorage):
    def __init__(self, disjointigs, aligner):
        # type: (DisjointigCollection, Aligner) -> None
        ContigStorage.__init__(self, [], True)
        self.disjointigs = disjointigs
        self.aligner = aligner
        self.items = dict() # type: Dict[str, NewLine]
        self.cnt = 1
        self.listeners = [] # type: List[LineStorageListener]
        self.name_printer = None

    def __iter__(self):
        # type: () -> Iterator[NewLine]
        return self.items.values().__iter__()

    def __getitem__(self, item):
        # type: (str) -> NewLine
        return self.items[item]

    def addListener(self, listener):
        self.listeners.append(listener)

    def removeListener(self, listener):
        self.listeners.remove(listener)

    def notifyMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        for listener in self.listeners:
            listener.FireMergedLines(al1, al2)

    def notifySplitLine(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        for listener in self.listeners:
            listener.FireSplitLine(al1, al2)

    def addNew(self, seq, name = None):
        # type: (str, Optional[str]) -> NewLine
        if name is None:
            name = "L" + str(self.cnt)
            self.cnt += 1
        else:
            if not basic.isCanonocal(name):
                name = basic.Reverse(basic.Reverse(name))
        new_line = NewLine(seq, str(name), ExtensionHandler(self.disjointigs, self.aligner))
        self.add(new_line)
        new_line.name_printer = self.name_printer
        new_line.rc.name_printer = self.name_printer
        return new_line

    def fillFromContigs(self, contigs):
        # type: (Iterable[Contig]) -> None
        for contig in UniqueList(contigs):
            line = self.addNew(contig.seq, "L" + contig.id)
            line.initial.add(AlignmentPiece.Identical(contig.asSegment(), line.asSegment()))

    def splitFromContigs(self, contigs, max_contig = 50000, cut_size = 20000):
        # type: (ContigStorage, int, int) -> None
        for contig in contigs.unique():
            if not basic.isCanonocal(contig.id):
                contig = contig.rc
            if len(contig) > max_contig:
                line1 = self.addNew(contig.seq[:cut_size], "L" + contig.id + "l")
                line2 = self.addNew(contig.seq[-cut_size:], "L" + contig.id + "r")
                line1.initial.add(AlignmentPiece.Identical(contig.asSegment().prefix(length=cut_size), line1.asSegment()))
                line2.initial.add(AlignmentPiece.Identical(contig.asSegment().suffix(length=cut_size), line2.asSegment()))
                line1.tie(line2, len(contig) - 2 * cut_size, contig.seq[cut_size:-cut_size])
            else:
                line = self.addNew(contig.seq, "L" + contig.id)
                line.initial.add(AlignmentPiece.Identical(contig.asSegment(), line.asSegment()))

    def alignDisjointigs(self):
        for line in self:
            line.disjointig_alignments.clean()
        for al in self.aligner.dotplotAlign(self.disjointigs.unique(), self):
            if len(al) > params.k - 500:
                line = al.seg_to.contig # type: NewLine
                line.disjointig_alignments.add(al)

    def mergeLines(self, alignment, k):
        # type: (AlignmentPiece, int) -> NewLine
        sys.stdout.trace("Line operation Merge", alignment.seg_from.contig, alignment.seg_to.contig, alignment)
        line1 = alignment.seg_from.contig #type: NewLine
        line2 = alignment.seg_to.contig #type: NewLine
        assert line1 != line2
        if len(alignment) < k + 100:
            sys.stdout.trace("Prolonging line to ensure alignment of at least k")
            seg = line2.segment(alignment.seg_to.right, alignment.seg_to.right + k + 100 - len(alignment))
            line1.extendRight(seg.Seq())
            alignment = alignment.mergeDistant(AlignmentPiece.Identical(line1.asSegment().suffix(length=len(seg)), seg))
        # Cutting hanging tips of both lines
        al_storage = AlignmentStorage()
        al_storage.add(alignment)
        storage = TwoLineAlignmentStorage(line1, line2)
        line2.addListener(storage)
        line1.addListener(storage.reverse)
        storage.add(alignment)
        if alignment.seg_from.right < len(line1):
            line1.cutRight(alignment.seg_from.right)
            sys.stdout.trace( "Cut right")
            sys.stdout.trace( list(storage.content)[0])
            sys.stdout.trace( "\n".join(list(storage.content)[0].asMatchingStrings()))
            sys.stdout.trace( list(storage.content)[0].cigar)
        if alignment.seg_to.left > 0:
            line2.rc.cutRight(len(line2) - alignment.seg_to.left)
            sys.stdout.trace( "Cut left")
            sys.stdout.trace( list(storage.content)[0])
            sys.stdout.trace( "\n".join(list(storage.content)[0].asMatchingStrings()))
            sys.stdout.trace( list(storage.content)[0].cigar)
        alignment = list(storage.content)[0] # type: AlignmentPiece
        line2.removeListener(storage)
        line1.removeListener(storage.reverse)

        # Making sure line sequences match on the overlap
        if alignment.seg_from.left > 0:
            new_seq = Contig(line1.asSegment().prefix(pos=alignment.seg_from.left).Seq() + line2.seq, "new_seq")
        else:
            new_seq = Contig(line2.seq, "new_seq")
        al2 = AlignmentPiece.Identical(line2.asSegment(), new_seq.asSegment().suffix(length=len(line2)))
        sys.stdout.trace( "Al2:", al2)
        alignment = alignment.compose(al2).reverse()
        sys.stdout.trace( "Composed alignment", alignment)
        sys.stdout.trace("\n".join(alignment.asMatchingStrings()))
        sys.stdout.trace( alignment.cigar)
        assert alignment.seg_to.right == len(line1)
        assert alignment.seg_from.left == al2.seg_to.left
        line1.correctSequence([alignment])

        # Now lines have exact match
        name = "(" + ",".join(basic.parseLineName(line1.id) + basic.parseLineName(line2.id))  + ")"
        line = self.addNew(new_seq.seq, name)
        assert line.seq.startswith(line1.seq)
        assert line.seq.endswith(line2.seq)
        al1 = AlignmentPiece.Identical(line1.asSegment(), line.asSegment().prefix(length=len(line1)))
        al2 = AlignmentPiece.Identical(line2.asSegment(), line.asSegment().suffix(length=len(line2)))

        line.initial.addAll(line1.initial.targetAsSegment(al1.seg_to).merge(line2.initial.targetAsSegment(al2.seg_to)))
        line.correct_segments.addAll(line1.correct_segments.contigAsSegment(al1.seg_to).
                                     merge(line2.correct_segments.contigAsSegment(al2.seg_to)))
        line.completely_resolved.addAll(line1.completely_resolved.contigAsSegment(al1.seg_to).
                                     merge(line2.completely_resolved.contigAsSegment(al2.seg_to), k))
        line.disjointig_alignments.addAll(line1.disjointig_alignments.targetAsSegment(al1.seg_to).
                                          merge(line2.disjointig_alignments.targetAsSegment(al2.seg_to)))
        for al in line1.read_alignments.targetAsSegment(al1.seg_to).merge(line2.read_alignments.targetAsSegment(al2.seg_to)):
            line.addReadAlignment(al)
        line1.cleanReadAlignments()
        line2.cleanReadAlignments()


        self.notifyMergedLines(al1, al2)
        knot_right = line2.knot
        knot_left = line1.rc.knot
        self.remove(line1)
        self.remove(line2)
        if knot_right is not None:
            if knot_right.line_right == line1:
                line.tie(line, knot_right.gap, knot_right.gap_seq)
            else:
                line.tie(knot_right.line_right, knot_right.gap, knot_right.gap_seq)
        if knot_left is not None and knot_left.line_right != line2.rc:
            line.rc.tie(knot_left.line_right, knot_left.gap, knot_left.gap_seq)
        return line

    def splitLine(self, seg):
        # type: (Segment) -> Tuple[NewLine, NewLine]
        sys.stdout.trace("Line operation Split", seg)
        line = seg.contig # type: NewLine
        seg1 = line.asSegment().prefix(pos=seg.right)
        line1 = self.addNew(seg1.Seq(), line.id + "l")
        seg2 = line.asSegment().suffix(pos=seg.left)
        line2 = self.addNew(seg2.Seq(), line.id + "r")
        al1 = AlignmentPiece.Identical(seg1, line1.asSegment())
        al2 = AlignmentPiece.Identical(seg2, line2.asSegment())
        line1.initial.addAll([al.embed(al1) for al in line.initial.allInter(seg1, params.min_alignment_size)])
        line2.initial.addAll([al.embed(al2) for al in line.initial.allInter(seg2, params.min_alignment_size)])
        line1.correct_segments.addAll(line.correct_segments.cap(seg=seg1, min_inter=params.k).map(al1))
        line2.correct_segments.addAll(line.correct_segments.cap(seg=seg2, min_inter=params.k).map(al2))
        line1.completely_resolved.addAll(line.completely_resolved.cap(seg=seg1, min_inter=params.k).map(al1).filterBySize(min=params.k))
        line2.completely_resolved.addAll(line.completely_resolved.cap(seg=seg2, min_inter=params.k).map(al2).filterBySize(min=params.k))

        line1.disjointig_alignments.addAll([al.embed(al1) for al in line.disjointig_alignments.allInter(seg1, params.k)])
        line2.disjointig_alignments.addAll([al.embed(al2) for al in line.disjointig_alignments.allInter(seg2, params.k)])
        for al in line.read_alignments:
            if al.seg_to.interSize(seg1) > params.k:
                line1.addReadAlignment(al.embed(al1))
        for al in line.read_alignments:
            if al.seg_to.interSize(seg2) > params.k:
                line2.addReadAlignment(al.embed(al2))
        line.cleanReadAlignments()
        self.notifySplitLine(al1, al2)
        self.remove(line)
        if line.knot is not None:
            line2.tie(line.knot.line_right, line.knot.gap, line.knot.gap_seq)
        if line.rc.knot is not None:
            line1.rc.tie(line.rc.knot.line_right, line.rc.knot.gap, line.rc.knot.gap_seq)
        return line1, line2


    def printToFile(self, handler):
        # type: (BinaryIO) -> None
        for line in self.items:# type: NewLine
            handler.write(line.__str__() + "\n")

    def printToFasta(self, handler):
        # type: (BinaryIO) -> None
        for line in UniqueList(self.items.values()):
            SeqIO.write(line, handler, "fasta")

    def printKnottedToFasta(self, handler):
        # type: (BinaryIO) -> None
        printed = set()
        cnt = 1
        for chain in self.chains():
            if chain[0].rc.id in printed:
                continue
            for line in chain:
                printed.add(line.id)
            seq = []
            id = []
            if chain[-1].knot is not None:
                id.append("Circular")
            for line in chain:
                id.append(line.id)
                if line.knot is not None:
                    id.append(str(line.knot.gap))
                    if line.knot.gap < 0:
                        seq.append(line.seq[:line.knot.gap])
                    else:
                        seq.append(line.seq)
                        seq.append(line.knot.gap_seq)
                else:
                    seq.append(line.seq)
            sys.stdout.trace( cnt, ":", ";".join(id))
            SeqIO.write(NamedSequence("".join(seq), "contig_" + str(cnt)), handler, "fasta")
            cnt += 1

    def chains(self):
        # type: () -> Generator[List[NewLine]]
        visited = set()
        for line in self.items.values():
            if line.id in visited or line.rc.knot is not None:
                continue
            res = [line]
            visited.add(line.id)
            while line.knot is not None:
                line = line.knot.line_right
                res.append(line)
                visited.add(line.id)
            yield res
        for line in self.items.values():
            if line.id in visited:
                continue
            res = [line]
            visited.add(line.id)
            line = line.knot.line_right
            while line != res[0]:
                res.append(line)
                visited.add(line.id)
                line = line.knot.line_right
            yield res

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(str(self.cnt))
        line_ids = map(lambda line: line.id, UniqueList(self.items.values()))
        handler.writeTokens(line_ids)
        for line_id in line_ids:
            line = self.items[line_id]
            handler.writeTokenLine(line.seq)
        for line_id in line_ids:
            line = self.items[line_id]
            line.save(handler)
        for line_id in line_ids:
            line = self.items[line_id]
            if line.knot is not None:
                handler.writeToken("Knot")
                line.knot.save(handler)
            else:
                handler.writeToken("None")
            if line.rc.knot is not None:
                handler.writeToken("Knot")
                line.rc.knot.save(handler)
            else:
                handler.writeToken("None")

    def load(self, handler, reads, contigs):
        # type: (TokenReader, ReadCollection, ContigCollection) -> None
        self.cnt = int(handler.readToken())
        keys = list(handler.readTokens())
        for key in keys:
            self.addNew(handler.readToken(), key)
        for key in keys:
            line = self.items[key]
            line.loadLine(handler, self.disjointigs, reads, contigs)
        for key in keys:
            line = self.items[key]
            if handler.readToken() is not None:
                knot = Knot.load(handler, self)
                line.knot = knot
                line.knot.line_right.rc.knot = knot.rc
            line = line.rc
            if handler.readToken() is not None:
                knot = Knot.load(handler, self)
                line.knot = knot
                line.knot.line_right.rc.knot = knot.rc

    def remove(self, line):
        if line.knot is not None:
            line.unTie()
        if line.rc.knot is not None:
            line.rc.unTie()
        line.cleanReadAlignments()
        del self.items[line.id]
        del self.items[line.rc.id]

    def removeLine(self, line):
        self.remove(line)
        for listener in self.listeners:
            listener.FireRemoveLine(line)

class LineStorageListener:
    def __init__(self):
        pass

    def FireMergedLines(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        pass

    def FireSplitLine(self, al1, al2):
        # type: (AlignmentPiece, AlignmentPiece) -> None
        pass

    def FireRemoveLine(self, line):
        # type: (NewLine) -> None
        pass
