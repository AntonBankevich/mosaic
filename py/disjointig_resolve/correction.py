import itertools

from typing import List, Iterable, Optional, Callable, Iterator, Generator

from common import params
from common.alignment_storage import MatchingSequence, AlignmentPiece
from common.line_align import Scorer
from common.seq_records import NamedSequence
from common.sequences import Segment, Contig


# this class stores global alignment between very long sequences.
# It only stores different parts explicitly. All the rest is expected to be the same in seq_from and seq_to
class Correction:
    def __init__(self, seq_from, seq_to, alignments):
        # type: (Contig, Contig, List[AlignmentPiece]) -> None
        self.seq_from = seq_from
        self.seq_to = seq_to
        self.alignments = alignments
        self.scorer = Scorer()

    def changeQT(self, contig_from, contig_to):
        self.alignments = [al.changeTargetContig(contig_to).changeQueryContig(contig_from) for al in self.alignments]

    def isSubstantial(self):
        # type: () -> bool
        for al in self.alignments:
            if len(list(al.splitRef())) > 1 or (len(al) > 150 and al.percentIdentity() <0.98) or (len(al.seg_from) - len(al.seg_to)) > 50:
                return True
        return False

    def mapSegmentsUp(self, segments):
        # type: (List[Segment]) -> List[Segment]
        left = self.mapPositionsUp([seg.left for seg in segments])
        right = self.mapPositionsUp([seg.right - 1 for seg in segments])
        return [Segment(self.seq_from, l, r + 1) for l, r in zip(left, right)]

    def mapSegmentsDown(self, segments):
        # type: (Iterable[Segment]) -> List[Segment]
        segments = list(segments)
        left = self.mapPositionsDown([seg.left for seg in segments])
        right = self.mapPositionsDown([seg.right - 1 for seg in segments])
        return [Segment(self.seq_to, l, r + 1) for l, r in zip(left, right)]

    def mapPositionsUp(self, positions, none_on_miss = False):
        # type: (List[int], bool) -> List[Optional[int]]
        tmp = [(pos, i) for i, pos in enumerate(positions)]
        tmp = sorted(tmp)
        res = [0] * len(positions)
        cur_pos = 0
        for al in self.alignments:
            while cur_pos < len(tmp) and tmp[cur_pos][0] <= al.seg_to.left:
                res[tmp[cur_pos][1]] = al.seg_from.left - (al.seg_to.left - tmp[cur_pos][0])
                cur_pos += 1
            for p1, p2 in al.matchingPositions(equalOnly=True):
                while cur_pos < len(positions) and tmp[cur_pos][0] <= p2:
                    if tmp[cur_pos][0] == p2 or not none_on_miss:
                        res[tmp[cur_pos][1]] = p1
                    else:
                        res[tmp[cur_pos][1]] = None
                    cur_pos += 1
        while cur_pos < len(positions):
            res[tmp[cur_pos][1]] = len(self.seq_from) - (len(self.seq_to) - tmp[cur_pos][0])
            cur_pos += 1
        return res

    def mapPositionsDown(self, positions, none_one_miss = False):
        # type: (List[int], bool) -> List[Optional[int]]
        tmp = [(pos, i) for i, pos in enumerate(positions)]
        tmp = sorted(tmp)
        res = [0] * len(positions)
        cur_pos = 0
        for al in self.alignments:
            while cur_pos < len(tmp) and tmp[cur_pos][0] <= al.seg_from.left:
                res[tmp[cur_pos][1]] = al.seg_to.left - (al.seg_from.left - tmp[cur_pos][0])
                cur_pos += 1
            for p1, p2 in al.matchingPositions(equalOnly=False):
                while cur_pos < len(positions) and tmp[cur_pos][0] <= p1:
                    if tmp[cur_pos][0] == p1 or not none_one_miss:
                        res[tmp[cur_pos][1]] = p2
                    else:
                        res[tmp[cur_pos][1]] = None
                    cur_pos += 1
        while cur_pos < len(positions):
            res[tmp[cur_pos][1]] = len(self.seq_to) - (len(self.seq_from) - tmp[cur_pos][0])
            cur_pos += 1
        return res

    def continuousMapping(self, map_function, iter):
        # type: (Callable[[List[int]], List[int]], Iterator[int]) -> Generator[int]
        chunk = []
        for item in iter:
            chunk.append(item)
            if len(chunk) > 100000:
                for res in map_function(chunk):
                    yield res
                chunk = []
        for res in map_function(chunk):
            yield res

    # This method may change the order of alignments. But they will be sorted by start.
    def composeQueryDifferences(self, als):
        # type: (List[AlignmentPiece]) -> List[AlignmentPiece]
        order = sorted(range(len(als)), key = lambda i: als[i].seg_to.left)
        # Sorting alignments into those that intersect corrections (complex) and those that do not (easy)
        easy = [] # type: List[int]
        complex = [] # type: List[int]
        cur = 0
        for al in self.alignments:
            while cur < len(als) and als[order[cur]].seg_to.left < al.seg_to.left:
                if als[order[cur]].seg_to.right >= al.seg_to.left:
                    complex.append(order[cur])
                else:
                    easy.append(order[cur])
                cur += 1
            while cur < len(als) and als[order[cur]].seg_to.left < al.seg_to.right:
                complex.append(order[cur])
                cur += 1
        while cur < len(als):
            easy.append(order[cur])
            cur += 1

        res = [None] * len(als) # type: List[AlignmentPiece]
        # Mapping alignments that do not intersect corrections
        new_easy_segs = self.mapSegmentsUp([als[i].seg_to for i in easy])
        for seg, i in zip(new_easy_segs, easy):
            res[i] = als[i].changeTargetSegment(seg)
        # Mapping alignments that intersect corrections
        func = lambda items: self.mapPositionsUp(items, True)
        matchings = [als[i].matchingSequence(True) for i in complex]
        positions = map(lambda matching: map(lambda pair: pair[1], matching), matchings)
        generator = self.continuousMapping(func, itertools.chain.from_iterable(positions))
        for i, matching in zip(complex, matchings):
            al = als[i]
            new_pairs = []
            for pos_from, pos_to in matching.matches:
                new_pos = generator.next()
                if new_pos is not None:
                    new_pairs.append((pos_from, new_pos))
            new_matching = MatchingSequence(matching.seq_from, self.seq_from.seq, new_pairs)
            corrected_matching = self.scorer.polyshMatching(new_matching, params.alignment_correction_radius)
            res[i] = corrected_matching.asAlignmentPiece(al.seg_from.contig, self.seq_from)
        return res


    @staticmethod
    def constructCorrection(alignments):
        # type: (List[AlignmentPiece]) -> Correction
        initial = alignments[0].seg_to.contig
        alignments = sorted(alignments, key = lambda al: al.seg_to.left)
        sb = []
        pos = initial.left()
        new_pos = 0
        for al in alignments:
            sb.append(initial.subSequence(pos, al.seg_to.left).seq)
            new_pos += al.seg_to.left - pos
            pos = al.seg_to.left
            sb.append(al.seg_from.Seq())
            new_pos += al.seg_from.__len__()
            pos = al.seg_to.right
        sb.append(initial.segment(alignments[-1].seg_to.right,initial.right()).Seq())
        new_pos += initial.right() - alignments[-1].seg_to.right
        new_seq = Contig("".join(sb), "TMP1_" + initial.id)
        new_als = []
        pos = initial.left()
        new_pos = 0
        for al in alignments:
            new_pos += al.seg_to.left - pos
            new_seg_from = Segment(new_seq, new_pos, new_pos + al.seg_from.__len__())
            new_als.append(al.changeQuerySegment(new_seg_from))
            pos = al.seg_to.right
            new_pos += al.seg_from.__len__()
        return Correction(new_seq, initial, new_als)