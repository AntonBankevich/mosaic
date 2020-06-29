import itertools
import sys

from typing import Iterator, Tuple, Optional, List

from alignment.polishing import Polisher
from common import params
from common.alignment_storage import AlignedRead, AlignmentPiece
from common.line_align import Scorer
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.line_storage import NewLineStorage


class LineMerger:
    def __init__(self, storage, polisher, dot_plot):
        # type: (NewLineStorage, Polisher, LineDotPlot) -> None
        self.storage = storage
        self.polisher = polisher
        self.dot_plot = dot_plot

    class Record:
        def __init__(self, al1, al2):
            # type: (AlignmentPiece, AlignmentPiece) -> None
            self.read = al1.seg_from.contig # type: AlignedRead
            line = al1.seg_to.contig # type: NewLine
            if al1.seg_from.inter(al2.seg_from):
                if al1.seg_from.interSize(al2.seg_from) > 100:
                    tmp = al1.composeTargetDifference(al2)
                    self.gap = -tmp.seg_to.left - tmp.rc.seg_from.right
                else:
                    self.gap = al1.matchingSequence().mapDown(al2.seg_from.left, roundDown=False) - len(al1.seg_to.contig) - al2.seg_to.left
            else:
                self.gap = len(self.read) - al1.seg_from.right - (len(line) - al1.seg_to.right) - (len(self.read) - al2.seg_from.left) - al2.seg_to.left
            self.al1 = al1
            self.al2 = al2
            self.other = al2.seg_to.contig # type: NewLine
            self.initial_gap = line.rc.initial[0].seg_to.left + self.gap + self.other.initial[0].seg_to.left

        def __repr__(self):
            return str([self.other, self.gap, self.al1, self.al2])

        def __str__(self):
            return str([self.other, self.gap, self.al1, self.al2])

    # Find connection of line to any other line using reads. Line is supposed to contain or precede the other line.
    def tryMergeRight(self, line):
        # type: (NewLine) -> Optional[NewLine]
        assert line.read_alignments.checkLine(line), str(line.read_alignments)
        if line.circular or line.knot is not None:
            return None
        read_alignments = line.read_alignments.allInter(line.asSegment().suffix(length=1000))
        candidates = [] # type: List[LineMerger.Record]
        for al1 in read_alignments:
            read = al1.seg_from.contig # type: AlignedRead
            if al1.contradictingRTC():
                continue
            for al2 in read.alignments:
                if al2.contradictingRTC() or (al1.canMergeTo(al2) and al1.deepInter(al2)):
                    continue
                new_rec = self.Record(al1, al2)
                if len(line) + new_rec.other.initial[0].seg_to.left + new_rec.gap > line.initial[0].seg_to.left:
                    candidates.append(new_rec)
        #             (al2.seg_from.contig, gap, read, al1, al2)
        candidates = sorted(candidates, key = lambda rec: rec.other.id)
        final_candidates = []
        for other_line, iter in itertools.groupby(candidates, lambda rec: rec.other): # type: NewLine, Iterator[LineMerger.Record]
            recs = list(iter)
            if recs[-1].gap - recs[0].gap > min(100, abs(recs[-1].gap) / 8):
                sys.stdout.trace("\n".join(map(str, candidates)))
                sys.stdout.warn("WARNING: Ambiguous knotting to the same line")
                if len(recs) >= 5 and recs[-2].gap - recs[1].gap < min(10, abs(recs[-2].gap) / 10):
                    recs = recs[1:-1]
                else:
                    return None
                # assert False, "Ambiguous knotting to the same line" + str(recs[0])
            avg = sum([rec.gap for rec in recs]) / len(recs)
            avg_initial = sum([rec.initial_gap for rec in recs]) / len(recs)
            if recs[0].other == line.rc:
                sys.stdout.warn( "WARNING: Ignoring connection to reverse-compliment line")
                continue
            final_candidates.append((avg, recs[0].other, recs, avg_initial))
        if len(final_candidates) == 0:
            return None
        final_candidates = sorted(final_candidates, key = lambda candidate: candidate[-1])
        final = final_candidates[0]
        if len(final_candidates) > 1:
            sys.stdout.warn("Extra candidates")
            sys.stdout.trace("\n".join(map(str, candidates)))
        # for candidate in final_candidates[1:]:
        #     if final[0] + len(final[1]) > candidate[0]:
        #         print "\n".join(map(str, candidates))
        #         assert False, "Contradicting candidates" + str(final[0]) + " " + str(final[1]) + " " + str(candidate[0]) + " " + str(candidate[1])
        if final[0] > 0:
            sys.stdout.trace("Positive gap. Can not merge line", line)
            sys.stdout.trace(final)
            return None
        elif len(final[2]) <= 1:
            sys.stdout.trace("Insufficient support to merge line", line)
            sys.stdout.trace(final)
            return None
        else:
            sys.stdout.info("Merging", line, "with", final[1], "with gap", final[0])
            sys.stdout.trace("Alignments:")
            sys.stdout.trace("\n".join(map(str, final[2])))
            other = final[1]
            assert line != other.rc
            assert other.rc.knot is None
            line_alignment = final[2][0].al1.composeTargetDifference(final[2][0].al2)
            sys.stdout.trace( "Alignment:", line_alignment)
            sys.stdout.trace( "\n".join(line_alignment.asMatchingStrings()))
            sys.stdout.trace( line_alignment.cigar)
            sys.stdout.trace( list(self.dot_plot.getAlignmentsToFrom(other, line)))
            tmp = None
            if final[0] < -params.k - 100:
                for al in self.dot_plot.getAlignmentsToFrom(other, line):
                    if len(list(al.matchingSequence().common(line_alignment.matchingSequence()))) > 0:
                        if tmp is None or len(tmp) < len(al):
                            tmp = al
                if tmp is None:
                    sys.stdout.warn("No good line alignment found. Alignment based on reads will be used.")
                else:
                    sys.stdout.trace("Switched to line alignment:", tmp)
                if (tmp.seg_to.left < 20 and tmp.rc.seg_to.left < 20) or \
                        (tmp.seg_from.left < 20 and tmp.rc.seg_from.left < 20) or \
                        (tmp.seg_from.left < 20 and tmp.rc.seg_to.left < 20):
                    sys.stdout.warn( "One line is substring of another.", str(line_alignment) + " " + str(tmp))
                elif tmp.seg_to.left > 30 or tmp.rc.seg_from.left > 30:
                    sys.stdout.warn("Line alignment is not overlap!", tmp)
                    if params.strict_merging_alignment:
                        assert tmp.seg_to.left < 30 and tmp.rc.seg_from.left < 30, str(line_alignment) + " " + str(tmp)
                line_alignment = tmp
            pref = line_alignment.seg_from.left
            suff = len(line_alignment.seg_to.contig) - line_alignment.seg_to.right
            line_alignment = Scorer().polyshAlignment(line_alignment, params.alignment_correction_radius)
            sys.stdout.trace("Polished alignment:", line_alignment)
            sys.stdout.trace("\n".join(line_alignment.asMatchingStrings()))
            sys.stdout.trace(line_alignment.cigar)
            if line == other:
                gap = -line_alignment.rc.seg_from.right - line_alignment.seg_to.left + line.correct_segments[0].left + line.rc.correct_segments[0].left
                if gap > 0:
                    sys.stdout.trace("Line is circular but not ready for completion. Skipping.")
                    return None
                line.cutRight(line.correct_segments[-1].right)
                line.rc.cutRight(line.rc.correct_segments[-1].right)
                line.tie(line, gap, "")
                sys.stdout.info(line, "is circular")
                return line
            new_line = self.storage.mergeLines(line_alignment, params.k)
            seg = new_line.segment(pref, len(new_line) - suff)
            correction = self.polisher.polishSegment(seg, list(new_line.read_alignments.allInter(seg)))
            new_line.correctSequence([correction])
            new_line.updateCorrectSegments(new_line.segment(pref, len(new_line) - suff).expand(100))
            return new_line


