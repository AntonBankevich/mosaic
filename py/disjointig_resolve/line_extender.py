import itertools
import sys

from typing import Optional, Tuple, List, Dict, Iterable, Set, Generator, Union, Iterator

from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic, params
from common.line_align import Scorer
from common.sequences import Segment, ContigStorage
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from disjointig_resolve.accurate_line import NewLine, LinePosition
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage



# LineExtender contains algorithms for resolving repeats by extending lines.
# The results of its work is reflected in lines themselves as they grow and knot to each other
class LineExtender:
    def __init__(self, aligner, knotter, disjointigs, dot_plot):
        # type: (Aligner, LineMerger, DisjointigCollection, LineDotPlot) -> None
        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        self.knotter = knotter
        self.disjointigs = disjointigs
        self.dot_plot = dot_plot
        self.scorer = Scorer()

    def processLine(self, line):
        # type: (NewLine) -> int
        line.completely_resolved.mergeSegments(params.k)
        bound = LinePosition(line, line.left())
        new_recruits = 0
        new_line = self.knotter.tryMergeRight(line)
        if new_line is not None:
            self.updateAllStructures(list(new_line.completely_resolved))
            return 1
        self.updateAllStructures(line.completely_resolved)
        while True:
            seg_to_resolve = line.completely_resolved.find(bound.suffix(), params.k)
            if seg_to_resolve is None:
                break
            if line.knot is not None and seg_to_resolve.right == len(line):
                break
            if seg_to_resolve.right <= line.initial[0].seg_to.left + params.k:
                bound = LinePosition(line, seg_to_resolve.right - params.k + 1)
                continue
            result = self.attemptCleanResolution(seg_to_resolve)
            total = sum([len(arr) for seg, arr in result])
            new_recruits += total
            if total == 0:
                bound = LinePosition(line, seg_to_resolve.right - params.k + 1)
                continue
            self.updateAllStructures([seg for seg, arr in result])
            new_line = self.knotter.tryMergeRight(line)
            if new_line is not None:
                self.updateAllStructures(list(new_line.completely_resolved))
                return new_recruits + 1
        return new_recruits

    # input: a collection of segments that had reads recruited to.
    def updateAllStructures(self, interesting_segments):
        # type: (Iterable[Segment]) -> None
        sys.stdout.trace("Updating structures:", interesting_segments)
        # Correct contig sequences, update correct segment storages. Return segments that were corrected.
        corrected = self.correctSequences(interesting_segments)
        # Collect all relevant contig segments, collect all reads that align to relevant segments.
        # Mark resolved bound for each read.
        sys.stdout.trace("Expanding resolved segments:")
        records = self.collectRecords(corrected)  # type: List[LineExtender.Record]
        for rec in records:
            sys.stdout.trace("Record:", rec.line, rec.correct, rec.resolved)
            sys.stdout.trace("Reads from record:")
            for al in rec:
                sys.stdout.trace(al, al.seg_from.contig.alignments)
            sys.stdout.trace(rec.reads)
            sys.stdout.trace(rec.potential_good)
        # Update resolved segments on all relevant contig positions
        self.updateResolved(records)

    def updateResolved(self, records):
        # type: (List[LineExtender.Record]) -> None
        ok = True
        while ok:
            sys.stdout.trace( "Good reads:")
            rec = records[0]  # type: LineExtender.Record
            for read_name in rec.good_reads:
                sys.stdout.trace( read_name, rec.read_bounds[read_name])
            ok = False
            for rec in records:
                if self.attemptProlongResolved(rec):
                    sys.stdout.trace( "Successfully prolonged resolved:", rec.line, rec.line.initial, rec.resolved, rec.line.completely_resolved)
                    ok = True
        for rec in records:
            line = rec.resolved.contig  # type: NewLine
            line.completely_resolved.add(rec.resolved)
            for seg in rec.old_resolved:
                line.completely_resolved.add(seg)
            line.completely_resolved.mergeSegments(params.k - 1)

    def collectRecords(self, corrected):
        # type: (List[Segment]) -> List[LineExtender.Record]
        sys.stdout.trace("Collecting records", corrected)
        read_bounds = dict()
        records = dict()  # type: Dict[Segment, LineExtender.Record]
        good_reads = set()
        for seg in corrected:
            sys.stdout.trace("Oppa initial:", seg)
            seg = seg.expandLeft(params.k)
            sys.stdout.trace("Alignments relevant for", seg, list(self.dot_plot.allInter(seg)))
            for al in self.dot_plot.allInter(seg):
                seg1 = al.matchingSequence().mapSegUp(al.seg_from.contig, seg)
                line = al.seg_from.contig  # type:NewLine
                for seg_correct in line.correct_segments.allInter(al.seg_from):
                    for seg_resolved in line.completely_resolved.allInter(seg_correct):
                        if seg_resolved in records:
                            continue
                        if seg_resolved.right == len(line):
                            next_start = len(line)
                        else:
                            next = line.completely_resolved.find(line.asSegment().suffix(pos=seg_resolved.right), 1)
                            if next is None:
                                next_start = len(line)
                            else:
                                next_start = next.left
                        next_start = min(next_start, len(line) - 200)
                        focus = line.segment(max(seg_resolved.left, min(seg_resolved.right - params.k, seg1.left)),
                                             min(seg_correct.right, next_start + params.k))
                        als = list(line.getRelevantAlignmentsFor(focus))
                        reads  = ContigStorage()
                        for al in als:
                            reads.add(al.seg_from.contig)
                        als = list(self.aligner.localAlign(reads.unique(), ContigStorage([line])))
                        final_als = []
                        sys.stdout.trace( "Focus:", focus, seg_resolved)
                        sys.stdout.trace( als)
                        for al in als:
                            if al.seg_to.contig == line.rc:
                                al = al.rc
                            if al.seg_to.interSize(focus) >= params.k - 100:
                                final_als.append(al)
                        sys.stdout.trace( final_als)
                        sys.stdout.trace("Finished realignment of reads")
                        records[seg_resolved] = self.createRecord(seg_resolved, next_start, seg_correct, final_als, good_reads, read_bounds)
        records = list(records.values())  # type: List[LineExtender.Record]
        return records

    def correctSequences(self, interesting_segments):
        # type: (Iterable[Segment]) -> List[Segment]
        interesting_segments = list(interesting_segments)
        to_correct = []
        for seg in interesting_segments:
            line = seg.contig # type: NewLine
            correct = line.correct_segments.find(seg)
            next = line.correct_segments.find(line.suffix(correct.right), 1)
            if next is None:
                right = len(line)
            else:
                right = min(len(line), next.left + params.k / 2)
            to_correct.append(line.segment(correct.right - params.k / 2, right))
        to_correct = sorted(to_correct, key=lambda seg: (basic.Normalize(seg.contig.id), seg.left))
        corrected = []
        for line_id, it in itertools.groupby(to_correct,
                                          key=lambda seg: basic.Normalize(seg.contig.id)):  # type: NewLine, Iterable[Segment]
            it = list(it)
            line = None # type: NewLine
            forward = SegmentStorage()
            backward = SegmentStorage()
            for seg in it:
                if seg.contig.id != line_id:
                    backward.add(seg)
                    line = seg.contig.rc
                else:
                    forward.add(seg)
                    line = seg.contig
            to_polysh = SegmentStorage()
            to_polysh.addAll(forward).addAll(backward.rc)
            to_polysh.mergeSegments()
            line.addListener(to_polysh)
            line.addListener(forward)
            line.rc.addListener(backward)
            sys.stdout.trace( "Polishing:", to_polysh)
            if to_polysh[-1].RC().left < 200:
                l = to_polysh[-1].right
                if self.attemptExtend(line):
                    to_polysh.add(line.asSegment().suffix(pos=l))
                    forward.add(line.asSegment().suffix(pos=l))
            if to_polysh[0].left < 200:
                l = to_polysh[0].RC().right
                if self.attemptExtend(line.rc):
                    to_polysh.rc.add(line.rc.asSegment().suffix(pos=l))
                    backward.add(line.rc.asSegment().suffix(pos=l))
            to_polysh.mergeSegments()
            forward.mergeSegments()
            backward.mergeSegments()
            line.removeListener(to_polysh)
            new_segments = self.polyshSegments(line, to_polysh)
            line.removeListener(forward)
            line.rc.removeListener(backward)
            corrected.extend(forward)
            corrected.extend(backward)
            self.updateCorrectSegments(line)
        return corrected

    def attemptCleanResolution(self, resolved):
        # type: (Segment) -> List[Tuple[Segment, List[AlignmentPiece]]]
        # Find all lines that align to at least k nucls of resolved segment. Since this segment is resolve we get all
        sys.stdout.trace("Attempting recruitment:", resolved, resolved.contig, resolved.contig.correct_segments)
        resolved = resolved.suffix(length = min(len(resolved), params.k * 2))
        sys.stdout.trace( "Considering resolved subsegment:", resolved)
        line_alignments = filter(lambda al: len(al.seg_to) >= params.k and resolved.interSize(al.seg_to) > params.k - 30,
                                 self.dot_plot.allInter(resolved)) # type: List[AlignmentPiece]
        line_alignments = [al for al in line_alignments if
                           (al.seg_from.right >= al.seg_from.contig.initial[0].seg_to.left + params.k + 20
                           and al.seg_to.right >= al.seg_to.contig.initial[0].seg_to.left + params.k + 20)
                           or al.isIdentical()]
        sys.stdout.trace( "Alternative lines:", map(str, line_alignments))
        for al in line_alignments:
            if not al.isIdentical():
                sys.stdout.trace( al)
                sys.stdout.trace( "\n".join(al.asMatchingStrings()))
        line_alignments = [al.reduce(target=resolved) for al in line_alignments]
        read_alignments = [] # type: List[Tuple[AlignmentPiece, Segment]]
        correct_segments = []
        active_segments = set()
        for ltl in line_alignments:
            line = ltl.seg_from.contig # type: NewLine
            new_copy = line.correct_segments.find(ltl.seg_from)
            # assert new_copy is not None and new_copy.interSize(ltl.seg_from) >= max(len(ltl.seg_from) - 20, params.k), str([ltl, new_copy, str(line.correct_segments)])
            # assert new_copy is not None, str([ltl, line.correct_segments])
            if new_copy is None:
                return []
            if not new_copy.contains(ltl.seg_from):
                sys.stdout.trace( "Warning: alignment of resolved segment to uncorrected segment")
                sys.stdout.trace( ltl, new_copy, line.correct_segments)
            correct_segments.append(new_copy)
            if ltl.percentIdentity() > 0.95:
                active_segments.add(new_copy)
            read_alignments.extend(zip(line.getRelevantAlignmentsFor(ltl.seg_from), itertools.cycle([correct_segments[-1]])))
        read_alignments = sorted(read_alignments, key=lambda al: al[0].seg_from.contig.id)
        alignments_by_read = itertools.groupby(read_alignments, lambda al: al[0].seg_from.contig.id)
        new_recruits = []
        sys.stdout.trace("Starting read recruitment to", map(str, line_alignments))
        for name, it in alignments_by_read:
            als = list(it) # type: List[Tuple[AlignmentPiece, Segment]]
            read = als[0][0].seg_from.contig # type: AlignedRead
            sys.stdout.trace( "Recruiting read:", read, als)
            ok = False
            for al in als:
                if al[0].seg_to.interSize(resolved) >= params.k:
                    ok = True
                    break
            if not ok:
                sys.stdout.trace( "Read does not overlap with resolved", resolved)
                continue

            skip = False
            for al1 in als:
                for al2 in read.alignments:
                    if al1[0].seg_to.inter(al2.seg_to):
                        sys.stdout.trace( "Read already recruited", al1, al2)
                        skip = True
                        break
                if skip:
                    break
            if skip:
                continue
            new_als = []
            for al in als:
                if not al[0].contradictingRTC(tail_size=params.bad_end_length):
                    new_als.append((self.scorer.polyshAlignment(al[0], params.alignment_correction_radius), al[1]))
            if len(new_als) == 0:
                sys.stdout.warn("No noncontradicting alignments of a read")
                winner = None
                seg = None
            else:
                winner, seg = self.tournament(new_als) #type: AlignmentPiece, Segment
            if winner is None:
                sys.stdout.trace( "No winner")
            else:
                sys.stdout.trace("Winner for", winner.seg_from.contig.id, ":", winner, seg)
            if winner is not None:
                if seg not in active_segments:
                    sys.stdout.trace("Winner ignored since winning segment is too different from investigated segment")
                elif winner.percentIdentity() < 0.85:
                    sys.stdout.trace( "Winner ignored since it is too different from winning line")
                else:
                    line = winner.seg_to.contig # type: NewLine
                    line.addReadAlignment(winner)
                    new_recruits.append((seg, winner))
        new_recruits = sorted(new_recruits, key = lambda rec: (rec[0].contig.id, rec[0].left, rec[0].right))
        sys.stdout.info("Recruited " + str(len(new_recruits)) + " new reads")
        return [(seg, [al for seg, al in it]) for seg, it in itertools.groupby(new_recruits, key = lambda rec: rec[0])]

    def fight(self, c1, c2):
        # type: (Tuple[AlignmentPiece, Segment], Tuple[AlignmentPiece, Segment]) -> Optional[Tuple[AlignmentPiece, Segment]]
        assert c1[0].seg_from.contig == c2[0].seg_from.contig
        s1, s2, s12 = self.scorer.scoreInCorrectSegments(c1[0], c1[1], c2[0], c2[1])
        if s1 is not None and s2 is not None:
            diff = abs(s1 - s2)
        else:
            diff = None
        if s12 is None:
            if s1 is None:
                winner = c2
            else:
                winner = c1
        else:
            if s12 < 25 or (s12 < 40 and abs(s1 - s2) < s12 * 0.8) or (s12 < 100 and abs(s1 - s2) < s12 * 0.5) or abs(s1 - s2) < s12 * 0.3:
                winner = None
            elif s1 > s2:
                winner = c2
            else:
                winner = c1
        if winner is None:
            sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", diff, s12, s1, s2, "No winner")
        else:
            sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", diff, s12, s1, s2, "Winner:", winner)
        return winner

    def tournament(self, candidates):
        # type: (List[Tuple[AlignmentPiece, Segment]]) -> Tuple[Optional[AlignmentPiece], Optional[Segment]]
        best = None
        best_id = None
        wins = []
        for i, candidate in enumerate(candidates):
            if best is None:
                best = candidate
                best_id = i
            else:
                best = self.fight(candidate, best)
                if best is None:
                    best_id = None
                    wins = []
                elif best == candidates[best_id]:
                    wins.append(i)
                else:
                    best_id = i
                    wins = []
        if best is None:
            return None, None
        if len(candidates) > 2:
            for i, candidate in enumerate(candidates):
                if i == best_id or i in wins:
                    continue
                fight_results = self.fight(candidate, best)
                if fight_results is None or fight_results != best:
                    return None, None
        return best

    def attemptExtend(self, line):
        # type: (NewLine) -> bool
        sys.stdout.trace( "Attempting to extend:", line)
        if line.knot is not None:
            sys.stdout.trace( "Blocked by knot")
            return False
        relevant_reads = list(line.read_alignments.allInter(line.asSegment().suffix(length=min(params.k, len(line) - 20))))
        sys.stdout.trace( "Relevent reads for extending", relevant_reads)
        if len(relevant_reads) == 0:
            return False
        new_contig, relevant_als = self.polisher.polishEnd(relevant_reads)
        if len(new_contig) == len(line):
            return False
        assert line.seq == new_contig.prefix(len=len(line)).Seq()
        tmp = len(new_contig) - len(line)
        sys.stdout.trace( "Extending", line, "for", tmp)
        line.extendRight(new_contig.suffix(pos = len(line)).Seq(), relevant_als)
        sys.stdout.info( "Extended contig", line, "for", tmp)
        sys.stdout.trace("Reads:")
        sys.stdout.trace( list(line.read_alignments.allInter(line.asSegment().suffix(length=min(len(line), 2000)))))
        sys.stdout.trace( "Sequence:")
        sys.stdout.trace( line.seq)
        return True

    def polyshSegments(self, line, to_polysh):
        # type: (NewLine, Iterable[Segment]) -> List[Segment]
        segs = SegmentStorage()
        corrections = AlignmentStorage()
        line.addListener(segs)
        segs.addAll(to_polysh)
        segs.mergeSegments()
        segs.sort()
        for seg in segs:
            corrections.add(self.polisher.polishSegment(seg, list(line.read_alignments.allInter(seg))))
        line.correctSequence(list(corrections))
        line.removeListener(segs)
        return list(segs)

    def updateCorrectSegments(self, line):
        # type: (NewLine) -> None
        line.updateCorrectSegments(line.asSegment())

    class Record:
        def __init__(self, resolved, next, correct, good_reads, read_bounds):
            # type: (Segment, int, Segment, Set[str], Dict[str, int]) -> None
            self.line = resolved.contig # type: NewLine
            self.resolved = resolved
            self.old_resolved = []
            self.next_resolved_start = next
            self.correct = correct
            self.good_reads = good_reads
            self.read_bounds = read_bounds
            self.reads = [] # type: List[AlignmentPiece]
            self.sorted = True
            self.potential_good = []

        def setResolved(self, seg):
            # type: (Segment) -> None
            if seg.interSize(self.resolved) >= params.k - 1:
                self.resolved = self.resolved.cup(seg)
            else:
                self.old_resolved.append(self.resolved)
                self.resolved = seg
            self.updateGood()

        def add(self, al):
            # type: (AlignmentPiece) -> None
            tmp = list(al.split(100))
            if len(tmp) > 1:
                al = Scorer().polyshAlignment(al, params.alignment_correction_radius)
                for al1 in al.split(100):
                    self.innerAdd(al)
            else:
                self.innerAdd(al)

        def innerAdd(self, al):
            # type: (AlignmentPiece) -> None
            if al.seg_from.left < params.bad_end_length:
                self.potential_good.append(al)
            else:
                self.reads.append(al)
            read = al.seg_from.contig # type: AlignedRead
            if read.id not in self.read_bounds:
                self.read_bounds[read.id] = len(read)
            if al.rc.seg_to.left < 50:
                self.read_bounds[read.id] = min(self.read_bounds[read.id], al.seg_from.right)
            self.sorted = False

        def addAll(self, als):
            # type: (Iterator[AlignmentPiece]) -> None
            for al in als:
                self.add(al)

        def sort(self):
            if not self.sorted:
                self.reads = sorted(self.reads, key = lambda al: -al.seg_to.left)
                self.potential_good = sorted(self.potential_good, key = lambda al: -al.seg_to.left)
                self.sorted = True

        def get(self, num = None, right = None, min_inter = 0):
            # type: (int, Segment, int) -> List[AlignmentPiece]
            self.sort()
            if num is None:
                num = len(self.reads)
            if right is None:
                right = self.resolved.right
            popped = []
            res = []
            while len(res) < num and len(self.reads) > 0 and self.reads[-1].seg_to.left < right:
                al = self.reads.pop()
                necessary_contig_support = min(len(al.seg_from.contig), al.seg_from.left + params.k + 100)
                if al.seg_from.contig.id not in self.good_reads or necessary_contig_support > self.read_bounds[al.seg_from.contig.id]:
                    popped.append(al)
                    if len(al.seg_to) >= min_inter:
                        res.append(al)
            self.reads.extend(popped[::-1])
            return res

        def __iter__(self):
            for al in self.reads[::-1]:
                necessary_contig_support = min(len(al.seg_from.contig), al.seg_from.left + params.k + 100)
                if al.seg_from.contig.id not in self.good_reads or necessary_contig_support > self.read_bounds[al.seg_from.contig.id]:
                    yield al

        def unsupportedAlignments(self, inter_size):
            for al in self.reads[::-1]:
                necessary_contig_support = min(len(al.seg_from.contig), al.seg_from.left + inter_size + 100)
                if al.seg_from.contig.id not in self.good_reads or necessary_contig_support > self.read_bounds[al.seg_from.contig.id]:
                    yield al

        def updateGood(self):
            self.sort()
            while len(self.reads) > 0 and self.reads[-1].seg_to.left <= self.resolved.right - params.k:
                al = self.reads.pop()
                if al.seg_to.interSize(self.resolved) >= params.k:
                    if al.seg_from.contig.id not in self.good_reads:
                        sys.stdout.trace( "New good read:", al)
                        self.good_reads.add(al.seg_from.contig.id)
                else:
                    sys.stdout.trace( "Read does not overlap resolved", al, self.resolved)
            while len(self.potential_good) > 0 and self.potential_good[-1].seg_to.left <= self.resolved.right - params.k:
                al = self.potential_good.pop()
                if al.seg_to.interSize(self.resolved) >= params.k:
                    if al.seg_from.contig.id not in self.good_reads:
                        sys.stdout.trace( "New good read from potential:", al)
                        self.good_reads.add(al.seg_from.contig.id)
                else:
                    sys.stdout.trace( "Read does not overlap resolved", al, self.resolved)

        def pop(self):
            return self.reads.pop()

        def __str__(self):
            return str([self.resolved, self.correct, self.next_resolved_start, self.reads])

    def createRecord(self, resolved, next_start, correct, als, good_reads, read_bounds):
        # type: (Segment, int, Segment, List[AlignmentPiece], Set[str], Dict[str, int]) -> Record
        line = resolved.contig # type: NewLine
        focus = line.segment(resolved.right - params.k, min(correct.right, next_start + params.k))
        res = self.Record(resolved, next_start, correct, good_reads, read_bounds)
        res.addAll(als)
        res.updateGood()
        return res

    def findResolvedBound(self, rec, inter_size):
        # type: (Record, int) -> int
        bad_reads = []
        for read in rec.unsupportedAlignments(inter_size):
            if len(read.seg_to) >= inter_size:
                bad_reads.append(read)
            if len(bad_reads) >= params.min_contra_for_break:
                if bad_reads[-1].seg_to.left - bad_reads[0].seg_to.left > 50:
                    bad_reads = bad_reads[1:]
                else:
                    break
        if len(bad_reads) < params.min_contra_for_break:
            sys.stdout.trace("No resolved bound for", rec.resolved)
            return len(rec.line)
        else:
            sys.stdout.trace( "Resolved bound for", rec.resolved, ":", bad_reads[0].seg_to.left)
            sys.stdout.trace( "Bound caused by read alignments:", map(str, bad_reads))
            return bad_reads[0].seg_to.left

    def attemptProlongResolved(self, rec):
        # type: (Record) -> bool
        sys.stdout.trace( "Working on prolonging", rec.resolved)
        res = self.findAndFilterResolvedBound(rec, params.k)
        if res <= rec.resolved.right:
            sys.stdout.trace( "No luck with", rec.resolved)
            return False
        sys.stdout.trace( "Prolonged", rec.resolved, "to", res)
        rec.setResolved(rec.resolved.contig.segment(rec.resolved.left, res))
        return True

    def findAndFilterResolvedBound(self, rec, sz):
        bound0 = self.findResolvedBound(rec, sz) + params.k * 9 / 10
        bound = min(rec.correct.right, rec.next_resolved_start + sz - 1, bound0)
        res = rec.resolved.right
        if bound > rec.resolved.right:
            sys.stdout.trace( "Checking resolved bound against known copies")
            candidates = self.segmentsWithGoodCopies(rec.resolved, rec.line.segment(max(0, rec.resolved.right - sz), bound), sz)
            sys.stdout.trace( "Candidates:", candidates)
            for candidate in candidates:
                if candidate.left == max(0, rec.resolved.right - sz) and candidate.right > rec.resolved.right:
                    res = candidate.right
        sys.stdout.trace( "Final resolved bound for", rec.resolved, " and k =", sz, ":", res)
        return res

    def attemptJump(self, rec):
        # type: (Record) -> bool
        bound = self.findAndFilterResolvedBound(rec, params.l)
        bad_segments = SegmentStorage()
        for al in rec:
            if al.seg_to.left > bound:
                break
            if al.seg_from.left > min(params.bad_end_length, params.k / 2) and \
                    al.rc.seg_from.left > min(params.bad_end_length, params.k / 2):
                bad_segments.add(al.seg_to)
        for al in self.dot_plot.allInter(rec.line.segment(rec.resolved.right - params.k, bound)):
            if al.seg_from.left > min(params.bad_end_length, params.k / 2):
                if al.rc.seg_from.left > min(params.bad_end_length, params.k / 2):
                    bad_segments.add(al.seg_to)
        bad_segments.mergeSegments(params.k - 200)
        sys.stdout.trace( "Bad segments:", bad_segments)
        good_segments = bad_segments.reverse(rec.line, params.k - 100).reduce(rec.line.segment(rec.resolved.right - params.k, bound))
        for seg in good_segments:
            seg = Segment(seg.contig, max(0, seg.left), seg.right)
            for seg1 in self.segmentsWithGoodCopies(rec.resolved, seg, params.k):
                if len(seg1) >= params.k and seg1.right > rec.resolved.right:
                    rec.setResolved(seg1)
                    return True
        return False

    def segmentsWithGoodCopies(self, resolved, seg, inter_size):
        # type: (Segment, Segment, int) -> List[Segment]
        als = [al for al in self.dot_plot.allInter(seg) if al.seg_from.left > 20 or al.rc.seg_to.left > 20 or al.isIdentical()]
        segs = SegmentStorage()
        for al in als:
            line = al.seg_from.contig # type: NewLine
            if len(al.seg_to) >= inter_size and al.seg_from.right > line.initial[0].seg_to.left:
                cap = al.seg_from.cap(line.suffix(pos=line.initial[0].seg_to.left))
                incorrect = line.correct_segments.reverse(line, inter_size - 1).reduce(cap)
                matching = al.matchingSequence()
                sys.stdout.trace( "Incorrect: ", line, cap, incorrect)
                for seg1 in incorrect:
                    seg2 = matching.mapSegDown(seg.contig, seg1, mapIn=False)
                    sys.stdout.trace( "Relevant unpolished k-mer segment alignment:", seg1, seg2)
                    segs.add(seg2)
                if al.rc.seg_from.left < 50 and al.seg_to.right >= resolved.right - 100:
                    segs.add(al.seg_to.contig.suffix(pos=al.seg_to.right).expand(inter_size + 100))
                    sys.stdout.trace( "Incoming line:", resolved, seg, al)
        segs.mergeSegments(inter_size - 1)
        return list(segs.reverse(seg.contig, inter_size - 1 - max(100, inter_size / 10)).reduce(seg))
