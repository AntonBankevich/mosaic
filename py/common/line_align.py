from collections import deque

from typing import Optional, Tuple, List, Dict

from common.scoring_model import ScoringModel
from common.sequences import Segment

import sys
if __name__ == "__main__":
    sys.path.append("py")
from common.alignment_storage import AlignmentPiece, MatchingSequence
from common import params

class EventCounter:
    def __init__(self):
        self.reset()

    def reset(self):
        self.ms = 0
        self.mm = 0
        self.indel = 0
        self.homo = 0

    def __str__(self):
        return " ".join(map(str, [self.ms, self.mm, self.indel, self.homo]))

events = EventCounter()

class Scorer:

    def accurateScore(self, alignment, radius): #This score is nonsymmetric!!! Insertions and deletions have different cost
        # type: (MatchingSequence, int) -> int
        prev = Storage(alignment[0][1], alignment[1][1] + radius, self.scores.inf)
        prev.set(alignment[0][1], 0)
        cur_del = 0
        for j in range(alignment[0][1] + 1, min(alignment[1][1] + radius + 1, len(alignment.seq_to))):
            cur_del += self.scores.scoreDel(alignment.seq_to[j - 1])
            prev.set(j, cur_del)
        cur = 0
        bounds = self.generateBounds(alignment, params.alignment_smoothing_radius)
        total = []
        for i in range(alignment[0][0] + 1, alignment[-1][0] + 1):
            j_min = max(min(bounds[i - alignment[0][0]][0], alignment[cur][1]) - radius, alignment[0][1])
            if alignment[cur + 1][0] == i and cur + 2 < len(alignment):
                cur += 1
                if alignment[cur + 1][1] - alignment[cur][1] > 10 and alignment[cur + 1][0] - alignment[cur][0] > 10:
                    sys.stdout.trace("Long gap:", alignment[cur], alignment[cur + 1])
            j_max = min(max(bounds[i - alignment[0][0]][1], alignment[cur + 1][1]) + radius, alignment[-1][1])
            total.append(j_max - j_min)
            ne = Storage(j_min, j_max + 1, self.scores.inf)
            c1 = alignment.seq_from[i]
            for j in range(j_min, j_max + 1):
                res = self.scores.inf
                c2 = alignment.seq_to[j]
                if c1 == c2:
                    res = min(res, prev.get(j - 1) + self.scores.scoreMatch(alignment.seq_from[i - 1]))
                    if i > 0 and j > 0 and alignment.seq_from[i - 1] == c1 and alignment.seq_to[j - 1] == c2:
                        if i > 1 and alignment.seq_from[i - 2] == alignment.seq_from[i - 1]:
                            res = min(res, prev.get(j) + self.scores.scoreHomo(c1))
                        if j > 1 and alignment.seq_to[j - 2] == alignment.seq_to[j - 1]:
                            res = min(res, ne.get(j - 1) + self.scores.scoreHomo(c2))
                else:
                    res = min(res, prev.get(j - 1) + self.scores.scoreMM(c1, c2))
                res = min(res, prev.get(j) + self.scores.scoreIns(c1))
                res = min(res, ne.get(j - 1) + self.scores.scoreDel(c2))
                ne.set(j, res)
            prev = ne
        return prev.get(alignment[-1][1])

    def __init__(self, scores = None):
        if scores is None:
            scores = params.scores
        self.scores = scores # type: ScoringModel

    def countHomo(self, seq, pos, step):
        cpos = pos + step
        while 0 <= cpos < len(seq):
            if seq[cpos] != seq[pos]:
                break
            cpos += step
        return (cpos - pos) * step
    #     # type: (MatchingSequence) -> int
    #     matches = alignment.matches
    #     res = 0
    #     for match1, match2 in zip(matches[:-1], matches[1:]):
    #         if match2[0] - match1[0] == 1 and match2[1] - match1[1] == 1:
    #             continue
    #         l = [match2[0] - match1[0] - 1, match2[1] - match1[1] - 1]
    #         homo = 0
    #         if l[1] > l[0]:
    #             homo += self.countHomo(alignment.seq_to, match1[1], 1) - 1
    #             homo += self.countHomo(alignment.seq_to, match2[1], -1) - 1
    #             homo = min(homo, (match2[1] - match1[1]) - (match2[0] - match1[0]))
    #             res += self.scores.scoreMM()csub_score * l[0] + homo * self.scores.homo_score + (l[1] - l[0] - homo) * self.scores.del_score
    #         else:
    #             homo += self.countHomo(alignment.seq_from, match1[0], 1) - 1
    #             homo += self.countHomo(alignment.seq_from, match2[0], -1) - 1
    #             homo = min(homo, l[0] - l[1])
    #             res += self.scores.sub_score * l[1] + homo * self.scores.homo_score + (l[0] - l[1] - homo) * self.scores.del_score
    #     return res

    def score(self, alignment):
        assert False

    def countEvents(self, alignment):
        # type: (MatchingSequence) -> Tuple[int, int, int, int]
        # Matches, mismatches, indels, homo
        matches = alignment.matches
        res = 0
        ms = 0
        mm = 0
        indels = 0
        homo = 0
        for match1, match2 in zip(matches[:-1], matches[1:]):
            assert match2[0] - match1[0] == 1 or match2[1] - match1[1] == 1
            if match2[0] - match1[0] == 1 and match2[1] - match1[1] == 1:
                continue
            if match2[0] - match1[0] == 1:
                t = 1
                seq = alignment.seq_to
            else:
                t = 0
                seq = alignment.seq_from
            tmp = match1[t] + 1
            while tmp < match2[t] and seq[match1[t]] == seq[tmp]:
                tmp += 1
                homo += 1
            tmp2 = match2[t] - 1
            while tmp2 > tmp and seq[match2[t]] == seq[tmp2]:
                tmp2 -= 1
                homo += 1
            indels += match2[0] - match1[0] + match2[1] - match1[1]  - 2
        for m in alignment.matches:
            if alignment.seq_from[m[0]] == alignment.seq_to[m[1]]:
                ms += 1
            else:
                mm += 1
        return ms, mm, indels, homo

    def maxInRange(self, vals, r):
        dec = deque()
        res = []
        curl = 0
        curr = 1
        dec.append(vals[0])
        for i in range(vals[0][1], vals[-1][1] + 1):
            while curr + 1 < len(vals) and i + r >= vals[curr + 1][1]:
                curr += 1
                while len(dec) > 0 and dec[-1][0] <= vals[curr][0]:
                    dec.pop()
                dec.append(vals[curr])
            while curl < len(vals) and len(dec) > 0 and i > dec[0][1]:
                dec.popleft()
                curl += 1
            if len(dec) > 0:
                res.append(dec[0][0])
            else:
                res.append(None)
        assert res[0] is not None and res[-1] is not None
        res1 = list(res)
        for i, val in enumerate(res):
            if val is None:
                res1[i] = res1[i - 1]
        for i in range(len(res)):
            if res[len(res) - 1 - i] is None:
                res1[len(res) - 1 - i] = max(res1[len(res) - 1 - i], res1[len(res) - i])
        res = res1
        for i in range(len(res))[::-1]:
            res[i] = max(res[i], res[max(0, i - r)])
        return res

    def minInRange(self, vals, r):
        res = self.maxInRange([(-a, b) for a, b in vals], r)
        return [-a for a in res]

    def generateBounds(self, alignment, r):
        # type: (MatchingSequence, int) -> List[Tuple[int, int]]
        upper = self.maxInRange([(b - a, a) for a, b in alignment.matches], r)
        lower = self.minInRange([(b - a, a) for a, b in alignment.matches], r)
        return [(a + i, b + i) for a, b, i in zip(lower, upper, range(alignment.matches[0][0], alignment.matches[-1][0] + 1))]

    def adjustMin(self, old_val, old_shift, new_val, new_shift):
        # type: (int, Tuple[int, int], int, Tuple[int, int]) -> Tuple[int, Tuple[int, int]]
        if new_val < old_val:
            return new_val, new_shift
        else:
            return old_val, old_shift

    def polyshAlignment(self, alignment, radius):
        # type: (AlignmentPiece, int) -> AlignmentPiece
        return self.polyshMatching(alignment.matchingSequence(), radius).asAlignmentPiece(alignment.seg_from.contig, alignment.seg_to.contig)

    def polyshMatching(self, alignment, radius):
        # type: (MatchingSequence, int) -> MatchingSequence
        assert  len(alignment) > 0
        if len(alignment) == 1:
            return alignment
        storage = RectStorage(alignment[0][0], alignment[-1][0])
        prev = Storage(alignment[0][1], alignment[1][1] + radius, self.scores.inf)
        storage.set(alignment[0][0], Storage(alignment[0][1], alignment[1][1] + radius))
        cur_del = 0
        prev.set(alignment[0][1], 0)
        bounds = self.generateBounds(alignment, params.alignment_smoothing_radius)
        for j in range(alignment[0][1] + 1, min(alignment[1][1] + radius + 1, len(alignment.seq_to))):
            cur_del += self.scores.scoreDel(alignment.seq_to[j - 1])
            prev.set(j, cur_del)
            storage.get(alignment[0][0]).set(j, (alignment[0][0], j - 1))
        cur = 0
        for i in range(alignment[0][0] + 1, alignment[-1][0] + 1):
            j_min = max(min(bounds[i - alignment[0][0]][0], alignment[cur][1]) - radius, alignment[0][1])
            if alignment[cur + 1][0] == i and cur + 2 < len(alignment):
                cur += 1
                if alignment[cur + 1][1] - alignment[cur][1] > 10 and alignment[cur + 1][0] - alignment[cur][0] > 10:
                    sys.stdout.trace("Long gap:", alignment[cur], alignment[cur + 1])
            j_max = min(max(bounds[i - alignment[0][0]][1], alignment[cur + 1][1]) + radius, alignment[-1][1])
            ne = Storage(j_min, j_max + 1, self.scores.inf)
            storage.set(i, Storage(j_min, j_max + 1, None))
            c1 = alignment.seq_from[i]
            for j in range(j_min, j_max + 1):
                c2 = alignment.seq_to[j]
                res = self.scores.inf
                res_shift = (0, 0)
                if c1 == c2:
                    res, res_shift = self.adjustMin(res, res_shift, prev.get(j - 1), (-1, -1))
                    if i > 0 and j > 0 and alignment.seq_from[i - 1] == alignment.seq_from[i] and alignment.seq_to[j - 1] == alignment.seq_to[j]:
                        if i > 1 and alignment.seq_from[i - 2] == alignment.seq_from[i - 1]:
                            res, res_shift = self.adjustMin(res, res_shift, prev.get(j) + self.scores.scoreHomo(c1), (-1, 0))
                        if j > 1 and alignment.seq_to[j - 2] == alignment.seq_to[j - 1]:
                            res, res_shift = self.adjustMin(res, res_shift, ne.get(j - 1) + self.scores.scoreHomo(c1), (0, -1))
                else:
                    res, res_shift = self.adjustMin(res, res_shift, prev.get(j - 1) + self.scores.scoreMM(c1, c2), (-1, -1))
                res, res_shift = self.adjustMin(res, res_shift, prev.get(j) + self.scores.scoreIns(c1), (-1, 0))
                res, res_shift = self.adjustMin(res, res_shift, ne.get(j - 1) + self.scores.scoreDel(c2), (0, -1))
                ne.set(j, res)
                storage.get(i).set(j, (i + res_shift[0], j + res_shift[1]))
            prev = ne
        cur_i = alignment[-1][0]
        cur_j = alignment[-1][1]
        matches = [(cur_i, cur_j)]
        prev_i = cur_i
        prev_j = cur_j
        tmp_i = prev_i
        tmp_j = prev_j
        events.reset()
        while cur_i != alignment[0][0] or cur_j != alignment[0][1]:
            cur_i, cur_j = storage.get(cur_i).get(cur_j)
            assert tmp_i - cur_i <= 1 and tmp_j - cur_j <= 1
            if tmp_i - cur_i == 1 and tmp_j - cur_j == 1:
                if alignment.seq_from[cur_i] == alignment.seq_to[cur_j]:
                    events.ms += 1
                else:
                    events.mm += 1
            else:
                if alignment.seq_from[cur_i] == alignment.seq_to[cur_j] and \
                        alignment.seq_from[cur_i] == alignment.seq_to[tmp_j] and \
                        alignment.seq_from[tmp_i] == alignment.seq_to[cur_j]:
                    events.homo += 1
                else:
                    events.indel += 1
            tmp_i = cur_i
            tmp_j = cur_j
            if cur_i != prev_i and cur_j != prev_j and alignment.seq_from[cur_i] == alignment.seq_to[cur_j]:
                matches.append((cur_i, cur_j))
                prev_i = cur_i
                prev_j = cur_j
        matches[-1] = (alignment[0][0], alignment[0][1])
        return MatchingSequence(alignment.seq_from, alignment.seq_to, matches[::-1])

    def oldCompare(self, piece1, piece2):
        # type: (AlignmentPiece, AlignmentPiece) -> Tuple[Optional[int], Optional[int], Optional[int]]
        pid1 = piece1.percentIdentity()
        pid2 = piece2.percentIdentity()
        if pid1 < params.min_expected_Pacbio_PI and pid2 < params.min_expected_Pacbio_PI:
            return None, None, None
        # we only consider contradictions with the right part of the line since reads were aligned to line center suffixes
        contra1 = piece1.contradicting(piece1.seg_to.contig.centerPos.suffix())
        contra2 = piece2.contradicting(piece2.seg_to.contig.centerPos.suffix())
        if pid1 < params.min_allowed_Pacbio_PI or (contra1 and piece2.seg_from.right > piece1.seg_from.right + 500):
            return None, self.score(piece2.matchingSequence()), None
        if pid2 < params.min_allowed_Pacbio_PI or (contra2 and piece1.seg_from.right > piece2.seg_from.right + 500):
            return self.score(piece1.matchingSequence()), None, None
        return self.scoreCommon(piece1, piece2)

    def scoreInCorrectSegments(self, al1, seg1, al2, seg2):
        # type: (AlignmentPiece, Segment, AlignmentPiece, Segment) -> Tuple[Optional[int], Optional[int], Optional[int]]
        invalid1 = al1.contradictingRTC(tail_size=params.bad_end_length)
        invalid2 = al2.contradictingRTC(tail_size=params.bad_end_length)
        if invalid1 and invalid2:
            return 0, 0, 0
        elif invalid1:
            return 1000, 0, 1000
        elif invalid2:
            return 0, 1000, 1000
        p1 = 0
        p2 = 0
        alignment_length_penalty = min(self.scores.avgInsScore(), self.scores.avgDelScore()) / 3
        # we penalize alignment that ends earlier on the left by the alignment length differenth but only up to the start of alignment target
        # if alignment is too close to contig end it is ignored.
        if al1.seg_from.left > al2.seg_from.left:
            if al1.seg_to.left > 50:
                p1 += min(al1.seg_from.left - al2.seg_from.left, al1.seg_to.left) * alignment_length_penalty
        else:
            if al2.seg_to.left > 50:
                p2 += min(al2.seg_from.left - al1.seg_from.left, al2.seg_to.left) * alignment_length_penalty
        # same for the right end
        if al1.seg_from.right > al2.seg_from.right:
            if al1.rc.seg_to.left > 50:
                p2 += min(al1.seg_from.right - al2.seg_from.right, len(al2.seg_to.contig) - al2.seg_to.right) * alignment_length_penalty
        else:
            if al2.rc.seg_to.left > 50:
                p1 += min(al2.seg_from.right - al1.seg_from.right, len(al1.seg_to.contig) - al1.seg_to.right) * alignment_length_penalty
        full_scores = self.scoreCommon(al1, al2)
        # On this segment both alignments go to correct sequences. We place larger weight on differences in this segment.
        q1 = al1.reduce(target=seg1).seg_from
        q2 = al2.reduce(target=seg2).seg_from
        if q1.inter(q2):
            seg = q1.cap(q2)
            correct_scores = self.scoreCommon(al1.reduce(query=seg), al2.reduce(query=seg))
        else:
            correct_scores = full_scores
        if correct_scores[0] > correct_scores[1]:
            p = p1 - p2
        else:
            p = p2 - p1
        s1 = (full_scores[0] + correct_scores[0] * 9) / 10 + p1
        s2 = (full_scores[1] + correct_scores[1] * 9) / 10 + p2
        s12 = (full_scores[2] + correct_scores[2] * 9) / 10 + abs(p)
        sys.stdout.trace((p1, p2), full_scores, correct_scores, (s1, s2, s12))
        return s1, s2, s12

    def scoreCommon(self, piece1, piece2):
        # type: (AlignmentPiece, AlignmentPiece) -> Tuple[int, int, int]
        if not piece1.seg_from.inter(piece2.seg_from):
            sys.stdout.trace("Nonintersecting alignments fighting", piece1, piece2)
            return 0, 0, 0
        matches1, matches2 = self.cutHomo(piece1.matchingSequence(), piece2.matchingSequence())
        if matches1 is None:
            sys.stdout.info(piece1, piece2)
            assert False
        composite = matches1.composeDifference(matches2)
        matches1 = matches1.reduceTarget(composite.matches[0][0], composite.matches[-1][0] + 1)
        matches2 = matches2.reduceTarget(composite.matches[0][1], composite.matches[-1][1] + 1)
        accurate1 = self.accurateScore(matches1, params.score_counting_radius)
        accurate2 = self.accurateScore(matches2, params.score_counting_radius)
        if accurate1 > accurate2:
            composite = composite.reverse()
        accurate12 = self.accurateScore(composite, params.score_counting_radius)
        if not (abs(accurate1 - accurate2) <= accurate12 <= accurate1 + accurate2):
            sys.stdout.trace("Triangle inequality failed: " + \
                  str(accurate1) + " " + str(accurate2) + " " + \
                  str(abs(accurate1 - accurate2)) + "<=" + str(accurate12) + "<=" + str(accurate1 + accurate2))
        return accurate1, accurate2, self.accurateScore(composite, params.score_counting_radius)

    def cutHomo(self, m1, m2):
        # type: (MatchingSequence, MatchingSequence) -> Tuple[MatchingSequence, MatchingSequence]
        if len(list(m1.common_from(m2))) == 0:
            return None, None
        for a, b in m1.common_from(m2):
            p1 = m1[a][1]
            p2 = m2[b][1]
            p = m1[a][0]
            if m1.seq_to[p1 + 1] != m1.seq_to[p1] and m2.seq_to[p2 + 1] != m2.seq_to[p2] and m1.seq_from[p] != m1.seq_from[p + 1]:
                m1.matches = m1.matches[a:]
                m2.matches = m2.matches[b:]
                break
        for a, b in list(m1.common_from(m2))[::-1]:
            p1 = m1[a][1]
            p2 = m2[b][1]
            p = m1[a][0]
            if m1.seq_to[p1 - 1] != m1.seq_to[p1] and m2.seq_to[p2 - 1] != m2.seq_to[p2] and m1.seq_from[p] != m1.seq_from[p - 1]:
                m1.matches = m1.matches[:a + 1]
                m2.matches = m2.matches[:b + 1]
                break
        return m1, m2


class Storage:
    def __init__(self, left, right, default = None):
        self.left = left
        self.right = right
        self.default = default
        self.vals = []
        for i in range(left, right + 1):
            self.vals.append(default)

    def set(self, a, val):
        assert self.left <= a <= self.right
        self.vals[a - self.left] = val

    def get(self, a):
        if self.left <= a <= self.right:
            return self.vals[a - self.left]
        else:
            return self.default

class RectStorage(Storage):
    def __init__(self, left, right, default = None):
        Storage.__init__(self, left, right, default)
        self.vals = self.vals # type: List[Storage]

class Tournament:
    def __init__(self):
        self.scorer = Scorer()

    def fight(self, c1, c2):
        # type: (AlignmentPiece, AlignmentPiece) -> Optional[AlignmentPiece]
        assert c1.seg_from.contig == c2.seg_from.contig
        s1, s2, s12 = self.scorer.oldCompare(c1, c2)
        if s12 is None:
            if s1 is None and s2 is not None:
                sys.stdout.trace("Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", c2)
                return c2
            elif s1 is not None and s2 is None:
                sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "Winner:", c1)
                return c1
            elif s1 is None and s2 is None:
                sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", None, s12, s1, s2, "No winner")
                return None
            assert False, "Strange comparison results"
        else:
            if s12 < 25 or (s12 < 100 and abs(s1 - s2) < s12 * 0.8) or abs(s1 - s2) < s12 * 0.65:
                sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "No winner")
                return None
            if s1 > s2:
                sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "Winner:", c2)
                return c2
            else:
                sys.stdout.trace( "Fight:", c1, c2, "Comparison results:", abs(s1 - s2), s12, s1, s2, "Winner:", c1)
                return c1

    def tournament(self, candidates):
        # type: (list[AlignmentPiece]) -> Optional[AlignmentPiece]
        best = None
        for candidate in candidates:
            if best is None:
                best = candidate
            else:
                best = self.fight(candidate, best)
        if best is None:
            return None
        if len(candidates) > 2:
            for candidate in candidates:
                if candidate == best:
                    continue
                fight_results = self.fight(candidate, best)
                if fight_results is None or fight_results != best:
                    return None
        return best


if __name__ == "__main__":
    tmp = MatchingSequence(sys.argv[1], sys.argv[2], [(0, 0), (len(sys.argv[1]) - 1, len(sys.argv[2]) - 1)])
    print Scorer().accurateScore(tmp, params.alignment_correction_radius)