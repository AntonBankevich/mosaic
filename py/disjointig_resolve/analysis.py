import itertools

from typing import List

from alignment.align_tools import Aligner
from common import basic, params
from common.alignment_storage import ReadCollection
from common.sequences import Segment, ContigStorage
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.line_storage import NewLineStorage


class CoverageAnalyser:
    def __init__(self, aligner, reads, seg_size = params.ss_for_kl_adjustment):
        # type: (Aligner, ReadCollection, int) -> None
        self.aligner = aligner
        self.reads = reads
        self.seg_size = seg_size
        self.recs = None #type: List[CoverageAnalyser.CoverageRecord]

    class CoverageRecord:
        def __init__(self, k, covs):
            self.k = k
            self.covs = covs
            for i in range(0, len(covs) - 1):
                self.covs[i + 1] += self.covs[i]

        def __str__(self):
            if len(self.covs) > 0 and self.covs[-1] > 0:
                return "k:" + str(self.k) + ": " + " ".join(map(lambda cov: "%0.3f" % (float(cov) / self.covs[-1]), self.covs))
            else:
                return "None"

    def covSegments(self, region, segs, k):
        tmp = []
        for seg in segs:
            if seg.interSize(region) >= k:
                seg = seg.cap(region)
                tmp.append((seg.left, -1))
                tmp.append((seg.right - k + 1, 1))
        tmp = sorted(tmp)
        cur_cov = 0
        prev = 1000
        for pos, diff in tmp:
            assert pos >= prev
            if pos > prev:
                yield cur_cov, pos - prev
            # covs[i][min(cur_cov, len(covs[i]) - 1)] += pos - prev
            cur_cov -= diff
            prev = pos
        assert cur_cov == 0
        if prev < region.right - k + 1:
            yield 0, region.right - k + 1 - prev

    def analyseSegments(self, segs):
        # type: (List[Segment]) -> None
        contigs = ContigStorage()
        contigs.addAll([seg.asContig() for seg in segs if len(seg) > 5000])
        res = [] # type: List[Segment]
        for al in self.aligner.overlapAlign(self.reads, contigs):
            if basic.isCanonocal(al.seg_to.contig.id):
                res.append(al.seg_to)
            else:
                res.append(al.seg_to.RC())
        res = sorted(res, key=lambda seg: (seg.contig.id, seg.left))
        covs = [[0] * params.maxCoverageThreshold for i in range(100)]
        for contig, it in itertools.groupby(res, key = lambda seg: seg.contig):
            segs = list(it)
            shrink = contig.asSegment().shrink(1000)
            bad_seg = False
            for cov, slen in self.covSegments(shrink, segs, 1):
                if cov < 3:
                    bad_seg = True
            if bad_seg:
                continue
            for i in range(len(covs)):
                k = 500 + i * 100
                for cov, slen in self.covSegments(shrink, segs, k):
                    covs[i][min(cov, len(covs[i]) - 1)] += slen
        self.recs = [CoverageAnalyser.CoverageRecord(500 + i * 100, covs[i]) for i in range(len(covs)) if covs[i] > 1000]

    def analyseSegmentCoverage(self, contigs):
        # type: (ContigStorage) -> None
        segs = []
        for contig in contigs.unique():
            if len(contig) < 3 * self.seg_size:
                segs.append(contig.prefix(self.seg_size))
                segs.append(contig.suffix(pos=len(contig) - self.seg_size))
            else:
                segs.append(contig.asSegment())
        self.analyseSegments(segs)


    def printAnalysis(self):
        for cov in self.recs:
            print cov

    def chooseK(self, threshold):
        assert threshold < params.maxCoverageThreshold
        res = params.k
        for rec in self.recs:
            if rec.covs[-1] > 0 and rec.covs[threshold] < params.uncoveredFractionForK * rec.covs[-1]:
                res = rec.k
            else:
                break
        return res

