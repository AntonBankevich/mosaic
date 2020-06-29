import itertools
import os
import sys

sys.path.append("py")
sys.path.append(".")
import common.log_params
from common.basic import CreateLog
from typing import Optional, List, Iterable, Tuple
from alignment.align_tools import Aligner, DirDistributor
from common import basic, SeqIO, params
from common.line_align import Scorer
from common.seq_records import NamedSequence
from common.sequences import Consensus, Contig, ContigCollection, Segment, Contig, ContigStorage
from common.alignment_storage import AlignmentPiece, AlignedRead, ReadCollection
from flye_tools.polish import PolishException
from flye_tools.polysh_job import JobPolishing
import flye.polishing.polish

class Polisher:
    def __init__(self, aligner, dir_distributor):
        # type: (Aligner, DirDistributor) -> Polisher
        self.aligner = aligner
        self.dir_distributor = dir_distributor

    def polishMany(self, reads, sequences):
        # type: (Iterable[AlignedRead], List[Contig]) -> List[Contig]
        dir, new_files, same = self.dir_distributor.fillNextDir([(list(sequences), "ref.fasta"), (reads, "reads.fasta")])
        consensus_file_name = new_files[0]
        reads_file_name = new_files[1]
        args = FakePolishingArgs()
        basic.ensure_dir_existance(os.path.join(dir, "work"))
        job = JobPolishing(args, os.path.join(dir, "work"), os.path.join(dir, "log.info"), [reads_file_name], consensus_file_name, "polish")
        polished_file = job.out_files["contigs"]
        if same and not params.clean and os.path.exists(polished_file):
            sys.stdout.trace("Polishing reused:", polished_file)
        else:
            sys.stdout.trace("Running polishing:", polished_file)
            job.run()
        return map(lambda rec: Contig(rec.seq, rec.id), SeqIO.parse_fasta(open(polished_file, "r")))

    def polish(self, reads, consensus):
        # type: (Iterable[NamedSequence], Contig) -> str
        dir, new_files, same = self.dir_distributor.fillNextDir([([consensus], "ref.fasta"), (reads, "reads.fasta")])
        consensus_file_name = new_files[0]
        reads_file_name = new_files[1]
        args = FakePolishingArgs()
        basic.ensure_dir_existance(os.path.join(dir, "work"))
        job = JobPolishing(args, os.path.join(dir, "work"), os.path.join(dir, "log.info"), [reads_file_name], consensus_file_name, "polish")
        polished_file = job.out_files["contigs"]
        if same and not params.clean and os.path.exists(polished_file):
            sys.stdout.trace("Polishing reused:", polished_file)
        else:
            sys.stdout.trace("Running polishing:", polished_file)
            job.run()
        return list(SeqIO.parse_fasta(open(polished_file, "r")))[0].seq

    def polishAndAnalyse(self, reads, polishing_base, reliable_start = None):
        # type: (ReadCollection, Contig, Optional[int]) -> Consensus
        if reliable_start is None:
            reliable_start = len(polishing_base)
        seq = Contig(self.polish(reads, polishing_base), "contig")
        res = [0] * (len(seq) + 1)
        alignment = ReadCollection().extendClean(reads)
        self.aligner.alignReadCollection(alignment, [seq])
        contra = 0
        ok = 0
        late = 0
        for read in alignment:
            for al in read.alignmentsTo(seq.asSegment()):# type: AlignmentPiece
                if al.contradicting(seq.asSegment()):
                    contra += 1
                elif al.seg_to.left > reliable_start:
                    late += 1
                else:
                    res[al.seg_to.left] += 1
                    res[al.seg_to.right] -= 1
                    ok += 1
        for i in range(1, len(res)):
            res[i] += res[i - 1]
        sys.stdout.trace("Polyshed and analysed using", len(alignment), "reads. Ok:", ok, "late:", late, "contra:", contra)
        # if contra > 10 or contra > ok / 2:
        #     for read in alignment:
        #         print read
        #         for al in read.alignmentsTo(seq.asSegment()):
        #             if al.contradictingRTC(seq.asSegment()):
        #                 print "contra_al:", al
        #             elif al.seg_to.left > reliable_start:
        #                 print "late_al:", al
        #             else:
        #                 print "ok_al:", al
        return Consensus(seq.seq, res)

    # Returns alignment of polished sequence to old sequence
    def polishSegment(self, seg, als):
        # type: (Segment, List[AlignmentPiece]) -> AlignmentPiece
        sys.stdout.trace("Polishing segment", seg)
        w = max(900, params.k)
        r = 50
        first = seg.left / w
        last = min(seg.right + w - 1, len(seg.contig)) / w
        segs = []
        for i in range(first, last + 1):
            segs.append(Segment(seg.contig, max(0, i * w - r), min(len(seg.contig), (i + 1) * w + r)))
        als_by_segment = [[] for i in range(last - first + 1)]
        for al in als:
            l = al.seg_to.left / w
            r = al.seg_to.right / w + 1
            for i in range(max(0, l - first - 1), min(last - first + 1, r - first + 2)):
                if al.seg_to.inter(segs[i]):
                    als_by_segment[i].append(al)
        res_als = []
        for seg1, seg_als in zip(segs, als_by_segment):
            if seg1.inter(seg):
                res_als.append(self.polishSmallSegment(seg1, seg_als))
        res = AlignmentPiece.GlueOverlappingAlignments(res_als)
        return res.reduce(target=seg)

    # Returns alignment of polished sequence to old sequence
    def polishSmallSegment(self, seg, als):
        # type: (Segment, List[AlignmentPiece]) -> AlignmentPiece
        ok = False
        for al in als:
            if al.seg_to.contains(seg):
                ok = True
        if not ok:
            sys.stdout.log(common.log_params.LogPriority.warning, "Warning", seg, "has no covering reads")
            return AlignmentPiece.Identical(seg.asContig().asSegment(), seg)
        reads = []
        start = basic.randomSequence(200)
        end = basic.randomSequence(200)
        for al in als:
            new_seq = ""
            al = al.reduce(target=seg)
            if al.seg_to.left < seg.left + 20:
                new_seq += start
            new_seq += al.seg_from.Seq()
            if al.seg_to.right > seg.right - 20:
                new_seq += end
            reads.append(NamedSequence(new_seq, al.seg_from.contig.id))
        base = Contig(start + seg.Seq() + end, "base")
        polished = None
        try:
            polished = Contig(self.polish(reads, base), "polished")
        except PolishException:
            sys.stdout.log(common.log_params.LogPriority.warning, "Warning", seg, "has a sequence very different from reads. Using reads to correct.")
            for al, read in zip(als, reads):
                if al.seg_to.contains(seg):
                    try:
                        polished = Contig(self.polish(reads, Contig(read.seq, read.id)), "polished")
                        break
                    except PolishException:
                        pass
        if polished is None:
            sys.stdout.log(common.log_params.LogPriority.warning, "Warning", seg, "could not be corrected even though some reads cover it.")
            polished = seg.asContig()
        als = list(self.aligner.overlapAlign([polished], ContigStorage([base])))
        for al in als:
            if al.seg_from.left < 10 and al.rc.seg_from.left < 10:
                mapping = AlignmentPiece.Identical(base.segment(len(start), len(base) - len(end)), seg)
                return al.compose(mapping)
        assert False, "No alignment from polished to base: " + str(als)

    def allEnds(self, als, min_cov = 5):
        contig = als[0].seg_to.contig
        res_contigs = []
        res_als = []
        undefined = list(als)
        while True:
            new_contig, new_als = self.polishEnd(undefined, min_cov, 0)
            tmp_als = []
            read_ids = set()
            for al in new_als:
                if al.seg_to.right > len(contig) + 500:
                    tmp_als.append(al)
                    read_ids.add(al.seg_from.contig.id)
            if len(tmp_als) == 0:
                break
            if len(tmp_als) >= 5:
                res_contigs.append(new_contig)
                res_als.append(list(tmp_als))
            tmp_als = []
            for al in undefined:
                if al.seg_from.contig.id not in read_ids:
                    tmp_als.append(al)
            undefined = tmp_als
        for contig, tmp_als in zip(res_contigs, res_als):
            al = self.polishSmallSegment(contig.asSegment(), tmp_als)
            yield al.seg_to.contig, [al1.compose(al) for al1 in tmp_als]

    def polishEnd(self, als, min_cov = 4, min_cov_frac = 0.8, max_extension = None):
        # type: (List[AlignmentPiece], int, int, int) -> Tuple[Contig, List[AlignmentPiece]]
        if max_extension is None:
            max_extension = 10000000000
        scorer = Scorer()
        contig = als[0].seg_to.contig
        max_len = max_extension + len(contig)
        sys.stdout.trace("Polishing end of", als[0].seg_to.contig)
        new_contig = contig.asSegment().asContig()
        relevant_als = [al.changeTargetContig(new_contig) for al in als if al.rc.seg_to.left < 100]
        finished_als = []
        while True:
            tmp = []
            for al in relevant_als:
                if al.seg_to.inter(new_contig.asSegment().suffix(length=100)) and al.rc.seg_from.left > 100:
                    tmp.append(al)
                else:
                    finished_als.append(al)
            relevant_als = tmp
            if len(relevant_als) < min_cov:
                break
            start = "ACGTTCGA" + basic.randomSequence(params.flanking_size) + new_contig.asSegment().suffix(length=min(params.flanking_size, len(new_contig))).Seq()
            reduced_read_list = [
                AlignedRead.new(start + al.seg_from.contig.asSegment().suffix(pos=al.seg_from.right).Seq(), str(i) + "_" + al.seg_from.contig.id)
                for i, al in enumerate(relevant_als)]
            reduced_reads = ReadCollection(reduced_read_list)
            found = False
            for base_al in relevant_als:
                if base_al.rc.seg_from.left < params.flanking_size:
                    continue
                # Base consists 500 random nucleotides and 500 last nucls from the polished sequence a segment of read of length at most 500
                base_segment = base_al.seg_from.contig.segment(base_al.seg_from.right,
                                                     min(len(base_al.seg_from.contig), base_al.seg_from.right + max(params.window_size, params.k)))
                base = Contig(start + base_segment.Seq(), "base")
                for read in reduced_read_list:
                    read.clean()
                polished_base = Contig(self.polish(reduced_reads, base), "polished_base")
                for al in self.aligner.localAlign(reduced_reads, ContigStorage().addAll([polished_base])):
                    reduced_reads.reads[al.seg_from.contig.id].addAlignment(al)
                candidate_alignments = []
                for read in reduced_read_list:
                    candidate_alignments.append(None)
                    for al in read.alignmentsTo(polished_base.asSegment()):
                        if al.seg_to.left == 0 and ((candidate_alignments[-1] is None or candidate_alignments[-1].seg_to.right < al.seg_to.right)):
                            candidate_alignments[-1] = al
                trimmedAlignments = []
                for i, al in enumerate(candidate_alignments):
                    assert al is not None, reduced_read_list[i]
                    trimmedAlignments.append(al.trimByQuality(0.4, 100))
                contra_index = 0
                contra = []
                support = len(trimmedAlignments)
                cutoff_pos = len(start)
                for al in sorted(trimmedAlignments, key = lambda al: al.seg_to.right):
                    while contra_index < len(contra) and contra[contra_index].seg_to.right < al.seg_to.right - 50:
                        contra_index += 1
                    if support >= min_cov and len(contra) - contra_index <= (1 - min_cov_frac) * support:
                        cutoff_pos = al.seg_to.right
                        support -= 1
                        if al.contradictingRTCRight():
                            contra.append(al)
                    else:
                        sys.stdout.trace("Stopped at:", support, contra_index, (1 - min_cov_frac) * support)
                        break
                sys.stdout.trace("Positions:", [al.seg_to.right for al in trimmedAlignments])
                sys.stdout.trace("Contra:", contra)
                if cutoff_pos > len(start) + 100:
                    sys.stdout.trace("Chose to use read", base_al.__repr__(), "Extended for", cutoff_pos - len(start), "Alignments:")
                    sys.stdout.trace(map(str, reduced_read_list))
                    found = True
                    new_contig_candidate = Contig(new_contig.seq + polished_base[len(start):cutoff_pos], "candidate")
                    embedding = AlignmentPiece.Identical(polished_base.segment(len(start), cutoff_pos), new_contig_candidate.asSegment().suffix(pos=len(new_contig)))
                    read_mappings = []
                    for al1, al2 in zip(candidate_alignments, relevant_als):
                        seg_from = al2.seg_from.contig.asSegment().suffix(length = len(al1.seg_from.contig) - len(start))
                        seg_to = al1.seg_from.contig.asSegment().suffix(length = len(al1.seg_from.contig) - len(start))
                        read_mappings.append(AlignmentPiece.Identical(seg_from, seg_to))
                    embedded_alignments = []
                    for al1, al2 in zip(candidate_alignments, read_mappings):
                        if al1.seg_to.right <= len(start) + 10:
                            embedded_alignments.append(None)
                        else:
                            tmp = al2.compose(al1)
                            if tmp.seg_to.left > embedding.seg_from.right - 10:
                                embedded_alignments.append(None)
                            else:
                                embedded_alignments.append(tmp.compose(embedding))
                    corrected_relevant_alignments = [al.targetAsSegment(new_contig_candidate.asSegment().prefix(len(new_contig))) for al in relevant_als]
                    relevant_als = []
                    for al1, al2 in zip(corrected_relevant_alignments, embedded_alignments):
                        if al2 is None:
                            al = al1
                        else:
                            al = al1.mergeDistant(al2)
                            if al is None:
                                al = al1
                            elif al1.seg_from.dist(al2.seg_from) >= 10 or al1.seg_to.dist(al2.seg_to) >= 10:
                                al = scorer.polyshAlignment(al, params.alignment_correction_radius)
                        relevant_als.append(al)
                    finished_als = [al.targetAsSegment(new_contig_candidate.asSegment().prefix(len(new_contig))) for al in finished_als]
                    new_contig = new_contig_candidate
                    break
                else:
                    sys.stdout.trace("Could not prolong with read", base_al, "Alignments:")
                    sys.stdout.trace(map(str, reduced_read_list))
            if len(new_contig) >= max_len:
                break
            if not found:
                break
        return new_contig, relevant_als + finished_als

class FakePolishingArgs:
    def __init__(self):
        self.num_iters = params.num_iters
        self.platform = params.technology
        self.threads = params.threads

if __name__ == "__main__":
    reads_file = sys.argv[2]
    consensus_file = sys.argv[3]
    dir = sys.argv[1]
    extra_params = sys.argv[4:]
    CreateLog(dir)
    dd = DirDistributor(dir)
    aligner = Aligner(dd)
    polisher = Polisher(aligner, dd)
    reads = ContigStorage().loadFromFasta(open(reads_file, "r"), num_names=False)
    ref = ContigStorage().loadFromFasta(open(consensus_file, "r"), num_names=False)
    if "accurate" in extra_params:
        res = []
        als = sorted(aligner.overlapAlign(reads, ref), key = lambda al: al.seg_to.contig.id)
        for rid, rals in itertools.groupby(als, key = lambda al: al.seg_to.contig.id):
            if basic.isCanonocal(rid):
                contig = ref[rid]
                corrected_seq = polisher.polishSegment(contig.asSegment(), list(rals)).seg_from.Seq()
                res.append(Contig(corrected_seq, rid))
    else:
        res = polisher.polishMany(reads, list(ref.unique()))
    res_file = os.path.join(dir, "res.fasta")
    rf = open(res_file, "w")
    for c in res:
        SeqIO.write(c, rf, "fasta")
    rf.close()
    aligner.align_files(res_file, [reads_file], 16, "pacbio", "overlap", os.path.join(dir, "res.sam"))
