import itertools
import random
import sys
import inspect
import traceback
from StringIO import StringIO

from typing import Dict, List

import common.log_params
from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import params
from common.alignment_storage import AlignmentPiece
from common.line_align import Scorer
from common.save_load import TokenReader, TokenWriter
from common.sequences import Contig
from disjointig_resolve.dataset_simulation import TestDataset
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.correction import Correction
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import DotPlot, LineDotPlot
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.line_extender import LineExtender
from disjointig_resolve.smart_storage import SegmentStorage, AlignmentStorage
from disjointig_resolve.unique_marker import UniqueMarker


class Tester:

    def __init__(self, aligner):
        # type: (Aligner) -> None
        # params.scores = ComplexScores()
        # params.scores.load(open("flye/config/bin_cfg/pacbio_substitutions.mat", "r"))

        self.aligner = aligner
        self.polisher = Polisher(aligner, aligner.dir_distributor)
        testList = []
        for name, obj in inspect.getmembers(sys.modules[__name__]):
            if inspect.isclass(obj) and name.endswith("Test"):
                testList.append(obj)
        self.tests = dict([(c.__name__, c) for c in testList])
        params.redo_alignments = True
        params.k = 500
        params.l = 1500
        params.min_k_mer_cov = 5
        sys.stdout.level = common.log_params.LogPriority.alignment_files - 1
        # sys.stdout.level = params.LogPriority.main_parts

    def testAll(self, fname):
        params = self.readParams(fname)
        fail = 0
        for name, cases in params.items():
            result = self.tests[name]().testAll(cases, self.aligner)
            if not result:
                fail += 1
        print len(params) - fail, "tests passed"
        print fail, "tests failed"
        # self.testLineExtension()
        # self.testLineMerging()
        # self.testLineKnotting()
        # self.testExtensionConsensus()
        # self.testSaves()

    def readParams(self, fname):
        # type: (str) -> Dict[str, List[List[str]]]
        params = dict()  # type: Dict[str, List[List[str]]]
        handler = TokenReader(open(fname, "r"))
        for i in range(handler.readInt()):
            name = handler.readToken()
            assert name is not None
            params[name] = []
            for j in range(handler.readInt()):
                params[name].append(list(handler.readTokens()))
        return params

class SimpleTest:
    def __init__(self):
        self.aligner = None # type: Aligner
        self.scorer = Scorer()

    def testAll(self, params, aligner):
        # type: (List[List[str]], Aligner) -> bool
        self.aligner = aligner
        print "Starting test", self.__class__.__name__
        fails = []
        try:
            self.testManual()
        except AssertionError as e:
            _, _, tb = sys.exc_info()
            fails.append(("Manual", tb, e.message))
        for tn, instance in enumerate(params):
            try:
                self.testCase(instance)
            except AssertionError as e:
                _, _, tb = sys.exc_info()
                fails.append([str(tn), tb, e.message])
        if len(fails) == 0:
            print "Finished test", self.__class__.__name__ + ": Passed"
            return True
        else:
            print "Finished test", self.__class__.__name__ + ": Failed"
            for tn, tb, message in fails:
                print "Failed test " + tn + ":"
                traceback.print_tb(tb)
                print "Message:", message
            return False

    # def runTestSilently(self, testf):
    #     # type: (Callable[[], None]) -> None
    #     tmp1 = sys.stdout # type: OStreamWrapper
    #     tmp2 = sys.stderr # type: OStreamWrapper
    #     tmp1.block()
    #     tmp2.block()
    #     testf()
    #     tmp1.release()
    #     tmp2.release()


    def assertResult(self, res, ethalon):
        # type: (str, str) -> None
        assert res.replace(" ", "") == ethalon, res.replace(" ", "") + " " + ethalon
        pass

    def testCase(self, instance):
        # type: (list[str]) -> None
        pass

    def testManual(self):
        pass


class SegmentStorageTest(SimpleTest):
    def testManual(self):
        contig = Contig("ACGT", "test")
        storage = SegmentStorage()
        storage.add(contig.segment(0,1))
        storage.add(contig.segment(1,2))
        storage.add(contig.segment(2,3))
        storage.add(contig.segment(3,4))
        assert str(storage) == "ReadStorage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        assert str(storage.rc) == "ReadStorage-:[-test[0:1], -test[1:2], -test[2:4-1], -test[3:4-0]]", str(storage.rc)
        storage.mergeSegments(1)
        assert str(storage) == "ReadStorage+:[test[0:1], test[1:2], test[2:4-1], test[3:4-0]]", str(storage)
        storage.mergeSegments()
        assert str(storage) == "ReadStorage+:[test[0:4-0]]", str(storage)
        assert str(storage.rc) == "ReadStorage-:[-test[0:4-0]]", str(storage.rc)
        contig = Contig("ACGTACGTACGTACGT", "test")
        storage = SegmentStorage()
        storage.add(contig.segment(0,5))
        storage.add(contig.segment(10,15))
        assert storage.find(contig.segment(5, 10)) == contig.segment(0,5), str(storage.find(contig.segment(5, 10)))
        assert storage.find(contig.segment(6, 10)) == contig.segment(10,15), str(storage.find(contig.segment(6, 10)))
        assert storage.find(contig.segment(5, 9)) == contig.segment(0,5), str(storage.find(contig.segment(5, 9)))
        assert storage.find(contig.segment(0, 16)) == contig.segment(0,5), str(storage.find(contig.segment(0, 16)))


class AlignmentPieceTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTA", "from")
        contig2 = Contig("ACTACGTACGTACAT", "to")
        al1 = AlignmentPiece(contig1.asSegment(), contig2.segment(0, 8), "2M1I6M")
        al2 = AlignmentPiece(contig1.segment(0, 8), contig2.segment(7, 15), "8M")
        glued = AlignmentPiece.GlueOverlappingAlignments([al1, al2])
        assert glued.cigar == "2M1I5M8M", str(glued) + " " + glued.cigar
        assert  glued.seg_from.Seq() == "ACGTACGTACGTACGT", str(glued) + " " + glued.cigar
        assert al1.reduce(query=contig1.segment(0, 2)).cigar == "2M"
        assert al1.reduce(query=contig1.segment(0, 3)).cigar == "2M"
        assert al1.reduce(query=contig1.segment(0, 4)).cigar == "2M1I1M"


class AlignmentPolishingTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTTAAACGT", "from")
        contig2 = Contig("ACGTTTAACGT", "to")
        al = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al1 = self.scorer.polyshAlignment(al, params.alignment_correction_radius)
        assert al1.cigar == "4M1D2M1I4M", str(al1.asMatchingStrings())
        contig1 = Contig("ACATGATCACT", "from")
        contig2 = Contig("ACGTGAAACGT", "to")
        al = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al1 = self.scorer.polyshAlignment(al, params.alignment_correction_radius)
        assert al1.cigar == "6M1I3M1D1M", str(al1.asMatchingStrings())


class AlignmentStorageTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTACGT", "from")
        contig2 = Contig("ACGTACGTACGT", "to")
        al1 = AlignmentPiece.Identical(contig1.segment(0, 4), contig2.segment(0, 4))
        al2 = AlignmentPiece.Identical(contig1.segment(0, 4), contig2.segment(4, 8))
        al3 = AlignmentPiece.Identical(contig1.segment(4, 8), contig2.segment(8, 12))
        storage = AlignmentStorage()
        storage.addAll([al1, al2, al3])
        assert str(list(storage)) == "[(from[0:4]->to[0:4]:1.000), (from[0:4]->to[4:12-4]:1.000), (from[4:12-4]->to[8:12-0]:1.000)]"
        assert str(list(storage.rc)) == "[(-from[4:12-4]->-to[0:4]:1.000), (-from[8:12-0]->-to[4:12-4]:1.000), (-from[8:12-0]->-to[8:12-0]:1.000)]"
        assert str(list(storage.calculateCoverage())) == "[(to[0:12-0], 1)]"
        assert str(list(storage.filterByCoverage(0, 1))) == "[]"
        assert str(list(storage.filterByCoverage(1, 2))) == "[to[0:12-0]]"
        assert str(list(storage.filterByCoverage(2))) == "[]"
        storage.addAndMergeRight(al3)
        assert str(list(storage)) == "[(from[0:4]->to[0:4]:1.000), (from[0:4]->to[4:12-4]:1.000), (from[4:12-4]->to[8:12-0]:1.000)]"
        al4 = AlignmentPiece.Identical(contig1.segment(2, 8), contig2.segment(2, 8))
        al5 = AlignmentPiece.Identical(contig1.segment(4, 10), contig2.segment(4, 10))
        storage.addAll([al4, al5])
        assert str(list(storage.calculateCoverage())) == "[(to[0:2], 1), (to[2:4], 2), (to[4:12-4], 3), (to[8:12-2], 2), (to[10:12-0], 1)]"
        assert str(list(storage.filterByCoverage(2,3))) == "[to[2:4], to[8:12-2]]"
        assert str(list(storage.filterByCoverage(2))) == "[to[2:12-2]]"
        assert str(list(storage.getAlignmentsTo(contig2.segment(2, 3)))) == "[(from[0:4]->to[0:4]:1.000), (from[2:12-4]->to[2:12-4]:1.000)]"
        assert str(list(storage.getAlignmentsTo(contig2.segment(2, 6)))) == "[(from[2:12-4]->to[2:12-4]:1.000)]"


class AlignmentCompositionTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTACGTACGT", "c1")
        contig2 = Contig("ACGTAGGTACGT", "c2")
        contig3 = Contig("ACTTACGTACGT", "c3")
        al1 = AlignmentPiece.Identical(contig1.asSegment(), contig2.asSegment())
        al2 = AlignmentPiece.Identical(contig2.asSegment(), contig3.asSegment())
        al3 = al1.compose(al2)
        assert al3.__repr__() == "(c1[0:12-0]->c3[0:12-0]:0.92)"
        assert al3.cigar == "12M"
        al4 = al1.reverse()
        al5 = al4.composeTargetDifference(al2)
        assert al5.__repr__() == "(c1[0:12-0]->c3[0:12-0]:0.92)"
        assert al5.cigar == "12M"


class CorrectionMappingTest(SimpleTest):
    def testManual(self):
        contig1 = Contig("ACGTAAAAGGGTACGT", "c1")
        contig2 = Contig("ACGTAAGGGGGTACGT", "c2")
        al = self.scorer.polyshAlignment(AlignmentPiece.Identical(contig1.segment(5, 12), contig2.segment(5, 12)), params.alignment_correction_radius)
        corr = Correction(contig1, contig2, [al])
        assert corr.mapPositionsUp(range(len(contig2))) == [0, 1, 2, 3, 4, 5, 8, 9, 9, 9, 10, 11, 12, 13, 14, 15]
        assert corr.mapPositionsDown(range(len(contig1))) == [0, 1, 2, 3, 4, 5, 6, 6, 6, 9, 10, 11, 12, 13, 14, 15]
        al2 = AlignmentPiece.Identical(contig2.segment(0, 4))
        al3 = AlignmentPiece.Identical(contig2.segment(6, 8))
        al4 = AlignmentPiece.Identical(contig2.segment(6, 16))
        al5 = AlignmentPiece.Identical(contig2.segment(7, 16))
        assert str(corr.composeQueryDifferences([al2, al3, al4, al5])) == "[(c2[0:4]->c1[0:4]:1.000), (c2[6:7]->c1[8:9]:1.000), (c2[6:16-0]->c1[8:16-0]:0.80), (c2[9:16-0]->c1[9:16-0]:1.000)]"


class DotPlotConstructionTest(SimpleTest):
    def testCase(self, instance):
        data = TokenReader(StringIO(" ".join(instance)))
        dataset = TestDataset.loadStructure(data)
        disjointigs = DisjointigCollection()
        for dis in dataset.disjointigs:
            disjointigs.addNew(dis.seq, dis.id)
        dp = DotPlot(disjointigs)
        dp.construct(self.aligner)
        save = StringIO()
        save_handler = TokenWriter(save)
        dp.save(save_handler)
        tmp = save.getvalue()
        test_result = tmp.replace(" ", "").replace("\n", "")
        ethalon = data.readToken()
        if test_result != ethalon:
            for dis in disjointigs:
                print list(dp.allInter(dis.asSegment()))
        assert test_result == ethalon, "\n" + test_result + "\n" + ethalon

class DotPlotModificationTest(SimpleTest):
    def testManual(self):
        self.test1()
        self.test2()
        self.test3()
        self.test4()
        self.test5()

    def test1(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line1 = lines.addNew("ACGTAAAAGGGTACGT", "c1")
        line2 = lines.addNew("ACGTAAGGGGGTACGT", "c2")
        al = self.scorer.polyshAlignment(AlignmentPiece.Identical(line1.asSegment(), line2.asSegment()), params.alignment_correction_radius)
        dp = LineDotPlot(lines, self.aligner)
        dp.addAlignment(al)
        alignment = AlignmentPiece.Identical(Contig("AGG", "tmp").asSegment(), line2.segment(0, 3))
        line2.correctSequence([alignment])
        assert str(list(dp.alignmentsToFrom[line2.id][line1.id])) == "[(c1[0:16-0]->c2[0:16-0]:0.81)]"

    def test2(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line = lines.addNew("ACGTACGTACGT", "c")
        dp = LineDotPlot(lines, self.aligner)
        al1 = AlignmentPiece.Identical(line.segment(0, 8), line.segment(4, 12))
        al2 = AlignmentPiece.Identical(line.segment(0, 4), line.segment(8, 12))
        dp.addAlignment(al1)
        dp.addAlignment(al2)
        alignment = AlignmentPiece.Identical(Contig("AGG", "tmp").asSegment(), line.segment(4, 7))
        line.correctSequence([alignment])
        assert str(list(dp.auto_alignments[
                            "c"])) == "[(c[0:12-4]->c[4:12-0]:0.75), (c[0:4]->c[8:12-0]:1.000), (c[4:12-0]->c[0:12-4]:0.75), (c[8:12-0]->c[0:4]:1.000), (c[0:12-0]->c[0:12-0]:1.000)]"

    def test3(self):
        lines = NewLineStorage(DisjointigCollection(), self.aligner)
        line = lines.addNew("ACGTACGTACGT", "c")
        dp = LineDotPlot(lines, self.aligner)
        al1 = AlignmentPiece.Identical(line.segment(0, 8), line.segment(4, 12))
        al2 = AlignmentPiece.Identical(line.segment(0, 4), line.segment(8, 12))
        dp.addAlignment(al1)
        dp.addAlignment(al2)
        alignment = AlignmentPiece.Identical(Contig("TCC", "tmp").asSegment(), line.segment(3, 6))
        line.correctSequence([alignment])
        assert str(list(dp.auto_alignments["c"])) == "[(c[1:12-4]->c[5:12-0]:0.86), (c[0:4]->c[8:12-0]:1.000), (c[5:12-0]->c[1:12-4]:0.86), (c[8:12-0]->c[0:4]:1.000), (c[0:12-0]->c[0:12-0]:1.000)]"

    def test4(self):
        dataset = TestDataset("abcABC")
        name = dataset.addContig("abcAB")
        lines, dp, reads = dataset.genAll(self.aligner)
        line = lines[name]
        line.extendRight(dataset.alphabet["C"].seq)
        assert str(list(dp.auto_alignments[line.id])) == "[(C0_abcAB[0:1650]->C0_abcAB[1650:3302-0]:0.995), (C0_abcAB[1650:3302-0]->C0_abcAB[0:1650]:0.995), (C0_abcAB[0:3302-0]->C0_abcAB[0:3302-0]:1.000)]", str(list(dp.auto_alignments[line.id]))

    def test5(self):
        dataset = TestDataset("abcABC")
        name1 = dataset.addContig("abc")
        name2 = dataset.addContig("ABC")
        lines, dp, reads = dataset.genAll(self.aligner)
        line = lines[name1]
        sa = dataset.alphabet["a"].seq
        sb = dataset.alphabet["b"].seq
        tmp = Contig(sa + "ACGACAGTAACTTGAACGACAGTAACTTGAACGACAGTAACTTGAACGACAGTAACTTGAACGACAGTAACTTGAACGACAGTAACTTGAACGACAGTAACTTGA" + sb, "tmp")
        al1 = AlignmentPiece.Identical(tmp.prefix(len=len(sa)), line.prefix(len=len(sa)))
        al2 = AlignmentPiece.Identical(tmp.asSegment().suffix(length=len(sb)), line.segment(len(sa), len(sa) + len(sb)))
        al = AlignmentPiece.MergeFittingAlignments([al1, al2])
        line.correctSequence([al])
        assert str(list(dp.allInter(line.asSegment()))) == "[(C0_abc[0:1755-0]->C0_abc[0:1755-0]:1.000), (C1_ABC[0:1652-0]->C0_abc[0:1755-0]:0.94)]"


class PotentialReadsTest(SimpleTest):
    def testManual(self):
        dataset = TestDataset("abcdefgh")
        dataset.addDisjointig("abcdefghabcdefgh")
        name = dataset.addContig("abcdefgh")
        dataset.generateReads(4, 2, True)
        lines, dp, reads = dataset.genAll(self.aligner)
        line = lines[name]
        assert str(list(line.getRelevantAlignmentsFor(line.asSegment()))) == "[(R0_abcd[0:2189-0]->C0_abcdefgh[0:2200]:0.96), (R1_bcde[0:2189-0]->C0_abcdefgh[550:2750]:0.97), (R2_cdef[3:2198-0]->C0_abcdefgh[1102:3300]:0.97), (R3_defg[0:2198-0]->C0_abcdefgh[1650:3850]:0.97), (R4_efgh[1:2205-0]->C0_abcdefgh[2200:4400-0]:0.97), (R5_fgha[0:1642]->C0_abcdefgh[2750:4400-0]:0.97), (R5_fgha[1642:2186-2]->C0_abcdefgh[0:547]:0.97), (R6_ghab[0:1112]->C0_abcdefgh[3300:4400-0]:0.96), (R6_ghab[1112:2210-0]->C0_abcdefgh[0:1100]:0.97), (R7_habc[543:2201-0]->C0_abcdefgh[0:1650]:0.97), (R7_habc[0:543]->C0_abcdefgh[3850:4400-0]:0.96)]", str(list(line.getRelevantAlignmentsFor(line.asSegment())))


class UniqueRegionMarkingTest(SimpleTest):
    def testCase(self, instance):
        data = TokenReader(StringIO(" ".join(instance)))
        dataset = TestDataset.loadStructure(data)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        ethalon1 = data.readToken()
        ethalon2 = data.readToken()
        line = lines[dataset.contigs[0].id]
        self.assertResult(str(line.correct_segments), ethalon1)
        self.assertResult(str(line.completely_resolved), ethalon2)


class ReadRecruitmentTest(SimpleTest):
    def testManual(self):
        dataset = TestDataset("abcdefghijklmCDEFGHInopqr")
        dname = dataset.addDisjointig("abcdefghijklmCDEFGHInopqr".upper())
        name1 = dataset.addContig("abcde")
        name2 = dataset.addContig("klmCDE")
        dataset.generateReads(4, 5, True)
        # dataset.saveStructure(TokenWriter(sys.stdout))
        lines, dp, reads = dataset.genAll(self.aligner)
        # UniqueMarker().markAllUnique(lines, dp)
        line1 = lines[name1]
        line1.correct_segments.add(line1.asSegment())
        line1.completely_resolved.add(line1.asSegment())
        line2 = lines[name2]
        line2.correct_segments.add(line2.asSegment())
        line2.completely_resolved.add(line2.asSegment())
        extender = LineExtender(self.aligner, None, lines.disjointigs, dp)
        res = extender.attemptCleanResolution(line1.asSegment())
        assert str(res[0][1]) == "[(R2_bcde[0:2200-0]->C0_abcde[550:2750-0]:0.97), (R3_bcde[4:2192-0]->C0_abcde[553:2750-0]:0.97), (R4_cdef[0:1657]->C0_abcde[1100:2750-0]:0.96), (R5_cdef[0:1656]->C0_abcde[1100:2750-0]:0.96)]", str(res[0][1])
        assert str(res[1][1]) == "[(R24_mCDE[0:2201-0]->C1_klmCDE[1100:3298-0]:0.97), (R25_mCDE[0:2194-0]->C1_klmCDE[1100:3298-0]:0.97), (R27_CDEF[0:1658]->C1_klmCDE[1651:3298-0]:0.96)]", str(res[1][1])


class KnottingTest(SimpleTest):
    def testManual(self):
        dataset = TestDataset("abcdefgabhiDEFjkl")
        dname = dataset.addDisjointig("abcdefgabhiCDEjklabcdefgabhiCDEjkl".upper())
        dataset.generateReads(5, 15, True)
        read1 = dataset.addRead("cdefg")
        read2 = dataset.addRead("cdefg")
        name1 = dataset.addContig("abcde")
        name2 = dataset.addContig("efgabhi")
        lines, dp, reads = dataset.genAll(self.aligner)
        read1 = reads[read1]
        read2 = reads[read2]
        line1 = lines[name1]
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        knotter = LineMerger(lines, Polisher(self.aligner, self.aligner.dir_distributor), dp)
        dp.printAll(sys.stdout)
        res = knotter.tryMergeRight(line1)
        assert res is not None
        assert str(list(dp.allInter(res.asSegment()))) == \
               "[((C0_abcde,C1_efgabhi)[0:1100]->(C0_abcde,C1_efgabhi)[3850:4950]:1.000!!!), ((C0_abcde,C1_efgabhi)[3850:4950]->(C0_abcde,C1_efgabhi)[0:1100]:1.000!!!), ((C0_abcde,C1_efgabhi)[0:6050-0]->(C0_abcde,C1_efgabhi)[0:6050-0]:1.000)]", str(list(dp.allInter(res.asSegment())))
               # "[((C0_abcde,C1_efgabhi)[0:1100]->(C0_abcde,C1_efgabhi)[3850:4948]:0.997!!!), ((C0_abcde,C1_efgabhi)[3850:4948]->(C0_abcde,C1_efgabhi)[0:1100]:0.997!!!), ((C0_abcde,C1_efgabhi)[0:6048-0]->(C0_abcde,C1_efgabhi)[0:6048-0]:1.000)]

class StructureUpdatingTest(SimpleTest):
    def testManual(self):
        self.test1()
        # Incorrect calculation of unique segments. Should not occur on real datasets
        # self.test2()

    def test1(self):
        dataset = TestDataset("abcdefghijklmCDEFGHInopqr")
        dname = dataset.addDisjointig("abcdefghijklmCDEFGHInopqrabcdefghijklmCDEFGHInopqr".upper())
        name1 = dataset.addContig("abcde")
        name2 = dataset.addContig("klmCDE")
        dataset.generateReads(4, 20, True)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        line1 = lines[name1]
        line2 = lines[name2]
        extender = LineExtender(self.aligner, None, lines.disjointigs, dp)
        extender.updateAllStructures(list(line1.correct_segments))
        print str(line1.correct_segments), str(line1.completely_resolved), str(line2.correct_segments), str(line2.completely_resolved)
        assert str(line1.correct_segments) == "ReadStorage+:[C0_abcde[0:2200]]", str(line1.correct_segments)
        assert str(line1.completely_resolved) == "ReadStorage+:[C0_abcde[0:2098]]", str(line1.completely_resolved)
        assert str(line2.correct_segments) == "ReadStorage+:[C1_klmCDE[0:2749]]", str(line2.correct_segments)
        assert str(line2.completely_resolved) == "ReadStorage+:[C1_klmCDE[0:2649]]", str(line2.completely_resolved)

    def test2(self):
        dataset = TestDataset("abcdefgcijklmCDEFGHInopqr")
        dname = dataset.addDisjointig("abcdefgcijklmCDEFGHInopqrabcdefgcijklmCDEFGHInopqr".upper())
        name1 = dataset.addContig("abcde")
        name2 = dataset.addContig("klmCDE")
        dataset.generateReads(4, 20, True)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        line1 = lines[name1]
        line2 = lines[name2]
        extender = LineExtender(self.aligner, None, lines.disjointigs, dp)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
        # extender.updateAllStructures(list(line1.correct_segments))
        print line1, line2
        print str(line1.correct_segments)
        print str(line1.completely_resolved)
        print str(line2.correct_segments)
        print str(line2.completely_resolved)
        assert str(line1.correct_segments) == "ReadStorage+:[C0_abcde[550:3850]]", str(line1.correct_segments)
        assert str(line1.completely_resolved) == "ReadStorage+:[C0_abcde[550:3000], C0_abcde[3300:3845]]", str(line1.completely_resolved)
        assert str(line2.correct_segments) == "ReadStorage+:[C1_klmCDE[550:4395]]", str(line2.correct_segments)
        assert str(line2.completely_resolved) == "ReadStorage+:[C1_klmCDE[851:3549], C1_klmCDE[3851:4395]]", str(line2.completely_resolved)


class LineExtensionTest(SimpleTest):
    def testCase(self, instance):
        # type: (list[str]) -> None
        dataset = TestDataset(instance[0], mutation_rate=0.01)
        dname = dataset.addDisjointig(instance[0] + instance[0].upper())
        dataset.generateReads(int(instance[1]), 25, True)
        ethalon = int(instance[2])
        for s in instance[3:]:
            dataset.addContig(s)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        knotter = LineMerger(lines, Polisher(self.aligner, self.aligner.dir_distributor), dp)
        extender = LineExtender(self.aligner, knotter, lines.disjointigs, dp)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
        while True:
            stop = True
            for line_id in list(lines.items.keys()):
                if line_id not in lines.items:
                    continue
                line = lines[line_id]
                dp.printAll(sys.stdout)
                extended = extender.processLine(line)
                if extended:
                    stop = False
            if stop:
                break
        print " ".join([str(dataset.translateBack(line, self.aligner)) for line in lines.unique()])
        print [line.circular for line in lines.unique()]
        breaks = 0
        for line in lines.unique():
            if not line.circular:
                breaks += 1
        assert breaks == ethalon, str(breaks) + " " + str(ethalon)

    def testManual(self):
        pass
        # self.test1()

    def test1(self):
        dataset = TestDataset("abcdefghijklmCDEFGHInopqr", mutation_rate=0.01)
        dname = dataset.addDisjointig("abcdefghijklmCDEFGHInopqrabcd".upper())
        name1 = dataset.addContig("abcde")
        name2 = dataset.addContig("klmCDE")
        dataset.generateReads(5, 25, True)
        lines, dp, reads = dataset.genAll(self.aligner)
        UniqueMarker(self.aligner).markAllUnique(lines, reads)
        line1 = lines[name1]
        line2 = lines[name2]
        knotter = LineMerger(lines, Polisher(self.aligner, self.aligner.dir_distributor), dp)
        extender = LineExtender(self.aligner, knotter, lines.disjointigs, dp)
        print "New iteration results"
        print dataset.translateBack(line1, self.aligner), dataset.translateBack(line2, self.aligner)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
        while True:
            stop = True
            for line_id in list(lines.items.keys()):
                if line_id not in lines.items:
                    continue
                line = lines[line_id]
                dp.printAll(sys.stdout)
                extended = extender.processLine(line)
                if extended:
                    stop = False
            if stop:
                break
        print " ".join([str(dataset.translateBack(line, self.aligner)) for line in lines.unique()])
        print [line.circular for line in lines.unique()]

class WindowMaxTest(SimpleTest):
    def testManual(self):
        random.seed(0)
        arr = []
        tmp_arr = []
        for i in range(1000):
            arr.append((random.randint(0, 100), i))
            tmp_arr.append(arr[-1][0])
        res = Scorer().maxInRange(arr, 50)
        for i, val in enumerate(res):
            assert val == max(tmp_arr[max(i - 50, 0): min(i + 50 + 1, len(arr))])