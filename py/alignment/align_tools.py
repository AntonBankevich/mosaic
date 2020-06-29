import itertools
import os
import shutil
import subprocess
import sys

sys.path.append("py")
import common.log_params
from common.save_load import TokenWriter, TokenReader
from common.seq_records import NamedSequence
from flye_tools.alignment import make_alignment
from common.sequences import ContigCollection, Segment, Contig, \
    ContigStorage
from common.alignment_storage import AlignmentPiece, ReadCollection
from typing import Iterable, Tuple, Generator, BinaryIO, List
from common import basic, sam_parser, SeqIO, params

class AlignmentException(Exception):
    pass

class DirDistributor:
    def __init__(self, dir):
        basic.ensure_dir_existance(dir)
        self.dir = dir
        self.cur_dir = 0

    def nextDir(self):
        # type: () -> str
        name = os.path.join(self.dir, str(self.cur_dir))
        if params.save_alignments:
            self.cur_dir += 1
        assert self.cur_dir <= 100000
        basic.ensure_dir_existance(name)
        return name

    def CheckSequences(self, reads, reads_file):
        # type: (Iterable[NamedSequence], str) -> bool
        if not os.path.exists(reads_file):
            return False
        try:
            for rec, read in itertools.izip_longest(SeqIO.parse_fasta(open(reads_file, "r")), reads):
                if str(rec.id) != str(read.id) or rec.seq != read.seq:
                    return False
            return True
        except:
            return False

    def WriteSequences(self, reads, reads_file):
        # type: (Iterable[NamedSequence], str) -> None
        f = open(reads_file, "w")
        for read in reads:
            SeqIO.write(read, f, "fasta")
        f.close()

    def hash(self, reads):
        # type: (Iterable[NamedSequence]) -> int
        res = 0
        for read in reads:
            res += read.seq.__hash__() + read.id.__hash__()
        return res

    def calculateHash(self, content):
        # type: (List[Tuple[Iterable[NamedSequence], str]]) -> Generator[Tuple[str, str, str]]
        for reads, f_name in content:
            yield f_name, str(self.hash(reads)), str(len(reads))

    def printHash(self, handler, hashs):
        # type: (BinaryIO, List[Tuple[str, str, str]]) -> None
        for rec in hashs:
            handler.write(" ".join(rec) + "\n")

    def compareHash(self, handler, hashs):
        # type: (BinaryIO, List[Tuple[str, str, str]]) -> bool
        lines = handler.readlines()
        if len(lines) != len(hashs):
            return False
        for l, rec in zip(lines, hashs):
            l = l.split()
            if len(l) != len(rec):
                return False
            for s1, s2 in zip(l, rec):
                if s1 != s2:
                    return False
        return True

    def fillNextDir(self, content):
        # type: (List[Tuple[Iterable[NamedSequence], str]]) -> Tuple[str, list[str], bool]
        same = True
        dir = self.nextDir()
        content_files = []
        for reads, f_name in content:
            content_files.append(os.path.join(dir, f_name))
        hash_file = os.path.join(dir, "hashs.txt")
        hashs = list(self.calculateHash(content))
        if os.path.isfile(hash_file) and self.compareHash(open(hash_file, "r"), hashs) and not params.clean:
            return dir, content_files, True
        basic.recreate(dir)
        for reads, f_name in content:
            f_name = os.path.join(dir, f_name)
            self.WriteSequences(reads, f_name)
        self.printHash(open(hash_file, "w"), hashs)
        return dir, content_files, False

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeToken(self.dir)
        handler.writeIntLine(self.cur_dir)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> DirDistributor
        res = DirDistributor(handler.readToken())
        res.cur_dir = handler.readInt()
        return res


class Aligner:
    def __init__(self, dir_distributor, threads = params.threads):
        # type: (DirDistributor, int) -> Aligner
        self.dir_distributor = dir_distributor
        self.threads = threads
        self.filters = dict()
        self.filters["overlap"] = lambda als: self.filterOverlaps(als)
        self.filters["dotplot"] = lambda als: self.filterLocal(als)
        self.filters["local"] = lambda als: als
        self.filters["reference"] = lambda als: als

    def alignReadCollection(self, reads_collection, contigs):
        # type: (ReadCollection, Iterable[Contig]) -> None
        contig_collection = ContigCollection(contigs)
        contig_ids = set()
        for contig in contigs:
            if contig.rc.id not in contig_ids:
                contig_ids.add(contig.id)
        read_ids = set()
        for read in reads_collection:
            if read.rc.id not in read_ids:
                read_ids.add(read.id)
        contigs = filter(lambda contig: contig.id in contig_ids, contigs)
        reads = filter(lambda read: read.id in read_ids, reads_collection)
        reads_collection.fillFromSam(self.align(reads, contigs, "local"), contig_collection)

    def alignReadsToSegments(self, reads, segments):
        # type: (ReadCollection, Iterable[Segment]) -> None
        segments = list(segments)
        seg_dict = dict()
        for i, seg in enumerate(segments):
            seg_dict[str(i + 1)] = seg
        contigs = map(lambda (i, seg): Contig(seg.Seq(), str(i + 1)), enumerate(segments))
        read_collection = ReadCollection().extendClean(reads)
        self.alignReadCollection(read_collection, ContigCollection(contigs))
        read_collection.contigsAsSegments(seg_dict)
        reads.mergeAlignments(read_collection)

    def realignCollection(self, reads_collection, contigs):
        # type: (ReadCollection, Iterable[Contig]) -> None
        for read in reads_collection:
            read.clean()
        self.alignReadCollection(reads_collection, contigs)

    # def fixExtendedLine(self, line):
    #     # type: (Line) -> None
    #     toFix = []
    #     for read in line.reads:
    #         if not read.noncontradicting(line.asSegment()):
    #             toFix.append(read)
    #     newAlignments = ReadCollection(line.reads.contigs).extend(toFix)
    #     self.alignReadCollection(newAlignments)
    #     for read in newAlignments:
    #         for al in line.reads[read.id].alignments:
    #             if not al.contradicting(line.asSegment()):
    #                 continue
    #             for new_al in read.alignments:
    #                 if new_al.contains(al) and len(new_al) > len(al):
    #                     al.seg_from = new_al.seg_from
    #                     al.seg_to = new_al.seg_to
    #                     al.cigar = new_al.cigar

    # def separateAlignments(self, reads, contigs):
    #     # type: (Iterable[NamedSequence], Iterable[Contig]) -> ReadCollection
    #     contigs_collection = ContigCollection(list(contigs))
    #     res = ReadCollection(contigs_collection)
    #     for read in reads:
    #         res.addNewRead(NamedSequence(read.seq, read.id)) # remove when all ids are str
    #     for contig in contigs:
    #         res.fillFromSam(self.align(res, [contig]), contigs_collection)
    #     return res

    def align_files(self, reference_file, reads_files, threads, platform, mode, out_alignment):
        out_file = out_alignment + "_initial"
        cmdline = [params.MINIMAP_BIN, reference_file]
        cmdline.extend(reads_files)
        cmdline.extend(["-N1000", "-a", "-Q", "-w5", "-m100", "-g10000", "--max-chain-skip",
                        "25", "-t", str(threads)])
        # cmdline.extend(["-x", "ava-pb"])
        if mode == "dotplot":
            cmdline.append("-p0.00")
        elif mode in ["overlap", "local"]:
            cmdline.append("-p0.1")
            cmdline.append("-DP")
        if mode == "reference":
            cmdline.extend(["-x", "asm5"])
        else:
            if platform == "nano":
                cmdline.append("-k15")
            else:
                cmdline.append("-Hk19")
        try:
            devnull = open(os.devnull, "w")
            sys.stdout.log(common.log_params.LogPriority.alignment_files, "Performing alignment:", " ".join(cmdline))
            subprocess.check_call(cmdline, stderr=devnull, stdout=open(out_file, "w"))
            env = os.environ.copy()
            env["LC_ALL"] = "C"
            subprocess.check_call(["sort", "-T", os.path.dirname(out_alignment), out_file],
                                  stdout=open(out_alignment, "w"), env=env)
            # os.remove(out_file)
        except (subprocess.CalledProcessError, OSError) as e:
            raise AlignmentException(str(e))

    def align(self, reads, reference, mode):
        # type: (Iterable[NamedSequence], Iterable[Contig], str) -> sam_parser.Samfile
        reference = list(reference)
        dir, new_files, same = self.dir_distributor.fillNextDir([(reference, "contigs.fasta"), (list(reads), "reads.fasta")])
        contigs_file = new_files[0]
        reads_file = new_files[1]
        alignment_dir = os.path.join(dir, "alignment")
        alignment_file = os.path.join(dir, "alignment.sam")
        basic.ensure_dir_existance(dir)
        basic.ensure_dir_existance(alignment_dir)
        if same and not params.clean and os.path.exists(alignment_file):
            sys.stdout.log(common.log_params.LogPriority.alignment_files, "Alignment reused:", alignment_file)
            pass
        else:
            if os.path.isfile(alignment_file):
                os.remove(alignment_file)
            self.align_files(contigs_file, [reads_file], self.threads, params.technology, mode, alignment_file)
        return sam_parser.Samfile(open(alignment_file, "r"))

    def alignAndFilter(self, reads, ref_storage, mode):
        # type: (Iterable[Contig], ContigStorage, str) -> Generator[AlignmentPiece]
        filter = self.filters[mode]
        read_storage = ContigStorage(reads, False)
        cnt = 0
        cnt_out = 0
        als = []
        for rec in self.align(read_storage, list(ref_storage.unique()), mode):
            if rec.is_unmapped:
                continue
            cnt += 1
            if len(als) > 0 and rec.query_name != als[0].seg_from.contig.id:
                cnt += len(als)
                res = list(filter(als))
                cnt_out += len(res)
                for al in res:
                    yield al
                als = []
            if len(als) > 0:
                seq_from = als[0].seg_from.contig
            else:
                seq_from = read_storage[rec.query_name]
            seq_to = ref_storage[rec.tname]
            tmp = AlignmentPiece.FromSamRecord(seq_from, seq_to, rec)
            if tmp is not None:
                if mode == "dotplot":
                    als.extend(tmp.splitRef())
                elif (mode == "local") and tmp.indelLength * 8 > tmp.matchingPositionsCount:
                    als.extend(tmp.splitRead())
                else:
                    als.append(tmp)
        if len(als) > 0:
            cnt += len(als)
            res = list(filter(als))
            cnt_out += len(res)
            for al in res:
                yield al

    def dotplotAlign(self, reads, ref_storage):
        # type: (Iterable[Contig], ContigStorage) -> Generator[AlignmentPiece]
        for al in self.alignAndFilter(reads, ref_storage, "dotplot"):
            yield al

    def overlapAlign(self, reads, ref_storage):
        # type: (Iterable[Contig], ContigStorage) -> Generator[AlignmentPiece]
        for al in self.alignAndFilter(reads, ref_storage, "overlap"):
            yield al

    def localAlign(self, reads, ref_storage):
        # type: (Iterable[Contig], ContigStorage) -> Generator[AlignmentPiece]
        for al in self.alignAndFilter(reads, ref_storage, "local"):
            yield al

    def referenceAlign(self, reads, ref_storage):
        # type: (Iterable[Contig], ContigStorage) -> Generator[AlignmentPiece]
        for al in self.alignAndFilter(reads, ref_storage, "reference"):
            yield al

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeIntLine(self.threads)
        self.dir_distributor.save(handler)

    @staticmethod
    def load(handler):
        # type: (TokenReader) -> Aligner
        threads = handler.readInt()
        return Aligner(DirDistributor.load(handler), threads)

    def filterLocal(self, als):
        # type: (List[AlignmentPiece]) -> List[AlignmentPiece]
        all_res = []
        if len(als) == 1:
            return als
        als = sorted(als, key = lambda al: al.seg_to.contig.id)
        for k, iter in itertools.groupby(als, key=lambda al: (al.seg_from.contig.id, al.seg_to.contig.id)):
            res = [] # type: List[AlignmentPiece]
            iter = sorted(iter, key = lambda al: al.seg_from.left)
            for al in iter:
                ok = True
                for i, al1 in enumerate(res):
                    if al.seg_from.inter(al1.seg_from) and al.seg_to.inter(al1.seg_to):
                        if al.deepInter(al1):
                            if al.matchingPositionsCount > al1.matchingPositionsCount + 50 or \
                                    (al.matchingPositionsCount > al1.matchingPositionsCount - 50 and al.indelLength < al1.indelLength):
                                res[i] = al
                            ok = False
                            break
                if ok:
                    res.append(al)
            all_res.extend(res)
        return all_res

    def filterOverlaps(self, als):
        # type: (List[AlignmentPiece]) -> List[AlignmentPiece]
        als = filter(lambda al: not al.contradictingRTC(tail_size=params.bad_end_length), als)
        return self.filterLocal(als)

if __name__ == "__main__":
    dir = sys.argv[1]
    query = sys.argv[2]
    target = sys.argv[3]
    extra_params = sys.argv[4:]
    contra = "contra" in extra_params
    over = "over" in extra_params
    long = "long" in extra_params
    start = "start" in extra_params
    forward = "forward" in extra_params
    aln = Aligner(DirDistributor(dir))
    basic.CreateLog(dir)
    contigs = ContigCollection().loadFromFasta(open(target, "r"), False)
    for al in aln.localAlign(ReadCollection().loadFromFile(query), contigs):
        if start:
            if al.seg_to.contig.id.startswith("-"):
                al = al.rc
            if al.seg_to.left > 50:
                continue
        if over and al.contradictingRTC():
            continue
        if forward:
            if al.seg_to.contig.id.startswith("-"):
                al = al.rc
        if contra and (len(al) < 8000 or not al.contradictingRTC()):
            continue
        if long and len(al) < 5000:
            continue
        sys.stdout.write(str(len(al)) + " ")
        sys.stdout.write(str(al))
        if len(al) > 40000:
            print ""
            continue
        sys.stdout.write("\n")
        s = list(al.asMatchingStrings())
        print s[0]
        for a, b in zip(s[0], s[1]):
            if a != b:
                sys.stdout.write("-")
            else:
                sys.stdout.write("|")
        print ""
        print s[1]



