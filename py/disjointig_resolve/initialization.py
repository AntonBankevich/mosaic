import itertools
import os
import subprocess
import sys

from typing import Dict, List, Tuple

from alignment.align_tools import Aligner
from alignment.polishing import Polisher
from common import params, basic
from common.SimpleGraph import SimpleGraph
from common.alignment_storage import ReadCollection, AlignmentPiece, AlignedRead
from common.dot_parser import DotParser
from common.save_load import TokenReader
from common.sequences import ContigCollection, ContigStorage, Contig
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.analysis import CoverageAnalyser
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.unique_marker import UniqueMarker

def adjustKL(aligner, reads, contigs):
    # type: (Aligner, ReadCollection, ContigStorage) -> None
    analyser = CoverageAnalyser(aligner, reads)
    sys.stdout.info("Analysing k-mer coverage by reads.")
    analyser.analyseSegmentCoverage(contigs)
    # analyser.printAnalysis()
    newK = analyser.chooseK(params.k_cov) - 2 * params.bad_end_length
    newL = analyser.chooseK(params.l_cov) - 2 * params.bad_end_length
    sys.stdout.info("Chosen k and l:", newK, newL)
    newK = min(newK, newL - 300)
    newL = min(newL, newK + 1000)
    sys.stdout.info("Adjusted k and l:", newK, newL)
    if newK > params.k:
        sys.stdout.info("k and l adjusted")
        params.k = newK
        params.l = newL
    else:
        sys.stdout.info("k and l not adjusted since k is too small. Returning to default values")


def CreateLineCollection(dir, aligner, contigs, disjointigs, reads, split):
    sys.stdout.info("Creating line collection")
    lines = NewLineStorage(disjointigs, aligner)
    if split:
        lines.splitFromContigs(contigs)
    else:
        lines.fillFromContigs(contigs)
    lines.writeToFasta(open(os.path.join(dir, "initial_lines.fasta"), "w"))
    lines.alignDisjointigs()
    # lines.fillFromDisjointigs()
    sys.stdout.info("Marking unique regions")
    UniqueMarker(aligner).markAllUnique(lines, reads)
    for line in lines.unique():
        sys.stdout.trace(line, line.completely_resolved)
    if not split:
        sys.stdout.trace("Splitting lines into parts")
        line_list = list(lines.unique()) # type: List[NewLine]
        while len(line_list) > 0:
            line = line_list.pop()
            sys.stdout.trace("Splitting line", line)
            if len(line.completely_resolved) > 1:
                sys.stdout.trace("Splitted line because of two resolved segments:", str(line.completely_resolved))
                left = line.completely_resolved[1].left
                right = line.completely_resolved[0].right
                if left >= right:
                    right = left + 1
                line1, line2 = lines.splitLine(line.segment(left, right))
                line_list.extend([line1, line2])
            else:
                if line.initial[-1].seg_to.right + 5000 < len(line):
                    sys.stdout.trace("Cut line on the right because too long unresolved segment:", line, str(line.completely_resolved))
                    line.cutRight(line.initial[-1].seg_to.right + 3000)
                if line.initial[0].seg_to.left > 5000:
                    sys.stdout.trace("Cut line of the left because too long unresolved segment:", line, str(line.completely_resolved))
                    line.rc.cutRight(len(line) - line.initial[0].seg_to.left + 3000)
                if len(line.completely_resolved[0]) > 40000:
                    sys.stdout.trace("Splitted line because it is too long", str(line.completely_resolved))
                    line12, line3 = lines.splitLine(line.completely_resolved[0].suffix(length=max(10000, params.k + params.bad_end_length)).prefix(length=1000))
                    line1, line2 = lines.splitLine(line12.completely_resolved[0].prefix(length=max(10000, params.k + params.bad_end_length)).suffix(length=1000))
                    line1.tie(line2, -1000, "")
                    line2.tie(line3, -1000, "")
        line_list = sorted(lines.unique(), key=lambda line: line.id)
        sys.stdout.trace("Final list of lines:")
        for line in line_list:
            sys.stdout.trace(line, line.completely_resolved)
    lines.writeToFasta(open(os.path.join(dir, "initial_prolonged_lines.fasta"), "w"))
    return lines


def LoadLineCollection(dir, lc_file, aligner, contigs, disjointigs, reads, polisher):
    # type: (str, str, Aligner, ContigStorage, DisjointigCollection, ReadCollection, Polisher) -> NewLineStorage
    sys.stdout.info("Initializing lines from init file", lc_file)
    lines = NewLineStorage(disjointigs, aligner)
    f = TokenReader(open(lc_file, "r"))
    n = f.readInt()
    for i in range(n):
        id = f.readToken()
        contig = contigs[id]
        assert contig.id == id
        line = lines.addNew(contig.seq, contig.id)
        read_ids = f.readTokens()
        for al in aligner.overlapAlign([reads[rid] for rid in read_ids], ContigStorage([line])):
            if len(al.seg_to) >= min(params.k, len(line) - 100):
                tmp_line = al.seg_to.contig # type: NewLine
                tmp_line.addReadAlignment(al)
        if len(line) < params.k + 200:
            new_contig, new_als = polisher.polishEnd(list(line.read_alignments), max_extension=params.k + 100 - len(line))
            line.extendRight(new_contig.suffix(pos=len(line)).Seq(), new_als)
        line.correct_segments.add(line.asSegment().shrink(100))
        line.completely_resolved.add(line.asSegment().shrink(100))
        line.initial.add(AlignmentPiece.Identical(line.asSegment().asContig().asSegment(), line.asSegment()))
    sys.stdout.trace("Final list of lines:")
    for line in lines.unique():
        sys.stdout.trace(line, line.completely_resolved)
    lines.writeToFasta(open(os.path.join(dir, "initial_lines.fasta"), "w"))
    lines.alignDisjointigs()
    sys.stdout.info("Constructing line dot plot")
    return lines


def CreateDisjointigCollection(d_files, dir, aligner, reads):
    sys.stdout.info("Creating disjointig collection")
    disjointigs = DisjointigCollection()
    for f in d_files:
        disjointigs.loadFromFasta(open(f, "r"))
    sys.stdout.info("Extending disjointig collection")
    clen = 5000000
    bad_reads = reads.cleanCopy()
    tlen0 = sum(map(len, bad_reads))
    good_reads = set()
    for al in aligner.localAlign(reads, disjointigs):
        if not al.contradictingRTC(al.seg_to.contig.asSegment(), params.bad_end_length) and len(al.seg_from.contig) > len(al) - 2 * params.bad_end_length:
            good_reads.add(al.seg_from.contig.id)
    sys.stdout.info("Fraction of reads without full alignment to disjointigs:", 1 - float(len(good_reads)) / len(reads))
    if len(good_reads) > 0.99 * len(bad_reads):
        sys.stdout.info("Alomst all reads have good alignments. Skipping disjointig collection extension.")
        disjointigs.addAlignments(aligner.localAlign(reads, disjointigs))
        return disjointigs
    rf = os.path.join(dir, "badreads.fasta")
    bad_reads = bad_reads.filter(lambda read: read.id not in good_reads)
    tlen = sum(map(len, bad_reads))
    bad_reads.print_fasta(open(rf, "w"))
    l = tlen * clen / tlen0
    assembly_dir = os.path.join(dir, "assembly0")
    disjointigs_file = constructDisjointigs(bad_reads, l, assembly_dir)
    code = 0

    if code == 0:
        disjointigs.loadFromFasta(open(disjointigs_file, "r"))
        sys.stdout.trace("Disjointigs:")
        for dis in disjointigs:
            sys.stdout.trace(dis.id, len(dis))
        disjointigs.writeToFasta(open(os.path.join(dir, "disjointigs.fasta"), "w"))
    else:
        sys.stdout.trace("Could not assemble new disjointigs")
    sys.stdout.info("Aligning reads to disjointigs")
    disjointigs.addAlignments(aligner.localAlign(reads, disjointigs))
    return disjointigs


def CreateContigCollection(graph_file, contigs_file, min_cov, aligner, polisher, reads, force_unique, all_unique):
    sys.stdout.info("Creating contig collection")
    if force_unique is None and not all_unique:
        graph = SimpleGraph().ReadDot(graph_file)
        graph.FillSeq(contigs_file)
        covs = []
        for e in graph.e.values():
            covs.append((e.len, e.cov))
        tmp_cov = []
        total = sum(l for c,l in covs) / 2
        for l, c in sorted(covs)[::-1]:
            if total < 0:
                break
            tmp_cov.append((l, c))
            total -= l
        avg_cov = float(sum([l * c for l, c in tmp_cov])) / sum(l for l, c in tmp_cov)
        sys.stdout.info("Average coverage determined:", avg_cov)
        nonunique = set()
        for edge in graph.e.values():
            if edge.unique and edge.len < 20000 and len(graph.v[edge.start].out) > 1:
                if edge.cov >= min_cov and (edge.cov < 0.8 * avg_cov or edge.len > 40000):
                    alter = ContigStorage()
                    for e in graph.v[edge.start].out:
                        if e != edge:
                            alter.add(Contig(e.seq, e.id))
                    for al in aligner.localAlign([Contig(edge.seq, edge.id)], alter):#type: AlignmentPiece
                        if al.percentIdentity() > 0.98 and (al.seg_from.left < 100 and al.seg_to.left < 100 and len(al) > min(500, edge.len)):
                            nonunique.add(edge.id)
                            nonunique.add(basic.Reverse(edge.id))
        contigs = ContigCollection()
        for edge in graph.e.values():
            if basic.isCanonocal(edge.id):
                if edge.unique and (edge.len > params.min_isolated_length or len(graph.v[edge.end].out) > 0 or len(graph.v[edge.start].inc) > 0):
                    if edge.cov >= min_cov and (edge.cov < 1.5 * avg_cov or edge.len > 40000):
                        if edge.id in nonunique:
                            sys.stdout.info("Edge removed based on alignment to alternative:", edge.id, edge.cov, edge.len)
                        else:
                            contigs.add(Contig(edge.seq, edge.id))
                    else:
                        sys.stdout.info("Edge removed based on coverage:", edge.id, edge.cov, edge.len)
                elif (edge.len > 100000 and edge.cov < 1.5 * avg_cov) or (edge.len > 40000 and 1.3 * avg_cov > edge.cov > 0.7 * avg_cov):
                    contigs.add(Contig(edge.seq, edge.id))
                    sys.stdout.info("Edge added based on length and coverage:", edge.id, edge.cov, edge.len)

    elif force_unique is not None:
        sys.stdout.info("Using forced unique edge set")
        sys.stdout.trace(force_unique)
        contigs = ContigCollection().loadFromFile(contigs_file).filter(lambda contig: contig.id in force_unique)
    else:
        sys.stdout.info("Considering all contigs unique")
        contigs = ContigCollection().loadFromFile(contigs_file)
    # contigs.loadFromFasta(open(contigs_file, "r"), num_names=True)
    # contigs = contigs.filter(lambda contig: contig.id not in nonunique and len(contig) > params.k + 20)
    sys.stdout.info("Created", len(contigs), "initial contigs")
    if not all_unique or force_unique is not None:
        sys.stdout.info("Polishing contigs")
        polished_contigs = polisher.polishMany(reads, list(contigs.unique()))
        contigs = ContigCollection().addAll(polished_contigs)
    else:
        sys.stdout.info("Skipping contig polishing step since manual unique contig initialization was used")
    return contigs


def CreateReadCollection(reads_file, cut_reads, downsample):
    sys.stdout.info("Creating read collection")
    num = params.downsample
    if downsample < 1:
        sys.stdout.info("Downsampling:", downsample)
        reads = ReadCollection()
        reads.loadFromFile(reads_file)
        num = int(reads.__len__() * downsample)
    reads = ReadCollection()
    reads.loadFromFile(reads_file, num, cut_reads)
    return reads

def RelevantReadsFromDump(read_dump, edges, reads):
    # type: (str, ContigStorage, ReadCollection) -> Dict[str, List[AlignedRead]]
    relevant_read_ids = dict()
    for edge in edges:
        relevant_read_ids[edge.id] = []
    for s in open(read_dump, "r").readlines():
        s = s.split()
        if s[0] != "Aln":
            continue
        eid = s[6].split("_")[1]
        rid = s[2][1:]
        if eid in relevant_read_ids:
            relevant_read_ids[eid].append(reads[rid])
    return relevant_read_ids

def ExtendShortContigs(contigs, reads, aligner, polisher, read_dump):
    # type: (ContigStorage, ReadCollection, Aligner, Polisher, str) -> None
    sys.stdout.info("Extending short lines")
    short_contigs = ContigStorage()
    als = dict() # type: Dict[str, List[AlignmentPiece]]
    for contig in contigs.unique():
        if len(contig) < params.k + 500:
            short_contigs.add(contig)
            als[contig.id] = []
            als[contig.rc.id] = []

    if read_dump is not None:
        sys.stdout.trace("Using flye read dump file to extend short contigs")
        relevant_reads = RelevantReadsFromDump(read_dump, short_contigs, reads)
        for contig in short_contigs:
            for al in aligner.overlapAlign(relevant_reads[contig.id], ContigStorage([contig])):
                als[al.seg_to.contig.id].append(al)
                als[al.seg_to.contig.rc.id].append(al.rc)
    else:
        sys.stdout.trace("Realigning all reads to extend short contigs")
        for al in aligner.overlapAlign(reads, short_contigs):
            if al.seg_to.left <= 20 and al.rc.seg_to.left <= 20:
                added = False
                for i, al1 in enumerate(als[al.seg_to.contig.id]):
                    if al1.seg_from.contig.id == al.seg_from.contig.id:
                        added = True
                        if al.percentIdentity() > al1.percentIdentity():
                            als[al.seg_to.contig.id][i] = al
                            als[al.seg_to.contig.rc.id][i] = al.rc
                        break
                if not added:
                    als[al.seg_to.contig.id].append(al)
                    als[al.seg_to.contig.rc.id].append(al.rc)
    for contig in short_contigs.unique():
        if len(als[contig.id]) > 0:
            tmp_contig, new_als = polisher.polishEnd(als[contig.id], params.reliable_coverage, max_extension=params.l - len(contig))
            r = len(tmp_contig) - len(contig)
            tmp_contig, new_als = polisher.polishEnd([al.rc for al in new_als], params.reliable_coverage, max_extension=params.l - len(contig))
            l = len(tmp_contig) - len(contig) - r
        else:
            tmp_contig, new_als = contig, als[contig.id]
            l = 0
            r = 0
#        if l > params.k / 2 and r > params.k / 2:
#            tmp_contig.seq = tmp_contig.seq[l - params.k / 2:-r + params.k / 2]
#        else:
#            tmp_contig.seq = tmp_contig.seq[max(0, l - params.k):-max(1, r - params.k)]
        if len(tmp_contig) > params.k + 500:
            sys.stdout.info("Prolonged contig", contig.id, "for", l, "and", r, "nucleotides from left and right")
            contigs.add(Contig(tmp_contig.rc.seq, contig.id))
        else:
            sys.stdout.warn("Could not prolong contig", contig.id, "enough. Removing it.")
            contigs.remove(contig)


def constructDisjointigs(reads, total_length, dir):
    # type: (ReadCollection, int, str) -> str
    basic.ensure_dir_existance(dir)
    reads_file = os.path.join(dir, "reads.fasta")
    disjointigs_file = os.path.join(dir, "disjointigs.fasta")
    log_file = os.path.join(dir, "log.txt")
    reads.print_fasta(open(reads_file, "w"))
    subprocess.check_call([os.path.join(params.bin_path, "flye-modules"), "assemble", "--reads", reads_file, "--out-asm", disjointigs_file, "--genome-size", str(total_length),
                           "--config", "flye/config/bin_cfg/asm_raw_reads.cfg", "--min-ovlp", "1500", "--threads", str(params.threads), "--log", log_file])
    return disjointigs_file
