import itertools
import os
import subprocess
import sys
import time
import traceback

from typing import Iterable


sys.path.append("py")
sys.path.append(".")

from common.sequences import ContigStorage
from disjointig_resolve.dot_plot import LineDotPlot
from common.basic import CreateLog
from disjointig_resolve.accurate_line import NewLine
from disjointig_resolve.smart_storage import AlignmentStorage
from alignment.align_tools import Aligner, DirDistributor
from alignment.polishing import Polisher
from common import basic, SeqIO, params
from disjointig_resolve.knotter import LineMerger
from disjointig_resolve.tests import Tester
from disjointig_resolve.line_extender import LineExtender
from common.save_load import TokenReader, SaveHandler
from common.alignment_storage import ReadCollection, AlignmentPiece
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.saves_io import loadAll, saveAll
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.cl_params import Params, parseFlyeDir
from disjointig_resolve.initialization import CreateLineCollection, CreateDisjointigCollection, CreateContigCollection, \
    CreateReadCollection, constructDisjointigs, LoadLineCollection, adjustKL, ExtendShortContigs


def prepare_disjointigs_file(disjointigs_file, disjointigs_file_list):
    recs = []
    for fn in disjointigs_file_list:
        for rec in SeqIO.parse_fasta(open(fn, "r")):
            recs.append(rec)
    h = open(disjointigs_file, "w")
    for rec in recs:
        SeqIO.write(rec, h, "fasta")
    h.close()


def assemble(args, bin_path):
    params.bin_path = bin_path
    start = time.time()
    cl_params = Params().parse(args)
    ref = ContigStorage()
    if cl_params.test:
        cl_params.reads_file = os.path.dirname(__file__)  + "/../../test_dataset/reads.fasta"
        cl_params.genome_size = 30000
        cl_params.dir = os.path.dirname(__file__)  + "/../../test_results"
        ref.loadFromFile(os.path.dirname(__file__)  + "/../../test_dataset/axbctbdy.fasta", False)
    if cl_params.debug:
        params.save_alignments = True
    cl_params.check()
    CreateLog(cl_params.dir)
    sys.stdout.info("Command line:", " ".join(cl_params.args))
    sys.stdout.info("Started")
    if cl_params.debug:
        sys.stdout.info("Version:", subprocess.check_output(["git", "rev-parse", "HEAD"]))
        sys.stdout.info("Modifications:")
        print subprocess.check_output(["git", "diff"])
    sys.stdout.info("Preparing initial state")
    if cl_params.debug:
        save_handler = SaveHandler(os.path.join(cl_params.dir, "saves"))
    else:
        save_handler = None
    if cl_params.load_from is not None:
        # tmp = cl_params.focus
        sys.stdout.info("Loading initial state from saves")
        cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot = loadAll(TokenReader(open(cl_params.load_from, "r")))
        cl_params.parse(args)
        # cl_params.focus = tmp
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        dot_plot.printAll(sys.stdout)
        printState(lines)
    else:
        aligner = Aligner(DirDistributor(cl_params.alignmentDir()))
        polisher = Polisher(aligner, aligner.dir_distributor)

        reads = CreateReadCollection(cl_params.reads_file, cl_params.cut_reads, cl_params.downsample)


        if cl_params.contigs_file is None:
            sys.stdout.info("Running Flye")
            assembly_dir = os.path.join(cl_params.dir, "assembly_initial")
            reads_file = os.path.join(cl_params.dir, "actual_reads.fasta")
            reads.print_fasta(open(reads_file, "w"))
            subprocess.check_call([os.path.join(params.bin_path, "flye"), "--meta", "-o", assembly_dir, "-t", str(cl_params.threads), "--" + params.technology + "-raw", reads_file, "--genome-size", str(cl_params.genome_size), "--min-overlap", str(params.k)])
            cl_params.set_flye_dir(assembly_dir, cl_params.mode)
        elif len(cl_params.disjointigs_file_list) == 0:
            assembly_dir = os.path.join(cl_params.dir, "assembly_initial")
            reads_file = os.path.join(cl_params.dir, "actual_reads.fasta")
            reads.print_fasta(open(reads_file, "w"))
            disjointigs_file = constructDisjointigs(reads, params.expected_size, assembly_dir)
            # graph_file, contigs_file, disjointigs_file, rep_dir, graph_file_after, contigs_file_after = parseFlyeDir(assembly_dir)
            cl_params.disjointigs_file_list.append(disjointigs_file)
            params.min_contra_for_break = 8

        disjointigs = CreateDisjointigCollection(cl_params.disjointigs_file_list, cl_params.dir, aligner, reads)

        all_unique = cl_params.init_file is not None
        contigs = CreateContigCollection(cl_params.graph_file, cl_params.contigs_file, cl_params.min_cov, aligner, polisher, reads, cl_params.force_unique, all_unique)

        if cl_params.autoKL:
            adjustKL(aligner, reads, contigs)

        if cl_params.init_file is None:
            ExtendShortContigs(contigs, reads, aligner, polisher, cl_params.read_dump)
            lines = CreateLineCollection(cl_params.dir, aligner, contigs, disjointigs, reads, cl_params.split)
        else:
            lines = LoadLineCollection(cl_params.dir, cl_params.init_file, aligner, contigs, disjointigs, reads, polisher)

        sys.stdout.info("Constructing dot plot")
        dot_plot = LineDotPlot(lines, aligner)
        dot_plot.construct(aligner)
        # dot_plot.printAll(sys.stdout)

        sys.stdout.info("Updating sequences and resolved segments.")
        knotter = LineMerger(lines, Polisher(aligner, aligner.dir_distributor), dot_plot)
        extender = LineExtender(aligner, knotter, disjointigs, dot_plot)
        extender.updateAllStructures(itertools.chain.from_iterable(line.completely_resolved for line in lines))
        for line in list(lines.unique()): # type: NewLine
            line.completely_resolved.mergeSegments()
            if len(line.completely_resolved) == 0:
                lines.removeLine(line)
        if cl_params.debug:
            sys.stdout.info( "Saving initial state")
            try:
                writer = save_handler.getWriter()
                sys.stdout.info("Save details:", writer.info)
                saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)
            except Exception as e:
                _, _, tb = sys.exc_info()
                sys.stdout.warn("Could not write save")
                traceback.print_tb(tb)
                sys.stdout.INFO( "Message:", e.message)

    sys.stdout.trace( "Disjointig alignments")
    for line in lines:
        sys.stdout.trace( line.disjointig_alignments)
    sys.stdout.info("Starting expanding alignment-consensus loop")

    EACL(aligner, cl_params, contigs, disjointigs, dot_plot, extender, lines, reads, save_handler)

    dot_plot.printAll(sys.stdout)

    sys.stdout.trace( "Final result:")
    lines.printToFasta(open(os.path.join(cl_params.dir, "lines.fasta"), "w"))
    lines.printKnottedToFasta(open(os.path.join(cl_params.dir, "assembly.fasta"), "w"))
    printState(lines)
    sys.stdout.info("Finished")
    secs = int(time.time() - start)
    days = secs / 60 / 60 / 24
    hours = secs / 60 / 60 % 24
    mins = secs / 60 % 60
    sys.stdout.info("Finished in %d days, %d hours, %d minutes" % (days, hours, mins))
    if cl_params.test:
        passed = False
        for al in aligner.dotplotAlign(lines, ref):
            if len(al) > len(al.seg_to.contig) - 3000:
                passed = True
                break
        if passed:
            sys.stdout.info("Test passed")
        else:
            sys.stdout.info("Test failed")


def printState(lines):
    print "Lines:"
    for line in lines:
        print line
    print "Chains:"
    for chain in lines.chains():
        if chain[-1].knot is not None:
            print "->" + "->".join([line.id for line in chain]) + "->"
        else:
            print "->".join([line.id for line in chain])


def EACL(aligner, cl_params, contigs, disjointigs, dot_plot, extender, lines, reads, save_handler):
    cnt = 0
    stop = False
    while not stop:
        stop = True
        if cl_params.focus is None:
            keys = list(lines.items.keys())
        else:
            keys = cl_params.focus
        for line_id in keys:
            if line_id not in lines.items:
                if cl_params.focus is None:
                    continue
                else:
                    for key in lines.items.keys():
                        if str(basic.parseNegativeNumber(basic.parseLineName(key)[-1])) == line_id and lines[key].knot is None:
                            line_id = key
                            break
                    if line_id not in lines.items.keys():
                        continue
            line = lines[line_id]
            sys.stdout.info("Investigating", line)
            extended = extender.processLine(line)
            if extended > 0:
                cnt += 1
                stop = False
            if cnt > 10:
                cnt = 0
                if save_handler is not None:
                    printState(lines)
                    sys.stdout.trace( "Saving current state")
                    writer = save_handler.getWriter()
                    sys.stdout.trace( "Save details:", writer.info)
                    saveAll(writer, cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot)
