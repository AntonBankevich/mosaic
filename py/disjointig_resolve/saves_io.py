import sys

from typing import Tuple

from alignment.align_tools import Aligner
from common import params
from common.save_load import TokenReader, TokenWriter
from common.sequences import ContigCollection
from common.alignment_storage import ReadCollection
from disjointig_resolve.initialization import CreateReadCollection
from disjointig_resolve.line_storage import NewLineStorage
from disjointig_resolve.disjointigs import DisjointigCollection
from disjointig_resolve.dot_plot import LineDotPlot
from disjointig_resolve.cl_params import Params


def loadAll(handler):
    # type: (TokenReader) -> Tuple[Params, Aligner, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage, LineDotPlot]
    cl_params = Params()
    cl_params.load(handler)
    aligner = Aligner.load(handler)
    sys.stdout.info("Loading contigs")
    contigs = ContigCollection()
    contigs.load(handler)
    sys.stdout.info("Loading reads")
    reads = CreateReadCollection(cl_params.reads_file, cl_params.downsample)
    reads.loadFromFasta(open(cl_params.reads_file, "r"), downsample=params.downsample)
    tmp_reads = reads.copy().addAllRC()
    sys.stdout.info("Loading disjointigs")
    disjointigs = DisjointigCollection()
    disjointigs.load(handler, tmp_reads)
    sys.stdout.info("Loading lines")
    lines = NewLineStorage(disjointigs, aligner)
    lines.load(handler, tmp_reads, contigs)
    sys.stdout.info("Loading dot plot")
    dot_plot = LineDotPlot(lines, aligner)
    dot_plot.load(handler)
    sys.stdout.info("Loading finished")
    return cl_params, aligner, contigs, reads, disjointigs, lines, dot_plot


def saveAll(handler, params, aligner, contigs, reads, disjointigs, lines, dot_plot):
    # type: (TokenWriter, Params, Aligner, ContigCollection, ReadCollection, DisjointigCollection, NewLineStorage, LineDotPlot) -> None
    params.save(handler)
    aligner.save(handler)
    contigs.save(handler)
    disjointigs.save(handler)
    lines.save(handler)
    dot_plot.save(handler)