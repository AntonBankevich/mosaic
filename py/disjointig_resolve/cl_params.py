import getopt
import os
import sys

from typing import List

from common import params
from common.save_load import TokenWriter, TokenReader


class Params:
    def __init__(self):
        self.reads_file = None
        self.contigs_file = None
        self.disjointigs_file_list = []
        self.disjointigs_file = None
        self.load_from = None
        self.graph_file = None
        self.flye_dir = None
        self.dir = None
        self.args = None
        self.threads = 8
        self.test = False
        self.init_file = None
        self.long_params = "test debug nostrict stats genome-size= force-unique= init-file= size= mode= nano cut-reads= homo-score= clean min-cov= nosplit flye-dir= graph= focus= nofocus downsample= output-dir= reads= contigs= disjointigs= load= help".split(" ")
        self.short_params = "o:t:hk:l:"
        self.min_cov = 0
        self.stats = False
        self.new_disjointigs = False
        self.focus = None
        self.split = True
        self.downsample = 1.
        self.mode = "before"
        self.cut_reads = None
        self.autoKL = True
        self.read_dump = None
        self.force_unique = None
        self.genome_size = None
        self.debug = False

    def check(self):
        if self.dir is None:
            self.print_usage_and_exit(1, "Output directory not defined")
        if os.path.exists(self.dir) and os.path.isfile(self.dir):
            self.print_usage_and_exit(1, "Incorrect output directory")
        if self.flye_dir is None and self.genome_size is None and self.contigs_file is None:
            self.print_usage_and_exit(1, "Please define one of the following: directory with Flye assembly, estimated genome size, file with assembled contigs")

    def parse(self, argv):
        # type: (List[str]) -> Params
        self.args = argv
        try:
            options_list, tmp = getopt.gnu_getopt(argv[1:], self.short_params, self.long_params)
            if len(tmp) != 0:
                self.print_usage_and_exit(1, "could not parse parameters")
        except getopt.GetoptError:
            _, exc, _ = sys.exc_info()
            sys.stderr.write(str(exc) + "\n")
            self.print_usage_and_exit(1, "could not parse parameters")
            return
        late_options = []
        for (key, value) in options_list:
            if key == "--output-dir" or key == "-o":
                self.dir = value
                # self.save_dir = os.path.join(self.dir, "saves")
            elif key == "--test":
                self.test = True
            elif key == "--flye-dir":
                self.set_flye_dir(value, self.mode)
            elif key == "--mode":
                self.mode = value
                self.set_flye_dir(self.flye_dir, self.mode)
            elif key == "--debug":
                self.debug = True
            elif key == "--stats":
                self.stats = True
            elif key == "--nostrict":
                params.strict_merging_alignment = False
            elif key == "--init-file":
                self.init_file = value
            elif key == "--nosplit":
                self.split = False
            elif key == "--cut-reads":
                self.cut_reads = int(value)
            elif key == "--homo-score":
                params.scores.sHomo = int(value)
            elif key == "--nano":
                params.technology = "nano"
            elif key == "--size":
                params.expected_size = int(value)
            elif key == "--clean":
                params.clean = True
            elif key == "--reads":
                self.reads_file = value
            elif key == "--min-cov":
                self.min_cov = float(value)
            elif key == "--load":
                self.load_from = value
            elif key == "--force-unique":
                self.force_unique = value.split(",")
            elif key == "--genome-size":
                self.genome_size = int(value)
            elif key == "-t":
                self.threads = int(value)
                params.threads = self.threads
            elif key == "-k":
                params.k = int(value)
                params.l = max(params.l, params.k + 100)
                self.autoKL = False
            elif key == "-l":
                params.l = int(value)
                self.autoKL = False
            elif key == "--focus":
                self.focus = value.split(",")
            elif key == "--nofocus":
                self.focus = None
            elif key == "--downsample":
                self.downsample = float(value)
            elif key == "--help" or key == "-h":
                self.print_usage_and_exit(0)
            else:
                late_options.append((key, value))
        if self.flye_dir is not None:
            self.set_flye_dir(self.flye_dir, self.mode)
        for (key, value) in late_options:
            if key == "--graph":
                self.graph_file = value
            elif key == "--contigs":
                self.contigs_file = value
            elif key == "--disjointigs":
                self.disjointigs_file_list.append(value)
                self.new_disjointigs = True
            else:
                self.print_usage_and_exit(1, "Unknown parameter " + key)
        self.disjointigs_file_list = list(set(self.disjointigs_file_list))
        return self

    def set_flye_dir(self, value, mode):
        self.flye_dir = value
        graph_file, contigs_file, disjointigs_file, rep_dir, read_dump, graph_file_after, contigs_file_after = parseFlyeDir(self.flye_dir)
        if mode == "after":
            graph_file = graph_file_after
            contigs_file = contigs_file_after
        if self.graph_file is None:
            self.graph_file = graph_file
        if self.contigs_file is None:
            self.contigs_file = contigs_file
        self.read_dump = read_dump
        self.disjointigs_file_list.append(disjointigs_file)

    def alignmentDir(self):
        return os.path.join(self.dir, "alignment")

    def saveDir(self):
        return os.path.join(self.dir, "saves")

    def print_usage_and_exit(self, code, message = None):
        if code != 0:
            if message is not None:
                print "Error in command line parameters: " + message
            else:
                print "Error in command line parameters"
        print "usage: mosaic --reads <file> -o <dir> (--genome-size <int> | --flye-dir <dir> | --contigs <contigs>) [--threads <int>]"
        sys.exit(code)

    def save(self, handler):
        # type: (TokenWriter) -> None
        handler.writeTokenLine(self.reads_file)
        handler.writeTokenLine(self.contigs_file)
        handler.writeTokenLine(self.disjointigs_file)
        handler.writeTokenLine(self.load_from)
        handler.writeTokenLine(self.dir)
        handler.writeTokenLine(self.save_dir)
        handler.writeToken(str(self.downsample))

    def load(self, handler):
        # type: (TokenReader) -> None
        self.reads_file = handler.readToken()
        self.contigs_file = handler.readToken()
        self.disjointigs_file = handler.readToken()
        self.load_from = handler.readToken()
        self.dir = handler.readToken()
        self.save_dir = handler.readToken()
        self.downsample = float(handler.readToken())

def parseFlyeDir(flye_dir):
    if "20-repeat" in os.listdir(flye_dir):
        res = os.path.join(flye_dir, "20-repeat", "graph_before_rr.gv"), os.path.join(flye_dir, "20-repeat", "graph_before_rr.fasta"), os.path.join(flye_dir, "10-consensus", "consensus.fasta"), os.path.join(flye_dir, "20-repeat"), os.path.join(flye_dir, "20-repeat", "read_alignment_dump")
    else:
        res = os.path.join(flye_dir, "2-repeat", "graph_before_rr.gv"), os.path.join(flye_dir, "2-repeat", "graph_before_rr.fasta"), os.path.join(flye_dir, "1-consensus", "consensus.fasta"), os.path.join(flye_dir, "2-repeat"), os.path.join(flye_dir, "2-repeat", "read_alignment_dump")
    if not os.path.isfile(os.path.join(flye_dir, "scaffolds.fasta")):
        res = res + (os.path.join(flye_dir, "assembly_graph.gv"), os.path.join(flye_dir, "assembly.fasta"))
    else:
        res = res + (os.path.join(flye_dir, "assembly_graph.gv"), os.path.join(flye_dir, "scaffolds.fasta"))
#       for f in res:
#           assert os.path.isfile(f), f
    return res
