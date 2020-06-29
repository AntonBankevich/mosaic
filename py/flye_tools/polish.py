#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs polishing binary in parallel and concatentes output
"""

import logging
import random
import subprocess
import os
import sys
from collections import defaultdict
from threading import Thread

import flye_tools.bubbles as bbl
import flye_tools.fasta_parser as fp
from common import params
from flye_tools.utils import which
import flye_tools.config as config


POLISH_BIN = "flye-polish"

logger = logging.getLogger()


class PolishException(Exception):
    pass


def check_binaries():
    if not which(POLISH_BIN):
        raise PolishException("polishing binary was not found. "
                              "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([os.path.join(params.bin_path, POLISH_BIN), "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def polish(bubbles_file, num_proc, err_mode, work_dir, iter_id, out_polished):
    _ROOT = os.path.dirname(__file__)

    subs_matrix = os.path.join(_ROOT, "resource",
                               config.vals["err_modes"][err_mode]["subs_matrix"])
    hopo_matrix = os.path.join(_ROOT, "resource",
                               config.vals["err_modes"][err_mode]["hopo_matrix"])

    consensus_out = os.path.join(work_dir, "consensus_{0}.fasta"
                                                .format(iter_id))
    _run_polish_bin(bubbles_file, subs_matrix, hopo_matrix,
                    consensus_out, num_proc)
    polished_fasta, polished_lengths = _compose_sequence([consensus_out])
    fp.write_fasta_dict(polished_fasta, out_polished)

    return polished_lengths


def _run_polish_bin(bubbles_in, subs_matrix, hopo_matrix,
                    consensus_out, num_threads):
    """
    Invokes polishing binary
    """
    cmdline = [os.path.join(params.bin_path, POLISH_BIN), "--threads", str(num_threads), "--bubbles", bubbles_in, "--subs-mat", subs_matrix, "--hopo-mat",
               hopo_matrix, "--out", consensus_out]
    try:
        subprocess.check_call(cmdline, stdout = open("/dev/null", "w"), stderr = open("/dev/null", "w"))
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise PolishException(str(e))
    except OSError as e:
        raise PolishException(str(e))


def _compose_sequence(consensus_files):
    """
    Concatenates bubbles consensuses into genome
    """
    consensuses = defaultdict(list)
    coverage = defaultdict(list)
    for file_name in consensus_files:
        with open(file_name, "r") as f:
            header = True
            for line in f:
                if header:
                    tokens = line.strip().split(" ")
                    ctg_id = tokens[0][1:]
                    ctg_pos = int(tokens[1])
                    coverage[ctg_id].append(int(tokens[2]))
                else:
                    consensuses[ctg_id].append((ctg_pos, line.strip()))
                header = not header

    polished_fasta = {}
    polished_stats = {}
    for ctg_id, seqs in consensuses.iteritems():
        sorted_seqs = map(lambda p: p[1], sorted(seqs, key=lambda p: p[0]))
        concat_seq = "".join(sorted_seqs)
        mean_coverage = sum(coverage[ctg_id]) / len(coverage[ctg_id])
        polished_fasta[ctg_id] = concat_seq
        polished_stats[ctg_id] = len(concat_seq)

    return polished_fasta, polished_stats
