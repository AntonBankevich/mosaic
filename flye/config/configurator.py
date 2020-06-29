#(c) 2013-2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Sets up some parameters for the run based on input
"""

from __future__ import absolute_import
from __future__ import division
import logging

import flye.utils.fasta_parser as fp
import flye.config.py_cfg as cfg
from flye.six import iteritems


logger = logging.getLogger()


def setup_params(args):
    logger.info("Configuring run")
    parameters = {}
    parameters["pipeline_version"] = cfg.vals["pipeline_version"]

    total_length = 0
    read_lengths = []
    for read_file in args.reads:
        for _, seq_len in iteritems(fp.read_sequence_lengths(read_file)):
            total_length += seq_len
            read_lengths.append(seq_len)

    _, reads_n50 = _calc_nx(read_lengths, total_length, 0.50)
    _, reads_n90 = _calc_nx(read_lengths, total_length, 0.90)

    #Selecting minimum overlap
    logger.info("Total read length: %d", total_length)

    coverage = total_length // args.genome_size
    logger.info("Input genome size: %d", args.genome_size)
    logger.info("Estimated coverage: %d", coverage)
    if coverage < 5 or coverage > 1000:
        logger.warning("Expected read coverage is " + str(coverage) +
                       ", the assembly is not " +
                       "guaranteed to be optimal in this setting." +
                       " Are you sure that the genome size " +
                       "was entered correctly?")

    logger.info("Reads N50/N90: %d / %d", reads_n50, reads_n90)
    if args.min_overlap is None:
        GRADE = 1000
        int_min_ovlp = int(round(reads_n90 / GRADE)) * GRADE

        MIN_OVLP = cfg.vals["min_overlap_range"][args.read_type][0]
        MAX_OVLP = cfg.vals["min_overlap_range"][args.read_type][1]
        if args.meta:
            MAX_OVLP = min(MAX_OVLP, cfg.vals["max_meta_overlap"])
        parameters["min_overlap"] = max(MIN_OVLP, min(MAX_OVLP, int_min_ovlp))
        logger.info("Minimum overlap set to %d", parameters["min_overlap"])
    else:
        parameters["min_overlap"] = args.min_overlap
        logger.info("Selected minimum overlap: %d", parameters["min_overlap"])

    #Selecting k-mer size
    if args.genome_size < cfg.vals["big_genome_kmer"]:
        parameters["kmer_size"] = cfg.vals["kmer_size"][args.read_type][0]
    else:
        parameters["kmer_size"] = cfg.vals["kmer_size"][args.read_type][1]
    logger.info("Selected k-mer size: %d", parameters["kmer_size"])

    #Downsampling reads for the first assembly stage to save memory
    target_cov = None
    if args.asm_coverage and args.asm_coverage < coverage:
        target_cov = args.asm_coverage
    #if not args.asm_coverage and args.genome_size >= 10 ** 9:
    #    target_cov = cfg.vals["reduced_asm_cov"]

    if target_cov:
        logger.info("Using longest %dx reads for contig assembly", target_cov)
        min_read = _get_downsample_threshold(read_lengths,
                                             args.genome_size * target_cov)
        logger.debug("Min read length cutoff: %d", min_read)
        parameters["min_read_length"] = min_read
    else:
        parameters["min_read_length"] = 0

    return parameters


def _calc_nx(scaffolds_lengths, assembly_len, rate):
    n50 = 0
    sum_len = 0
    l50 = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        l50 += 1
        if sum_len > rate * assembly_len:
            n50 = l
            break
    return l50, n50


def _get_downsample_threshold(read_lengths, target_len):
    sum_len = 0
    for l in sorted(read_lengths, reverse=True):
        sum_len += l
        if sum_len > target_len:
            return l

    return 0
