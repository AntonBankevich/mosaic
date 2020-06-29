#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Quick and dirty alignment consensus
"""

import logging
from collections import defaultdict
from itertools import izip
import multiprocessing
import signal

from flye_tools.alignment import shift_gaps, SynchronizedSamReader
import flye_tools.config as config
import flye_tools.fasta_parser as fp

logger = logging.getLogger()

class Profile:
    __slots__ = ("insertions", "matches", "nucl")

    def __init__(self):
        self.insertions = defaultdict(str)
        self.matches = defaultdict(int)
        self.nucl = "-"

def _thread_worker(aln_reader, contigs_info, platform, results_queue,
                   error_queue):
    try:
        aln_reader.init_reading()

        while not aln_reader.is_eof():
            ctg_id, ctg_aln = aln_reader.get_chunk()
            if ctg_id is None:
                break

            profile, aln_errors = _contig_profile(ctg_aln, platform,
                                                  contigs_info[ctg_id].length)
            sequence = _flatten_profile(profile)
            results_queue.put((ctg_id, sequence, aln_errors))

    except Exception as e:
        error_queue.put(e)


def get_consensus(alignment_path, contigs_path, contigs_info, num_proc,
                  platform, min_aln_length):
    """
    Main function
    """
    aln_reader = SynchronizedSamReader(alignment_path,
                                       fp.read_fasta_dict(contigs_path),
                                       min_aln_length)
    manager = multiprocessing.Manager()
    results_queue = manager.Queue()
    error_queue = manager.Queue()

    #making sure the main process catches SIGINT
    orig_sigint = signal.signal(signal.SIGINT, signal.SIG_IGN)
    threads = []
    for _ in xrange(num_proc):
        threads.append(multiprocessing.Process(target=_thread_worker,
                                               args=(aln_reader, contigs_info,
                                                     platform, results_queue,
                                                     error_queue)))
    signal.signal(signal.SIGINT, orig_sigint)

    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        for t in threads:
            t.terminate()

    if not error_queue.empty():
        raise error_queue.get()

    out_fasta = {}
    total_aln_errors = []
    while not results_queue.empty():
        ctg_id, ctg_seq, aln_errors = results_queue.get()
        total_aln_errors.extend(aln_errors)
        if len(ctg_seq) > 0:
            out_fasta[ctg_id] = ctg_seq

    mean_aln_error = float(sum(total_aln_errors)) / (len(total_aln_errors) + 1)
    logger.info("Alignment error rate: {0}".format(mean_aln_error))

    return out_fasta


def _contig_profile(alignment, platform, genome_len):
    """
    Computes alignment profile
    """
    max_aln_err = config.vals["err_modes"][platform]["max_aln_error"]
    aln_errors = []
    profile = [Profile() for _ in xrange(genome_len)]
    for aln in alignment:
        #if aln.err_rate > max_aln_err: continue
        aln_errors.append(aln.err_rate)

        #after gap shifting it is possible that
        #two gaps are aligned against each other
        qry_seq = shift_gaps(aln.trg_seq, aln.qry_seq)
        trg_seq = shift_gaps(qry_seq, aln.trg_seq)

        trg_pos = aln.trg_start
        for trg_nuc, qry_nuc in izip(trg_seq, qry_seq):
            if trg_nuc == "-":
                trg_pos -= 1
            if trg_pos >= genome_len:
                trg_pos -= genome_len

            prof_elem = profile[trg_pos]
            if trg_nuc == "-" and qry_nuc != "-":
                prof_elem.insertions[aln.qry_id] += qry_nuc
            else:
                prof_elem.nucl = trg_nuc
                prof_elem.matches[qry_nuc] += 1

            trg_pos += 1

    return profile, aln_errors


def _flatten_profile(profile):
    growing_seq = []
    ins_group = defaultdict(int)

    for elem in profile:
        pos_matches = elem.matches
        pos_insertions = elem.insertions
        pos_nucl = elem.nucl

        ins_group.clear()
        for ins_str in pos_insertions.values():
            ins_group[ins_str] += 1

        coverage = sum(pos_matches.values())

        max_match = pos_nucl
        if len(pos_matches):
            max_match = max(pos_matches, key=pos_matches.get)
        max_insert = None
        if ins_group:
            max_insert = max(ins_group, key=ins_group.get)

        if max_match != "-":
            growing_seq.append(max_match)
        if max_insert and max_insert != "-" and ins_group[max_insert] > coverage / 2:
            growing_seq.append(max_insert)

    return "".join(growing_seq)
