#(c) 2017 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Final output generator
"""

import sys
import logging

import flye_tools.fasta_parser as fp
import flye_tools.config as config

logger = logging.getLogger()


def generate_scaffolds(contigs_file, links_file, out_scaffolds):

    contigs_fasta = fp.read_fasta_dict(contigs_file)
    scaffolds_fasta = {}
    used_contigs = set()

    connections = {}
    with open(links_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            ctg_1, sign_1, ctg_2, sign_2 = line.split("\t")
            if ctg_1 in contigs_fasta and ctg_2 in contigs_fasta:
                connections[sign_1 + ctg_1] = sign_2 + ctg_2
                connections[rc(sign_2) + ctg_2] = rc(sign_1) + ctg_1

    scaffolds_fasta = {}
    scaffolds_seq = {}
    for ctg in contigs_fasta:
        if ctg in used_contigs: continue

        used_contigs.add(ctg)
        scf = ["-" + ctg]
        #extending right
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        for i, ctg in enumerate(scf):
            scf[i] = rc(ctg[0]) + unsigned(ctg)
        scf = scf[::-1]

        #extending left
        while (scf[-1] in connections and
               unsigned(connections[scf[-1]]) not in used_contigs):
            scf.append(connections[scf[-1]])
            used_contigs.add(unsigned(scf[-1]))

        #generating sequence interleaved by Ns
        if len(scf) == 1:
            scaffolds_fasta[unsigned(ctg)] = contigs_fasta[unsigned(ctg)]
            scaffolds_seq[unsigned(ctg)] = scf
        else:
            scf_name = "scaffold_" + unsigned(scf[0]).strip("contig_")
            scaffolds_seq[scf_name] = scf
            scf_seq = []
            for scf_ctg in scf:
                if scf_ctg[0] == "+":
                    scf_seq.append(contigs_fasta[unsigned(scf_ctg)])
                else:
                    scf_seq.append(fp.reverse_complement(
                                    contigs_fasta[unsigned(scf_ctg)]))
            gap = "N" * config.vals["scaffold_gap"]
            scaffolds_fasta[scf_name] = gap.join(scf_seq)

    fp.write_fasta_dict(scaffolds_fasta, out_scaffolds)
    return scaffolds_seq


class SeqStats:
    __slots__ = ("name", "length", "coverage", "circular",
                 "repeat", "mult", "telomere", "graph_path")

    def __init__(self, name="", length="", coverage="", circular="-",
                 repeat="-", mult="1", telomere="none",
                 graph_path=""):
        self.name = name
        self.length = length
        self.coverage = coverage
        self.circular = circular
        self.repeat = repeat
        self.mult = mult
        self.telomere = telomere
        self.graph_path = graph_path

    def print_out(self, handle):
        handle.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n"
                     .format(self.name, self.length, self.coverage,
                             self.circular, self.repeat, self.mult,
                             self.graph_path))


def generate_stats(repeat_file, polished_file, scaffolds, out_stats):
    """
    Compiles information fomr multiple stages
    """
    #contigs_length = {}
    #contigs_coverage = {}
    contigs_stats = {}
    header_line = "seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\tgraph_path"
    for line in open(repeat_file, "r").readlines()[1:]:
        tokens = line.strip().split("\t")
        contigs_stats[tokens[0]] = SeqStats(*tokens)
        #if polished_file is None:
            #contigs_length[tokens[0]] = int(tokens[1])
            #contigs_coverage[tokens[0]] = int(tokens[2])

    if polished_file is not None:
        for line in open(polished_file, "r").readlines()[1:]:
            tokens = line.strip().split("\t")
            contigs_stats[tokens[0]].length = tokens[1]
            contigs_stats[tokens[0]].coverage = tokens[2]

    scaffolds_stats = {}
    for scf, scf_seq in scaffolds.iteritems():
        scaffolds_stats[scf] = SeqStats(scf)
        scf_length = sum(map(lambda c: int(contigs_stats[unsigned(c)].length),
                             scf_seq))
        scf_length += (len(scf_seq) - 1) * config.vals["scaffold_gap"]
        scaffolds_stats[scf].length = str(scf_length)

        scf_cov = _mean(map(lambda c: int(contigs_stats[unsigned(c)].coverage),
                        scf_seq))
        scaffolds_stats[scf].coverage = str(scf_cov)

        scaffolds_stats[scf].repeat = contigs_stats[unsigned(scf_seq[0])].repeat
        scaffolds_stats[scf].circular = contigs_stats[unsigned(scf_seq[0])].circular

        scf_mult = min(map(lambda c: int(contigs_stats[unsigned(c)].mult),
                           scf_seq))
        scaffolds_stats[scf].mult = str(scf_mult)

        #telomere information
        telomere_left = contigs_stats[unsigned(scf_seq[0])].telomere
        telomere_right = contigs_stats[unsigned(scf_seq[-1])].telomere
        if scf_seq[0][0] == "+":
            scf_left = telomere_left in ["left", "both"]
        else:
            scf_left = telomere_left in ["right", "both"]
        if scf_seq[-1][0] == "+":
            scf_right = telomere_right in ["right", "both"]
        else:
            scf_right = telomere_right in ["left", "both"]
        #if scf_left and scf_right: scaffolds_stats[scf].telomere = "both"
        #elif scf_left and not scf_right: scaffolds_stats[scf].telomere = "left"
        #elif not scf_left and scf_right: scaffolds_stats[scf].telomere = "right"
        #else: scaffolds_stats[scf].telomere = "none"

        #graph path
        path = []
        for ctg in scf_seq:
            ctg_path = contigs_stats[unsigned(ctg)].graph_path
            if ctg[0] == "-":
                ctg_path = ",".join(map(lambda x: str(-int(x)),
                                        ctg_path.split(","))[::-1])
            path.append(ctg_path)
        prefix = "*," if scf_left else ""
        suffix = ",*" if scf_right else ""
        scaffolds_stats[scf].graph_path = prefix + ",??,".join(path) + suffix

    with open(out_stats, "w") as f:
        f.write(header_line + "\n")
        for scf in sorted(scaffolds_stats,
                          key=lambda x: int(x.rsplit("_", 1)[-1])):
            scaffolds_stats[scf].print_out(f)

    total_length = sum(map(lambda x: int(x.length), scaffolds_stats.values()))
    if total_length == 0: return

    num_scaffolds = len(scaffolds_stats)
    num_contigs = sum(map(lambda x: len(x), scaffolds.values()))

    scaffold_lengths = map(lambda s: int(s.length), scaffolds_stats.values())
    contig_lengths = []
    for scf in scaffolds.values():
        for ctg in scf:
            contig_lengths.append(int(contigs_stats[unsigned(ctg)].length))
    largest_scf = max(scaffold_lengths)

    ctg_n50 = _calc_n50(contig_lengths, total_length)
    scf_n50 = _calc_n50(scaffold_lengths, total_length)

    mean_read_cov = 0
    for scf in scaffolds_stats.values():
        mean_read_cov += int(scf.length) * int(scf.coverage)
    mean_read_cov /= total_length

    logger.info("Assembly statistics:\n\n"
                "\tTotal length:\t{0}\n"
                "\tContigs:\t{1}\n"
                "\tScaffolds:\t{3}\n"
                #"\tContigs N50:\t{2}\n"
                "\tScaffolds N50:\t{4}\n"
                "\tLargest scf:\t{5}\n"
                "\tMean coverage:\t{6}\n"
                .format(total_length, num_contigs, ctg_n50,
                        num_scaffolds, scf_n50, largest_scf, mean_read_cov))


def rc(sign):
    return "+" if sign == "-" else "-"


def unsigned(ctg):
    return ctg[1:]


def _mean(list):
    if not list: return 0
    return sum(list) / len(list)


def _calc_n50(scaffolds_lengths, assembly_len):
    n50 = 0
    sum_len = 0
    for l in sorted(scaffolds_lengths, reverse=True):
        sum_len += l
        if sum_len > assembly_len / 2:
            n50 = l
            break
    return n50

