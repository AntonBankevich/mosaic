import os
import logging
import json

import sys
sys.path.append("py")

import alignment as aln
import bubbles as bbl
import polish as pol
from flye_tools import config

logger = logging.getLogger()
logger.setLevel(logging.ERROR)
class Job(object):
    """
    Describes an abstract list of jobs with persistent
    status that can be resumed
    """
    run_description = {"stage_name" : ""}

    def __init__(self):
        self.name = None
        self.args = None
        self.work_dir = None
        self.out_files = {}
        self.log_file = None

    def run(self):
        pass

    def save(self, save_file):
        Job.run_description["stage_name"] = self.name

        with open(save_file, "w") as fp:
            json.dump(Job.run_description, fp)

    def load(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)
            Job.run_description = data

    def completed(self, save_file):
        with open(save_file, "r") as fp:
            data = json.load(fp)

            for file in self.out_files.values():
                if not os.path.exists(file):
                    return False

            return True

class JobPolishing(Job):
    # CHANGED SO THAT READS ARE NOT IN ARGS
    # Also changed the directory to an input
    def __init__(self, args, work_dir, log_file, reads, in_contigs, polish_dir):
        super(JobPolishing, self).__init__()

        self.args = args
        self.log_file = log_file
        self.in_contigs = in_contigs
        self.reads = reads
        self.polishing_dir = os.path.join(work_dir, polish_dir)

        self.name = "polishing"
        final_contigs = os.path.join(self.polishing_dir,
                                     "polished_{0}.fasta".format(args.num_iters))
        self.out_files["contigs"] = final_contigs
        self.out_files["stats"] = os.path.join(self.polishing_dir,
                                               "contigs_stats.txt")

    def run(self):
        if not os.path.isdir(self.polishing_dir):
            os.mkdir(self.polishing_dir)

        prev_assembly = self.in_contigs
        contig_lengths = None
        for i in xrange(self.args.num_iters):
            logger.info("Polishing genome ({0}/{1})".format(i + 1,
                                                            self.args.num_iters))

            alignment_file = os.path.join(self.polishing_dir,
                                          "minimap_{0}.sam".format(i + 1))
            logger.info("Running Minimap2")
            aln.make_alignment(prev_assembly, self.reads, self.args.threads,
                               self.polishing_dir, self.args.platform,
                               alignment_file)

            logger.info("Separating alignment into bubbles")

            contigs_info = aln.get_contigs_info(prev_assembly)
            bubbles_file = os.path.join(self.polishing_dir,
                                        "bubbles_{0}.fasta".format(i + 1))
            coverage_stats = \
                bbl.make_bubbles(alignment_file, contigs_info, prev_assembly,
                                 self.args.platform, self.args.threads,
                                 config.vals["min_aln_rate"], bubbles_file)

            logger.info("Correcting bubbles")
            polished_file = os.path.join(self.polishing_dir,
                                         "polished_{0}.fasta".format(i + 1))
            contig_lengths = pol.polish(bubbles_file, self.args.threads,
                                        self.args.platform, self.polishing_dir,
                                        i + 1, polished_file)
            prev_assembly = polished_file

        with open(self.out_files["stats"], "w") as f:
            f.write("seq_name\tlength\tcoverage\n")
            for ctg_id in contig_lengths:
                f.write("{0}\t{1}\t{2}\n".format(ctg_id,
                                                 contig_lengths[ctg_id], coverage_stats[ctg_id]))

class Args:
    def __init__(self):
        self.num_iters = 2
        self.platform = "pacbio"
        self.threads = 8


def polish_from_disk(dir, contigs_file, reads_file):
    args = Args()
    if not os.path.exists(dir):
        os.makedirs(dir)
    if not os.path.exists(os.path.join(dir, "work")):
        os.makedirs(os.path.join(dir, "work"))

    job = JobPolishing(args, os.path.join(dir, "work"), os.path.join(dir, "log.info"), [reads_file], contigs_file,
                             "polish")
    job.run()
    return job.out_files["contigs"]

if __name__ == "__main__":
    polish_from_disk(sys.argv[1], sys.argv[2], sys.argv[3])