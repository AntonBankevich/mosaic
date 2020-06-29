#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs assemble binary
"""

from __future__ import absolute_import
import subprocess
import logging
import os

from flye.utils.utils import which

ASSEMBLE_BIN = "flye-modules"
logger = logging.getLogger()


class AssembleException(Exception):
    pass


def check_binaries():
    if not which(ASSEMBLE_BIN):
        raise AssembleException("Assemble binary was not found. "
                                "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([ASSEMBLE_BIN, "assemble", "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))



def assemble(args, run_params, out_file, log_file, config_path):
    logger.info("Assembling disjointigs")
    logger.debug("-----Begin assembly log------")
    cmdline = [ASSEMBLE_BIN, "assemble", "--reads", ",".join(args.reads), "--out-asm", out_file,
               "--genome-size", str(args.genome_size), "--config", config_path,
               "--log", log_file, "--threads", str(args.threads)]
    if args.debug:
        cmdline.append("--debug")
    if args.meta:
        cmdline.append("--meta")
    cmdline.extend(["--min-ovlp", str(run_params["min_overlap"])])
    cmdline.extend(["--kmer", str(run_params["kmer_size"])])
    if run_params["min_read_length"] > 0:
        cmdline.extend(["--min-read", str(run_params["min_read_length"])])
    #if args.min_kmer_count is not None:
    #    cmdline.extend(["-m", str(args.min_kmer_count)])
    #if args.max_kmer_count is not None:
    #    cmdline.extend(["-x", str(args.max_kmer_count)])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))
