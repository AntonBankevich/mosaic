#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs assemble binary
"""

import subprocess
import logging
import os

from flye_tools.utils import which

ASSEMBLE_BIN = "flye-assemble"
logger = logging.getLogger()


class AssembleException(Exception):
    pass


def check_binaries():
    if not which(ASSEMBLE_BIN):
        raise AssembleException("Assemble binary was not found. "
                                "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([ASSEMBLE_BIN, "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))



def assemble(args, out_file, log_file, config_path):
    logger.info("Assembling reads")
    logger.debug("-----Begin assembly log------")
    cmdline = [ASSEMBLE_BIN, "-l", log_file, "-t", str(args.threads)]
    if args.debug:
        cmdline.append("-d")
    if args.min_overlap is not None:
        cmdline.extend(["-v", str(args.min_overlap)])
    #if args.min_kmer_count is not None:
    #    cmdline.extend(["-m", str(args.min_kmer_count)])
    #if args.max_kmer_count is not None:
    #    cmdline.extend(["-x", str(args.max_kmer_count)])
    cmdline.extend([",".join(args.reads), out_file,
                    str(args.genome_size), config_path])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))
