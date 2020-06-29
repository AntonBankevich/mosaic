#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
File with configurations
"""

vals = {
        "raw_cfg" : "asm_raw_reads.cfg",
        "corrected_cfg" : "asm_corrected_reads.cfg",
        "subasm_cfg" : "asm_subasm.cfg",

        "simple_kmer_length" : 4,
        "solid_kmer_length" : 10,
        "max_bubble_length" : 500,
        "max_bubble_branches" : 50,
        "min_aln_rate" : 0.01,

        "scaffold_gap" : 100,

        "err_modes" : {
            "pacbio" : {
                "subs_matrix" : "pacbio_substitutions.mat",
                "hopo_matrix" : "pacbio_homopolymers.mat",
                "solid_missmatch" : 0.2,
                "solid_indel" : 0.2,
                "max_aln_error" : 0.25
            },
            "nano" : {
                "subs_matrix" : "nano_substitutions.mat",
                "hopo_matrix" : "nano_homopolymers.mat",
                "solid_missmatch" : 0.3,
                "solid_indel" : 0.3,
                "max_aln_error" : 0.3
            },
            "pacbio_hi_err" : {
                "subs_matrix" : "p6c4_substitutions.mat",
                "hopo_matrix" : "p6c4_homopolymers.mat",
                "solid_missmatch" : 0.25,
                "solid_indel" : 0.25,
                "max_aln_error" : 0.3
            },
            "subasm" : {
                "subs_matrix" : "pacbio_substitutions.mat",
                "hopo_matrix" : "pacbio_homopolymers.mat",
                "solid_missmatch" : 0.2,
                "solid_indel" : 0.2,
                "max_aln_error" : 0.25
            },
        }
    }
