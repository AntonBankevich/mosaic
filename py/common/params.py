from common.scoring_model import SimpleScores

MINIMAP_BIN = "flye-minimap2"

bin_path = "bin"
clean = False
reliable_coverage = 9
k_cov = reliable_coverage
l_cov = 5
ss_for_kl_adjustment = 20000
maxCoverageThreshold = 20
uncoveredFractionForK = 0.001
radius = 8
alignment_correction_radius = 30
alignment_smoothing_radius = 100
score_counting_radius = 10
full_stats = False
num_iters = 2
min_reads_in_knot = 2
min_expected_Pacbio_PI = 0.78
min_allowed_Pacbio_PI = 0.7
max_jump = 6000
assert_pi = False
min_pi = 0.35
max_allowed_unaligned = 200
min_isolated_length = 15000
k = 1500
l = 2500
bad_end_length = 500
# k = 500
# l = 1500
threads = 16
min_k_mer_cov = 10
min_contra_for_break = 5 # minimal nunber of unaligned reads needed to break resolved segment
max_read_length = 80000
min_alignment_size = 100 # minimal size of alignment that will be considered. Important for composition of alignments.
technology = "pacbio"
expected_size = 5000000
strict_merging_alignment = True
save_alignments = False

scores = SimpleScores(scoreIns = 8, scoreDel = 7, scoreMM = 10, scoreHomo = 4)


downsample = 100000000


flanking_size = 500
window_size = 1500


