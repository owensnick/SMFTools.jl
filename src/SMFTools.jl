module SMFTools

using DataFrames, CSV
using BioAlignments, BioSequences, XAM, GenomicFeatures, GenomeFragments, FASTX
using ProgressMeter
using Statistics, StatsBase

export streambam_smf, validfrag, streamfrags, posmlp, bam_meth_stats

include("stream.jl")
include("readstats.jl")
end
