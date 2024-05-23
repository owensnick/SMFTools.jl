module SMFTools

using DataFrames, CSV
using BioAlignments, BioSequences, XAM, GenomicFeatures, GenomeFragments, FASTX
using ProgressMeter
using Statistics, StatsBase, OnlineStats
using Mmap, UnalignedVectors
using Parsers

using Blosc

export streambam_smf, validfrag, streamfrags, posmlp, bam_meth_stats, load_pos_matrix, getreadindex, getposview, modmeta, modheat

include("stream.jl")
include("readstats.jl")
include("posmatrix.jl")
include("processing.jl")
end
