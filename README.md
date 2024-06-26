# SMFTools

Tools for working with Nanopore alignments with methylation calls, particularly for single molecule footprinting with 6mA. For use with https://github.com/owensnick/GenomeFragments.jl.

## Installation
```julia
    ] # to enter Pkg mode
    add https://github.com/owensnick/SMFTools.jl
```

Convert BAM file to a binary file suitable for memory mapping:

```julia
    streambam_smf_pos("in.bam", "out.posbin")
```

Load memory mapped methylation file:
```julia
    PM = load_pos_matrix("out.posbin")
```

Calculate methylation summary statistics from `meth.bam` aligned to `genome.fa` and save in `meth_stats.tsv.gz`:

```julia
    bam_meth_stats("meth.bam", "genome.fa", "meth_stats.tsv.gz")
```
