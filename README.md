# SMFTools

Tools for working with Nanopore alignments with methylation calls, particularly for single molecule footprinting with 6mA. For use with https://github.com/owensnick/GenomeFragments.jl.


Convert BAM file to a binary file suitable for memory mapping

```julia
    streambam("in.bam", "out.bin")
```

