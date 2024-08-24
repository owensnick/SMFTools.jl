### functions for steaming modkit extract files


### record for a modifcation
struct Mod
    forward_read_pos::Int
    ref_pos::Int
    mod_qual::Float64
    mod_code::String1
    base_qual::Int
    query_kmer::String7
    ref_kmer::String7
end

 

### record for one line for modkit file
struct ModLine{T}
    read_id::T
    forward_read_position::Int
    ref_position::Int
    chrom::String31
    strand::String1
    fw_soft_clipped_start::Int
    fw_soft_clipped_end::Int
    read_length::Int
    mod_qual::Float64
    mod_code::String1
    base_qual::Int
    query_kmer::String7
    ref_kmer::String7
    canonical_base::String1
    flag::Int
end

"""
    parse_mod_line(line)

    Function to economically loop over modkit extract lines and parse information
"""
function parse_mod_line(line)
    read_id = ""
    forward_read_position = -1
    ref_position = -1
    chrom = String31("chrN")
    strand = String1(".")
    fw_soft_clipped_start = 0
    fw_soft_clipped_end = 0
    read_length = 0
    mod_qual = 0.0
    mod_code = String1("z")
    base_qual = 0
    query_kmer = String7("-----")
    ref_kmer = String7("-----")
    canonical_base = String1("N")
    flag = 0
    
    for (i, f) in enumerate(eachsplit(line, '\t'))
        if i == 0
        elseif i == 1   #"read_id"
            read_id = f
        elseif i == 2   #"forward_read_position"
            # @show f
            forward_read_position = Parsers.parse(Int, f)
        elseif i == 3   #"ref_position"
            ref_position = Parsers.parse(Int, f)
        elseif i == 4   #"chrom"
            chrom = String31(f)
        elseif i == 5   #"mod_strand"
            ## always +
            @assert f == "+"
        elseif i == 6   #"ref_strand"
            strand = String1(f) ## read strand
        elseif i == 7   #"ref_mod_strand"
            ## same as ref_Strand
        elseif i == 8   #"fw_soft_clipped_start"
            fw_soft_clipped_start = Parsers.parse(Int, f)
        elseif i == 9   #"fw_soft_clipped_end"
            fw_soft_clipped_end = Parsers.parse(Int, f)
        elseif i == 10   #"read_length"
            #read_length = Parsers.parse(Int, f) # commented out as not used 
        elseif i == 11   #"mod_qual"
            mod_qual = Parsers.parse(Float64, f)
        elseif i == 12   #"mod_code"
            mod_code = String1(f)
        elseif i == 13   #"base_qual"
            base_qual = Parsers.parse(Int, f)
        elseif i == 14   #"ref_kmer"
            ## always .
            if f != "."
                ref_kmer = String7(f)
            end
        elseif i == 15   #"query_kmer"
            query_kmer = String7(f)
        elseif i == 16   #"canonical_base"
            canonical_base = String1(f)
        elseif i == 17   #"modified_primary_base"
            ## same as canonical base
        elseif i == 18   #"inferred"
            ## always false
        elseif i == 19   #"flag"
            # flag = Parsers.parse(Int, f) ### commented out as not used
        end
    end


    ModLine(read_id,
    forward_read_position,
    ref_position,
    chrom,
    strand,
    fw_soft_clipped_start,
    fw_soft_clipped_end,
    read_length,
    mod_qual,
    extended_mod_code(mod_code, query_kmer),
    base_qual,
    query_kmer,
    ref_kmer,
    canonical_base,
    flag)
    
end

struct ReadModStats{T}
    read_id::T
    chrom::String31
    fw_soft_clipped_start::Int
    fw_soft_clipped_end::Int
end



"""
   extended_mod_code(modcode::T, kmer) 

   Extend modcode definitions from a, m, h

   See modcodedict for definitions
"""
function extended_mod_code(modcode::T, kmer) where T ### assumes modcode is a, m, h

    if modcode == "a"
        return modcode
    elseif modcode == "m"
        central_gc = kmer[2] == 'G' # NGCNN
        central_cg = kmer[4] == 'G' # NNCGN

        if central_gc && central_cg
            return T("x")
        elseif central_gc
            return T("i")
        elseif central_cg
            return modcode
        else
            return T("n")
        end
    elseif modcode == "h"
        modcode == "h"
        central_gc = kmer[2] == 'G' # NGCNN
        central_cg = kmer[4] == 'C' # NNCGN

        if central_gc && central_cg
            return T("z")
        elseif central_gc
            return T("g")
        elseif central_cg
            return modcode
        else
            return T("y")
        end

    else
        return modcode
    end
end

function modcodedict()
    Dict("a" => "6mA", "m" => "5mCG", "h" => "5hmCG",
        "x" => "ambig_5m_CG_GC", "i" => "5mGC", "n" => "5mC_noG",
        "z" => "ambig_5gm_CG_GC", "g" => "5hmGC", "h" => "5hmC_noG")
end


### type that allows different configurations of modifications
struct ModConfig{N}
    mods::NTuple{N, String1}
    modindex::Function
    validmodfun::Function
    ModConfig(nt::NTuple{N, String}, modindex::Function, validmodfun::Function) where {N} = new{N}(String1.(nt), modindex, validmodfun)
end
#### type that sets up ami detection
### useful function to process that the mod codes we are currently interested in are a, m, i
validmod_ami(m) = (m == "a") || (m == "m") || (m == "i")
function modcode_index_ami(modcode)
    if modcode == String1("a")
        return 1
    elseif modcode == String1("m")
        return 2
    elseif modcode == String1("i")
        return 3
    else
        error("unsupported modcode: $(modcode)")
    end
end
config_ami() = ModConfig(("a", "m", "i"), modcode_index_ami, validmod_ami)


#######################################################
## stats funs
### mismatchcode("AAAAA", "TAATT") |> x -> digits(x, base=2, pad=5) = [1 0 0 1 1]
function mismatchcode(kmerA, kmerB)
    collect(kmerA) .== collect(kmerB)
    s = 0
    v = 1
    for (i, (a, b)) in enumerate(zip(kmerA, kmerB))
        
        s += v*(a == b)
        v <<= 1
    end 
    
    2^5 - 1 - s
end

### online stats funs
struct ModOnlineStat{T, V}
    kind::T
    name::String
    fn::Function
    os::V

    ModOnlineStat(kind::T, name, fn, os::V) where {T, V} = new{T, V}(kind, name, fn, os)
end

struct ModQualityHist end
modqualhist() = ModOnlineStat(ModQualityHist(), "ModQuality", m -> methprob_to_uint8(m.mod_qual), CountMap(UInt8))
function statdf(stat::ModOnlineStat{ModQualityHist, V}, modcode) where V
    mod_quals = sort(collect(keys(stat.os.value)))
    counts = getindex.(Ref(stat.os.value), mod_quals)
    DataFrame(ModCode=modcode, Name=stat.name, ModQual=uint8_to_methprob.(mod_quals), Count=counts)
end

struct BaseQualityHist end
basequalhist() = ModOnlineStat(BaseQualityHist(), "BaseQuality", m -> UInt8(m.base_qual), CountMap(UInt8))
function statdf(stat::ModOnlineStat{BaseQualityHist, V}, modcode) where V
    mod_quals = sort(collect(keys(stat.os.value)))
    counts = getindex.(Ref(stat.os.value), mod_quals)
    DataFrame(ModCode=modcode, Name=stat.name, BaseQual=Int.(mod_quals), Count=counts)
end


struct ModBaseQualityHist end
modbasejointhist() = ModOnlineStat(ModBaseQualityHist(), "ModBaseQual", m -> (methprob_to_uint8(m.mod_qual), m.base_qual), HeatMap(0:255, 0:50))
function statdf(stat::ModOnlineStat{ModBaseQualityHist, V}, modcode) where V
    DataFrame(ModCode=modcode, Name=stat.name, BaseQual=repeat(stat.os.yedges[1:end-1], inner=length(stat.os.xedges)-1), ModQualInt=repeat(stat.os.xedges[1:end-1], length(stat.os.yedges)-1), Count=stat.os.counts[:])
end


struct QueryKmerQualityHist end
querykmerqualityhist() = ModOnlineStat(QueryKmerQualityHist(), "QueryKmerQuality", m -> (methprob_to_uint8(m.mod_qual), m.query_kmer), CountMap(Tuple{UInt8, String7}))
function statdf(stat::ModOnlineStat{QueryKmerQualityHist, V}, modcode) where V
    qual_kmer = collect(keys(stat.os.value))
    counts = getindex.(Ref(stat.os.value), qual_kmer)
    df = DataFrame(ModCode=modcode, Name=stat.name, BaseQual=uint8_to_methprob.(first.(qual_kmer)), Kmer=last.(qual_kmer), Count=counts)
    sort!(df, [:BaseQual, :Kmer])
end


struct RefKmerQualityHist end
refkmerqualityhist() = ModOnlineStat(RefKmerQualityHist(), "RefKmerQuality", m -> (methprob_to_uint8(m.mod_qual), m.ref_kmer), CountMap(Tuple{UInt8, String7}))
function statdf(stat::ModOnlineStat{RefKmerQualityHist, V}, modcode) where V
    qual_kmer = collect(keys(stat.os.value))
    counts = getindex.(Ref(stat.os.value), qual_kmer)
    df = DataFrame(ModCode=modcode, Name=stat.name, BaseQual=uint8_to_methprob.(first.(qual_kmer)), Kmer=last.(qual_kmer), Count=counts)
    sort!(df, [:BaseQual, :Kmer])
end

struct KmerMismatchHist end
kmermistmatchhist() = ModOnlineStat(KmerMismatchHist(), "KmerMismatchHist", m -> UInt8(mismatchcode(m.query_kmer, m.ref_kmer)), CountMap(UInt8))
function statdf(stat::ModOnlineStat{KmerMismatchHist, V}, modcode) where V
    kmer_mm_codes = sort(collect(keys(stat.os.value)))
    counts = getindex.(Ref(stat.os.value), kmer_mm_codes)
    DataFrame(ModCode=modcode, Name=stat.name, KmerMMCode=Int.(kmer_mm_codes), Count=counts)
end

# struct ModIntervalHist end
# intervalhist() = BamMethStat(ModIntervalHist(), "Interval Hist", identity, KHist(100))

mod_online_stats() = [modqualhist(), basequalhist(), modbasejointhist(), querykmerqualityhist(), refkmerqualityhist(), kmermistmatchhist()]


### read level stats
struct ModStats
    total::Int
    mismatchrate::Float64
    mean::Float64
    std::Float64
    median::Float64
    min::Float64
    max::Float64
    cor::Float64
    corspearman::Float64
end
modmismatch(mod) = mod.query_kmer[3] == mod.ref_kmer[3]

function modstats(mods, totalmods=length(mods))

    if totalmods == 0
        return ModStats(0, 0, 0, 0, 0, 0, 0, 0, 0)
    end
    
    quals = [mods[i].mod_qual for i = 1:totalmods]
    base_quals = [mods[i].base_qual for i = 1:totalmods]
    mismatches = [modmismatch(mods[i]) for i = 1:totalmods]
    
    ModStats(totalmods, mean(mismatches),
        mean(quals), std(quals), median(quals), minimum(quals), maximum(quals), cor(quals, base_quals), corspearman(quals, base_quals))
end


rename_moddf!(mdf, mod) = rename!(mdf, string.("mod_", mod, "_", names(mdf)))

function modstatdataframe_config(stats, config::ModConfig{N}) where {N}
    fields = fieldnames(ModStats)
    mapreduce(i -> rename_moddf!(DataFrame(Dict(f => [getfield(s, f) for s in stats[i]] for f in fields))[!, collect(fields)], config.mods[i]), hcat, 1:N)
end

function joinmodstat(seqsum, mdf)

    if size(seqsum, 1) == size(mdf, 1)
        [seqsum[!, Not(names(mdf))] mdf]
    else
        seqsum.Index = 1:size(seqsum, 1)
        mdf.Index = 1:size(mdf, 1)
        leftjoin(seqsum[!, Not(names(mdf))], mdf, on=:Index)
    end
end

function writemodstats_config(file, config::ModConfig{N}, readmodstats) where {N}

    path = joinpath(dirname(file), "modstats")
    mkpath(path)

    for i = 1:N
        for ms in readmodstats[i]
            file = joinpath(path, string("modstat_", config.mods[i], "_", ms.name, ".tsv.gz"))
            df = statdf(ms, config.mods[i])
            println("Writing $file")
            CSV.write(file, df, delim='\t', compress=true)
        end
    end
    
end

#######################################################
## main loop

### useful function to set a value of array or push if necessary
addpush!(x, v, i) = (i <= length(v)) ? (v[i] = x) : push!(v, x) 

#### function used to read all the modifications associated with a read
function readmods_config!(io, config::ModConfig{N}, modline=parse_mod_line(readline(io)) ;  mods, totalmods) where {N}
    
    readid = modline.read_id
    chrom = modline.chrom
    fw_soft_clipped_start = modline.fw_soft_clipped_start
    fw_soft_clipped_end = modline.fw_soft_clipped_end

    totalmods .= 0
    
    while (readid == modline.read_id) && !eof(io) 
        if config.validmodfun(modline.mod_code) && (modline.ref_position != -1)
            mod = Mod(modline.forward_read_position, modline.ref_position, modline.mod_qual, modline.mod_code, modline.base_qual, modline.query_kmer, modline.ref_kmer)
            mi = config.modindex(mod.mod_code)
            totalmods[mi] += 1
            addpush!(mod, mods[mi], totalmods[mi])
        end
        modline = parse_mod_line(readline(io))
    end
    readstats = ReadModStats(readid, chrom, fw_soft_clipped_start, fw_soft_clipped_end)
    readstats, modline
end


#### main function for streaming mod extract file
## N gives number of modifications in configuration e.g. N = 3 for a, m, i
function streamblockreadsconfig(file, config::ModConfig{N}=config_ami(); seqsumfile="", seqsum=CSV.read(seqsumfile, DataFrame), nr = -1) where {N}

    ### open file
    fio = open(file)
    io = GzipDecompressorStream(fio)
    readline(io) ### header

    ### remove unmapped reads from sequence summary file
    ind = seqsum.alignment_genome .!= "*"
    seqsum = seqsum[ind, :]
    # seqsum = @subset(seqsum, :alignment_genome .!= "*")
    totalreads = size(seqsum, 1)

    ## set up progress meter
    p = (nr == -1) ? Progress(totalreads) : Progress(nr)

    
    ### setup mod buffers for each modification underconsideration
    mods = [Mod[] for i = 1:N]
    totalmods = zeros(Int, N) ## length of buffer
   
    ### vector of readlevel stats 
    readmodstats = [Vector{ModStats}(undef, totalreads) for i = 1:N]
    ### setup online stats
    global_modstats = [mod_online_stats() for i = 1:N]
    
    ### annotate seqsum with soft clipping details
    seqsum.fw_soft_clipped_start = zeros(Int, totalreads)
    seqsum.fw_soft_clipped_end = zeros(Int, totalreads)
    seqsum.HasModData = falses(totalreads)

    
    n = 0
    modline = parse_mod_line(readline(io))
    while !eof(io)
        readstats, modline = readmods_config!(io, config, modline, mods=mods, totalmods=totalmods)
        (readstats.chrom == ".") && continue ## if read is ummaped then skip

        n += 1
        ### check that read in mods matches with next read_id
        ### this can happen if seqsum.read_id[n] has no modifications
        ### conduct forward sweep
        while (n <= totalreads) && (readstats.read_id != seqsum.read_id[n])
            n += 1
            next!(p)
        end
        if n > totalreads
            display("================================")
            display(readstats)
            display(n)
            display("next line:")
            display(modline)
            # display(seqsum[(n-1):(n+1), :])
            error("Cannot match $(readstats.read_id)")
        end

        seqsum.HasModData[n] = true
        seqsum.fw_soft_clipped_start[n] = readstats.fw_soft_clipped_start
        seqsum.fw_soft_clipped_end[n] = readstats.fw_soft_clipped_end
            
        ### calculate mod stats
        ### loop over mods and build mod stats vector
        for i = 1:N
            # @show modstats(mods[i], totalmods[i])
            readmodstats[i][n] = modstats(mods[i], totalmods[i])
            for j = 1:totalmods[i]
                for stat in global_modstats[i]
                    fit!(stat.os, stat.fn(mods[i][j]))
                end
            end
        end
        
        # (length(totalmods) == 3) && all(x -> totalmods[x] > 100, keys(totalmods)) && return mods
        ### stats of reads to annotate onto the dataframe
               
        next!(p)
        next!(p)
        
        (n == nr) && break
    end
    close(io)
    writemodstats_config(file, config, global_modstats)
    @time "Building Mod Stat DF" mdf = modstatdataframe_config(readmodstats, config)
    @time "Joining onto seq sum" seqsum_mdf = joinmodstat(seqsum, mdf)
    if seqsumfile != ""
        nsfile = replace(seqsumfile, ".tsv.gz" => ".modstats.tsv.gz")
        println("Writing $nsfile...")
        @time "Writing stat file" CSV.write(nsfile, seqsum_mdf, compress=true)
    end
    n, seqsum_mdf
end


methprob_to_int(p, offset=0.5) = Int(round(256.0.*p .+ offset))
methprob_to_uint8(p) = UInt8(methprob_to_int(p) - 1)

uint8_to_methprob(u) = (u + 1 - 0.5)/256