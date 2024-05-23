
using OnlineStats

np_frag_bins() = range(0, 10_1000, length=100)

alignlength(record)  = length(BAM.leftposition(record):BAM.rightposition(record))

getrec(data) = data.record
1+2

abstract type  BAMStat end

struct BamMethStat{T, V}
    kind::T
    name::String
    fn::Function
    os::V

    BamMethStat(kind::T, name, fn, os::V) where {T, V} = new{T, V}(kind, name, fn, os)
end

function plotstat(bamstat::BamMethStat{T, V}) where {T, V}
    error("Not implemented")
end

struct ChromCounter end
chromcounter() = BamMethStat(ChromCounter(), "Chromosome counts", BAM.refname ∘ getrec, CountMap(String))
function plotstat(bamstat::BamMethStat{ChromCounter, V}) where {V}
    h = bamstat.os.value
    k = sort(collect(keys(h))) 
    barplot(1:length(k), getindex.(Ref(h), k), axis=(xticks=(1:length(k), k),))
end

struct MappingQualityHist end
mapqualhist() = BamMethStat(MappingQualityHist(), "Mapping Quality", BAM.mappingquality ∘ getrec, CountMap(UInt8))
function plotstat(bamstat::BamMethStat{MappingQualityHist, V}) where {V}
   k = sort(collect(keys(bamstat.os.value)))
   barplot(k, getindex.(Ref(bamstat.os.value), k), axis=(title=string(bamstat.name, " histogram"), xlabel="Mapping Quality", ylabel="Count"))
end

struct TemplateLengthHist end
templengthhist() = BamMethStat(TemplateLengthHist(), "Read Length Hist", Float64 ∘ BAM.templength ∘ getrec, KHist(100))
function plotstat(bamstat::BamMethStat{TemplateLengthHist, V}) where {V}
    paircount = bamstat.os.bins
    
    f = lines(first.(paircount), last.(paircount), axis=(title=string(bamstat.name, " histogram"), xlabel="Read Length", ylabel="Count"), color=:black)
    band!(first.(paircount), zeros(length(paircount)), last.(paircount), color=(:black, 0.1))
    f
end
struct AlignLengthHist end
alignlengthhist() = BamMethStat(AlignLengthHist(), "Align length hist", Float64 ∘ alignlength ∘ getrec, KHist(100))
function plotstat(bamstat::BamMethStat{AlignLengthHist, V}) where {V}
    paircount = bamstat.os.bins
    
    f = lines(first.(paircount), last.(paircount), axis=(title=string(bamstat.name, " histogram"), xlabel="Align Length", ylabel="Count"), color=:black)
    band!(first.(paircount), zeros(length(paircount)), last.(paircount), color=(:black, 0.1))
    f
end
struct TempAlignJointHist end
tempalignjointhist() = BamMethStat(TempAlignJointHist(), "Frag/Genome heatmap", r -> (Float64(alignlength(r.record)), Float64(BAM.templength(r.record))), HeatMap(np_frag_bins(), np_frag_bins()))
function plotstat(bamstat::BamMethStat{TempAlignJointHist, V}) where {V}
    jh = bamstat.os
    heatmap(jh.xedges, jh.yedges, log10.(jh.counts'), axis=(title=string(bamstat.name, " heatmap"), xlabel="Align Length", ylabel="Read Length"))
    
end
struct MismatchHist end
mismatchhist() = BamMethStat(MismatchHist(), "Mismatch", r -> mean(r.matchtemplate), OnlineStats.Hist(range(0, 1, length=100)))
function plotstat(bamstat::BamMethStat{MismatchHist, V}) where {V}
    h = bamstat.os

    barplot(h.edges[1:end-1], h.counts, axis=(title=string(bamstat.name, " histogram"), xlabel="Mismatch rate", ylabel="Count"))
end
# struct SingleBaseMismatchHist end
# singlebasemmhist() = BamMethStat(SingleBaseMismatchHist(), "Single base Mismatch hist", r -> mean(r.seq_match[r.match_ind]), OnlineStats.Hist(range(0, 1, length=100)))

struct MethylatedAHist end
methylatedahist() = BamMethStat(MethylatedAHist(), "mA meth rate", r -> r.ml_A[r.atind]/256, OnlineStats.Hist(range(0, 1, length=100)))
function plotstat(bamstat::BamMethStat{MethylatedAHist, V}) where {V}
    h = bamstat.os
    barplot(h.edges[1:end-1], log10.(h.counts .+ 1), axis=(title=string(bamstat.name, " histogram"), xlabel="mA rate", ylabel="Count"))
end
struct MethylatedCHist end
methylatedchist() = BamMethStat(MethylatedCHist(), "mC meth rate", r -> r.ml_C[r.gcind]/256, OnlineStats.Hist(range(0, 1, length=100)))
function plotstat(bamstat::BamMethStat{MethylatedCHist, V}) where {V}
    h = bamstat.os
    barplot(h.edges[1:end-1], log10.(h.counts .+ 1), axis=(title=string(bamstat.name, " histogram"), xlabel="mC rate", ylabel="Count"))
end

onlinestats() = (chromcounter(),
                mapqualhist(),
                templengthhist(),
                alignlengthhist(),
                tempalignjointhist(),
                mismatchhist(),
                methylatedahist(),
                methylatedchist())


function char_to_dna(c)
    if c == 'A'
         return DNA_A
     elseif c == 'C'
         return DNA_C
     elseif c == 'G'
         return DNA_G
     elseif c == 'T'
         return DNA_T
     elseif c == 'N'
         return DNA_N
     else
         error("$c not recognised")
     end
 end
 

function correctmd!(md, seq)
    
    mode = 1
    modestart = 1
    seqindex = 1
    for i = 2:length(md)
        
        if mode == 1
            if !isnumeric(md[i])
                seqindex += Parsers.parse(Int, md[modestart:(i-1)])
                if md[i] == '^'
                    mode = 2
                else
                    mode = 3
                    seq[seqindex] = char_to_dna(md[i])
                    seqindex += 1
                end
                modestart = i
            end
        elseif mode == 2
            if isnumeric(md[i])
                for k = (modestart+1):(i - 1)
                    seq[seqindex] = char_to_dna(md[k])
                    seqindex += 1
                end
                
                if isnumeric(md[i])
                    mode = 1
                else
                    mode = 3
                    seq[seqindex] = char_to_dna(md[i])
                    seqindex += 1
                end
                modestart = i
            end
        else mode == 3
            if isnumeric(md[i])
                mode = 1
            elseif md[i] == '^'
                mode = 2
            else
                seq[seqindex] = char_to_dna(md[i])
                seqindex += 1
            end
            modestart = i
        end
        
    end
    if mode == 2
        for k = (modestart+1):length(md)
            seq[seqindex] = char_to_dna(md[k])
            seqindex += 1
        end
     end
    seq
end

"""
    bammatchtemplate(record)

    Calculate the match positions of the read to the reference genome
"""
function bammatchtemplate(record)
    match = falses(BAM.templength(record))
    # aln = BAM.alignment(record)
    # refmatches = ref2seq.(Ref(aln), BAM.leftposition(record):BAM.rightposition(record)); # this costs half the time!
    refindex = refseqindex(record)
    afm = falses(BAM.alignlength(record))
    # @time md = record["MD"]::String
    md = get_md_str(record)
    markmdmatch!(md, afm)
 
    for (i, t) in zip(refindex, afm)
       match[i] = match[i] | t 
    end
    
    match    
end

function get_md_str(record)
    start = BAM.auxdata_position(record)
    stop = BAM.data_size(record)    
    pos = BAM.findauxtag(record.data, start, stop, UInt8('M'), UInt8('D'))
     _, val = BAM.loadauxvalue(record.data, pos+3, String)
    val
end

function bammatchtemplate_old(record)
    match = falses(BAM.templength(record))
    aln = BAM.alignment(record)
    refmatches = ref2seq.(Ref(aln), BAM.leftposition(record):BAM.rightposition(record)); # this costs half the time!
    # refindex = refseqindex(record)
    afm = falses(BAM.alignlength(record))
    md = record["MD"]::String
    markmdmatch_orig!(md, afm)
 
    for ((i, o), t) in zip(refmatches, afm)
       match[i] = match[i] | t 
    end
    
    match    
end


# function refseqindex(record)
#     n = BAM.alignlength(record)
#     r2s = Vector{Int}(undef, n)
#     # rop = Vector{Operation}(undef, n)
#     ops, pos = BAM.cigar_rle(record)
#     tindex = 1
#     rindex = 1
#     for (o, p) in zip(ops, pos)
#         if (o == OP_SOFT_CLIP) || (o == OP_INSERT)
#             tindex += p
#         elseif o == OP_MATCH 
#             @inbounds for i = 1:p
#                 # rop[rindex] = OP_MATCH
#                 r2s[rindex] = tindex
#                 tindex += 1
#                 rindex += 1
#             end
#         elseif o == OP_DELETE
        
#             @inbounds for i = 1:p
#                 # rop[rindex] = OP_DELETE
#                 r2s[rindex] = tindex - 1
#                 rindex += 1
#             end

#         end
#     end
#     # r2s, rop
#     r2s
# end

# function cigar_rle(record::Record, checkCG::Bool = true)::Tuple{Vector{BioAlignments.Operation},Vector{Int}}
#     checkfilled(record)
#     idx, nops = cigar_position(record, checkCG)
#     ops, lens = extract_cigar_rle(record.data, idx, nops)
#     return ops, lens
# end

# function extract_cigar_rle(data::Vector{UInt8}, offset, n)
#     ops = Vector{BioAlignments.Operation}()
#     lens = Vector{Int}()
#     for i in offset:4:offset + (n - 1) * 4
#         x = unsafe_load(Ptr{UInt32}(pointer(data, i)))
#         op = BioAlignments.Operation(x & 0x0F)
#         push!(ops, op)
#         push!(lens, x >> 4)
#     end
#     return ops, lens
# end

"""
    refseqindex(record)

    Calculate the reference positions of the read to the reference genome

    This uses low level access to cigar_rle from XAM.jl to avoid allocating vector of operations and lengths.
"""
function refseqindex(record, checkCG::Bool = true)
    n = BAM.alignlength(record)
    r2s = Vector{Int}(undef, n)
    # rop = Vector{Operation}(undef, n)
    offset, m = BAM.cigar_position(record, checkCG)
    # ops, pos = BAM.cigar_rle(record)
    tindex = 1
    rindex = 1
    for i in offset:4:offset + (m - 1) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        o = BioAlignments.Operation(x & 0x0F)
        p = x >> 4
    # for (o, p) in zip(ops, pos)
        if (o == OP_SOFT_CLIP) || (o == OP_INSERT)
            tindex += p
        elseif o == OP_MATCH 
            for i = 1:p
                # rop[rindex] = OP_MATCH
                r2s[rindex] = tindex
                tindex += 1
                rindex += 1
            end
        elseif o == OP_DELETE
        
            for i = 1:p
                # rop[rindex] = OP_DELETE
                r2s[rindex] = tindex - 1
                rindex += 1
            end

        end
    end
    # r2s, rop
    r2s
end

function seqrefindex(record, checkCG::Bool = true)
    n = BAM.seqlength(record)
    # r2s = Vector{Int}(undef, n)
    s2r = Vector{Int}(undef, n)
    matchind = falses(n)
    # rop = Vector{Operation}(undef, n)
    offset, m = BAM.cigar_position(record, checkCG)
    # ops, pos = BAM.cigar_rle(record)
    tindex = 1
    rindex = 1
    for i in offset:4:offset + (m - 1) * 4
        x = unsafe_load(Ptr{UInt32}(pointer(record.data, i)))
        o = BioAlignments.Operation(x & 0x0F)
        p = x >> 4
    # for (o, p) in zip(ops, pos)
        if (o == OP_SOFT_CLIP) || (o == OP_INSERT)
            tindex += p
        elseif o == OP_MATCH 
            for i = 1:p
                # rop[rindex] = OP_MATCH
                r2s[rindex] = tindex
                tindex += 1
                rindex += 1
            end
        elseif o == OP_DELETE
        
            for i = 1:p
                # rop[rindex] = OP_DELETE
                r2s[rindex] = tindex - 1
                rindex += 1
            end

        end
    end
    # r2s, rop
    r2s
end


function markmdmatch!(md, match)
    
    mode = 1
    modestart = 1
    seqindex = 1
    for i = 2:length(md)
        if mode == 1
            if !isnumeric(md[i])
                # matchlen = Parsers.parse(Int, md[modestart:(i-1)])
                matchlen = Parsers.parse(Int, md, Parsers.OPTIONS, modestart, i-1)
                for i = 1:matchlen
                    match[seqindex] = true
                    seqindex += 1
                end
            
                if md[i] == '^'
                    mode = 2
                else
                    mode = 3
                    match[seqindex] = false
                    seqindex += 1
                end
                modestart = i
 
            end
        elseif mode == 2
            if isnumeric(md[i])
                for k = (modestart+1):(i - 1)
                    match[seqindex] = false
                    seqindex += 1
                end
                mode = 1
                modestart = i
            end    
        else mode == 3
            if isnumeric(md[i])
                mode = 1
            elseif md[i] == '^'
                mode = 2
            else
                match[seqindex] = false
                seqindex += 1
            end
            modestart = i
        end
        
    end
    if mode == 1
        matchlen = Parsers.parse(Int, md[modestart:end])
        for i = 1:matchlen
            match[seqindex] = true
            seqindex += 1
        end

    elseif mode == 2
        for k = (modestart+1):length(md)
            match[seqindex] = false
            seqindex += 1
        end
        
    end
    match
    
end


function markmdmatch_orig!(md, match)
    
    mode = 1
    modestart = 1
    seqindex = 1
    for i = 2:length(md)
        
        if mode == 1
            if !isnumeric(md[i])
                matchlen = Parsers.parse(Int, md[modestart:(i-1)])
                match[seqindex:(seqindex + matchlen - 1)] .= true
                seqindex += matchlen
                
                if md[i] == '^'
                    mode = 2
                else
                    mode = 3
                    match[seqindex] = false
                    seqindex += 1
                end
                modestart = i
            end
        elseif mode == 2
            if isnumeric(md[i])
                for k = (modestart+1):(i - 1)
                    match[seqindex] = false
                    seqindex += 1
                end
                
                if isnumeric(md[i])
                    mode = 1
                else
                    mode = 3
                    match[seqindex] = false
                    seqindex += 1
                end
                modestart = i
            end
        else mode == 3
            if isnumeric(md[i])
                mode = 1
            elseif md[i] == '^'
                mode = 2
            else
                match[seqindex] = false
                seqindex += 1
            end
            modestart = i
        end
        
    end
    if mode == 1
        matchlen = Parsers.parse(Int, md[modestart:end])
        match[seqindex:(seqindex + matchlen - 1)] .= true

    end
    if mode == 2
        for k = (modestart+1):length(md)
            match[seqindex] = false
            seqindex += 1
        end
        
    end
    match
    
end

struct MLIter{T}
    ml::Vector{T}
    start::Int
    stop::Int
    match::BitArray{1}
    baseind::BitArray{1}
    MLIter(ml::Vector{T}, start::Int, stop::Int, match::BitArray{1}, baseind::BitArray{1}) where {T} = new{T}(ml, start, stop, match, baseind)
end
Base.eltype(::Type{MLIter{T}}) where {T} = T
Base.length(iter::MLIter) = sum(iter.match .& iter.baseind)
function Base.iterate(iter::MLIter, state=(iter.start, 1))
    # @show state
    mli, biti = state
    @inbounds while (biti <= length(iter.match)) && (mli <= iter.stop)
        if iter.baseind[biti]
            if iter.match[biti]
                return iter.ml[mli], (mli + 1, biti + 1)
            end
            mli += 1

        end
        biti += 1
        # if iter.match[state] && iter.baseind[state]
        #     return iter.ml[state], state + 1
        # end
        # state += 1
    end
    return nothing
end

function expandrecordtemplate_orig(record)
    matchtemplate = bammatchtemplate(record)
    
    seq = BAM.sequence(record);
    # matchtemplate = falses(length(seq))
    ml = get_ml_data(record)
    # ml = record["ML"] #::Vector{UInt8}

    atind_called = BAM.ispositivestrand(record) ? (seq .== DNA_A) : (seq .== DNA_T) 
    gcind_called = BAM.ispositivestrand(record) ? (seq .== DNA_C) : (seq .== DNA_G) 

    totalat = sum(atind_called)
    # totalgc = sum(gcind_called)

    # ml_A = zeros(UInt8, length(seq))
    # ml_C = zeros(UInt8, length(seq))

    ml_A = Vector{UInt8}(undef, length(seq))
    ml_C = Vector{UInt8}(undef, length(seq))

    ml_A[atind_called] .= ml[1:totalat]
    ml_C[gcind_called] .= ml[(totalat + 1):end]

    atind = atind_called .& matchtemplate
    gcind = gcind_called .& matchtemplate

    # match_ind = matchtemplate
    # seq_match = matchtemplate
    # valid_ml_A = atind
    # valid_ml_C = gcind

    (; record, matchtemplate, ml_A, ml_C, atind, gcind)
    # (; record, match_ind, seq_match, ml_A, ml_C, valid_ml_A, valid_ml_C)

end

function zero_ml_mismatch_slow!(ml, ind, match)
    i = findall(ind)
    mi = findall(.!match)
    ml[i .∈ Ref(Set(mi))] .= 0
    ml
end

function zero_ml_mismatch_fast!(ml, ind, match)
    mli = 0

    for (i, m) in zip(ind, match)
        if i
            mli += 1
            if !m
                ml[mli] = 0
            end
        end
    end
        
    ml
    # i = findall(ind)
    # mi = findall(match)
    # ml[i[mi]] .= 0
    # ml
end


"""

    ind - whether basepair is called as methylated
    sum(ind) = total_ml

"""
function ml_match_ind(total_ml, ind, match)
    match_ml = falses(total_ml)
    mli = 0
    for (i, m) in zip(ind, match)
        if i 
            mli += 1
            m && (match_ml[mli] = true)
        end
    end
    
    match_ml
end

function expandrecordtemplate_new(record)
    matchtemplate = bammatchtemplate(record)
    
    seq = BAM.sequence(record);
    # matchtemplate = falses(length(seq))
    ml = record["Ml"] #::Vector{UInt8}

    atind_called = BAM.ispositivestrand(record) ? (seq .== DNA_A) : (seq .== DNA_T) 
    gcind_called = BAM.ispositivestrand(record) ? (seq .== DNA_C) : (seq .== DNA_G) 

    totalat = sum(atind_called)
    # totalgc = sum(gcind_called)
    

    # ml_A = view(ml, 1:totalat)
    # ml_C = view(ml, (totalat + 1):length(ml))
    # match_ml_A = ml_match_ind(totalat, atind_called, matchtemplate)
    # match_ml_C = ml_match_ind(length(ml) - totalat, gcind_called, matchtemplate)
    
    # ml_A[atind_called] .= ml[1:totalat]
    # ml_C[gcind_called] .= ml[(totalat + 1):end]


    # atind = atind_called .& matchtemplate
    # gcind = gcind_called .& matchtemplate

    # match_ind = matchtemplate
    # seq_match = matchtemplate
    # valid_ml_A = atind
    # valid_ml_C = gcind

    (; record, matchtemplate, ml, totalat, atind_called, gcind_called)
    # (; record, match_ind, seq_match, ml_A, ml_C, valid_ml_A, valid_ml_C)

end


function expandrecordtemplate_ind(record)
    matchtemplate = bammatchtemplate(record)
    
    seq = BAM.sequence(record);
    # matchtemplate = falses(length(seq))
    ml = record["Ml"] #::Vector{UInt8}

    atind_called = BAM.ispositivestrand(record) ? (seq .== DNA_A) : (seq .== DNA_T) 
    gcind_called = BAM.ispositivestrand(record) ? (seq .== DNA_C) : (seq .== DNA_G) 

    totalat = sum(atind_called)
    # totalgc = sum(gcind_called)
    

    ml_A = view(ml, 1:totalat)
    ml_C = view(ml, (totalat + 1):length(ml))

    
    ml_A[atind_called] .= ml[1:totalat]
    ml_C[gcind_called] .= ml[(totalat + 1):end]


    atind = atind_called .& matchtemplate
    gcind = gcind_called .& matchtemplate

    # match_ind = matchtemplate
    # seq_match = matchtemplate
    # valid_ml_A = atind
    # valid_ml_C = gcind

    (; record, matchtemplate, ml_A, ml_C, atind, gcind)
    # (; record, match_ind, seq_match, ml_A, ml_C, valid_ml_A, valid_ml_C)

end

function expandrecord(record, genomereader)

    aln = BAM.alignment(record)
    seq = BAM.sequence(record);
    auxdict = BAM.auxdata(record)
    ml = auxdict["Ml"]::Vector{UInt8}
    refseq = FASTA.extract(genomereader, BAM.refname(record), BAM.leftposition(record):BAM.rightposition(record));


    posmatch = seq2ref.(Ref(aln), 1:length(seq))
    pos = first.(posmatch)
    ops = last.(posmatch)
    relpos = pos .- leftposition(record) .+ 1
    
    match_ind = ops .== Ref(OP_MATCH)
    seq_match = falses(length(match_ind))
    seq_match[match_ind] .= LongDNA{4}(refseq[relpos[match_ind]]) .== seq[match_ind]

    atind = ifelse(BAM.ispositivestrand(record), seq .== DNA_A, seq .== DNA_T)
    gcind = ifelse(BAM.ispositivestrand(record), seq .== DNA_C, seq .== DNA_G)

    ml_A = ml[1:sum(atind)]
    ml_C = ml[sum(atind) .+ (1:sum(gcind))]

    valid_ml_A = findall(atind) .∈ Ref(Set(findall(seq_match)))
    valid_ml_C = findall(gcind) .∈ Ref(Set(findall(seq_match)))

    (; record, pos, ops, match_ind, seq_match, refseq, relpos, ml_A, ml_C, valid_ml_A, valid_ml_C)


end




function bamstatstemplate(bamfile, totalreads=-1; oss = onlinestats(), filtfun= r-> true)
   
    reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
    total_mapped_reads = mapreduce(d -> d[end].n_mapped, +, reader.index.index.data)
    @show total_mapped_reads
    if totalreads == -1
        p = Progress(total_mapped_reads)
    else
        p = Progress(totalreads)
    end
    
    record = BAM.Record()
    reads = 0
    
    while !eof(reader) && reads != totalreads
        read!(reader, record)
        !filtfun(record) && continue
        exrec = expandrecordtemplate_orig(record)
        for bamstat in oss
           fit!(bamstat.os, bamstat.fn(exrec))
        end
        reads += 1
        next!(p)
    end
    
    oss

end

function bamstats(bamfile, genomereader, totalreads=-1; oss = onlinestats(), filtfun= r-> true)

    genomereader = open(FASTA.Reader, genomefile, index=string(genomefile, ".fai"))
    # reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
    
    reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
    total_mapped_reads = mapreduce(d -> d[end].n_mapped, +, reader.index.index.data)
    if totalreads == -1
        p = Progress(total_mapped_reads)
    else
        p = Progress(totalreads)
    end
    
    record = BAM.Record()
    reads = 0

    while !eof(reader) && reads != totalreads
        read!(reader, record)
        !filtfun(record) && continue
        exrec = expandrecord(record, genomereader)
        for bamstat in oss
            fit!(bamstat.os, bamstat.fn(exrec))
        end
        reads += 1
        next!(p)
    end
    
    oss

end


count_chrom!(record, chroms) = push!(chroms, BAM.refname(record))
count_mapqual!(record, mapquals) = push!(mapquals, BAM.mappingquality(record))
count_templength(record, alignlengths) = push!(alignlengths, len)



"""
    bamrecord_meth_stats(record, genomereader, mtA=0, mtC=0)

    Calculate methylation stats for BAM `record` using `genomereader` to comapre to genome sequence

    stats of interest
    - mapping to chromosomes
    - read length
    - read context of methylation
    - mean meth bases per read
"""

function bamrecord_meth_stats(record, genomereader, mtA=0, mtC=0)
    aln = BAM.alignment(record)
    seq = BAM.sequence(record);
    auxdict = BAM.auxdata(record)
    ml = auxdict["Ml"]::Vector{UInt8}

    ### get ref sequence spanning the read
    refseq = FASTA.extract(genomereader, BAM.refname(record), leftposition(record):rightposition(record));

    posmatch = seq2ref.(Ref(aln), 1:length(seq))
    pos = first.(posmatch)
    ops = last.(posmatch)
    relpos = pos .- leftposition(record) .+ 1

    ## match base pairs to the genome
    @assert all(ops .∈ Ref(Set([OP_SOFT_CLIP, OP_MATCH, OP_INSERT])))
    match_ind = ops .== Ref(OP_MATCH)
    seq_match = falses(length(match_ind))
    seq_match[match_ind] = LongDNA{4}(refseq[relpos[match_ind]]) .== seq[match_ind]

    ### get methlation calls
    mmA = fill(-1, length(seq))
    mmC = fill(-1, length(seq))
    atind = ifelse(BAM.ispositivestrand(record), seq .== DNA_A, seq .== DNA_T)
    gcind = ifelse(BAM.ispositivestrand(record), seq .== DNA_C, seq .== DNA_G)


    mmA[atind] .= ml[1:sum(atind)]
    mmC[gcind] .= ml[sum(atind) .+ (1:sum(gcind))]
    
    valid_mmA = ifelse.(seq_match, mmA, -1)
    valid_mmC = ifelse.(seq_match, mmC, -1)
    
    ### methylation run length stats, run lengths identify consequetive methylated basepairs
    meth_symb_A, run_A = rle(filter(m -> m != -1, valid_mmA) .> mtA)
    meth_symb_C, run_C = rle(filter(m -> m != -1, valid_mmC) .> mtC)
    
    meth_run_stats_A = methrunstats(meth_symb_A, run_A)
    meth_run_stats_C = methrunstats(meth_symb_C, run_C)
 
    
    total_A = sum(mmA .!= -1)
    total_C = sum(mmC .!= -1)
    valid_A = sum(valid_mmA .!= -1)
    valid_C = sum(valid_mmC .!= -1)
    
    total_mA = sum(mmA .> mtA)
    total_mC = sum(mmC .> mtC)
    valid_mA = sum(valid_mmA .> mtA)
    valid_mC = sum(valid_mmC .> mtC)
    
    mA_rate  = total_mA./total_A
    mC_rate  = total_mC./total_C
    
    valid_mA_rate = valid_mA./valid_A
    valid_mC_rate = valid_mC./valid_C

    data = (seqname=BAM.seqname(record), total_matches=sum(seq_match), 
    total_A  = total_A ,
    total_C  = total_C ,
    valid_A  = valid_A ,
    valid_C  = valid_C ,
    total_mA = total_mA,
    total_mC = total_mC,
    valid_mA = valid_mA,
    valid_mC = valid_mC,
    mA_rate  = mA_rate ,
    mC_rate  = mC_rate ,

    valid_mA_rate = valid_mA_rate,
    valid_mC_rate = valid_mC_rate, 
    total_meth_runs_A = meth_run_stats_A.total_meth_runs,
    meth_rl_mu_A = meth_run_stats_A.meth_rl_mu,
    meth_rl_std_A = meth_run_stats_A.meth_rl_std,
    nmeth_rl_mu_A = meth_run_stats_A.nmeth_rl_mu,
    nmeth_rl_std_A = meth_run_stats_A.nmeth_rl_std,
    meth_min_run_A = meth_run_stats_A.meth_min_run,
    meth_max_run_A = meth_run_stats_A.meth_max_run,
    nmeth_min_run_A = meth_run_stats_A.nmeth_min_run,
    nmeth_max_run_A = meth_run_stats_A.nmeth_max_run,
    
    total_meth_runs_C = meth_run_stats_C.total_meth_runs,
    meth_rl_mu_C = meth_run_stats_C.meth_rl_mu,
    meth_rl_std_C = meth_run_stats_C.meth_rl_std,
    nmeth_rl_mu_C = meth_run_stats_C.nmeth_rl_mu,
    nmeth_rl_std_C = meth_run_stats_C.nmeth_rl_std,
    meth_min_run_C = meth_run_stats_C.meth_min_run,
    meth_max_run_C = meth_run_stats_C.meth_max_run,
    nmeth_min_run_C = meth_run_stats_C.nmeth_min_run,
    nmeth_max_run_C = meth_run_stats_C.nmeth_max_run)

    data
end

"""
    methrunstats(symb, run)

    Calculate methylation run length statistics
"""
function methrunstats(symb, run)
    total_meth_runs = sum(symb .== 1)
    total_nmeth_runs = sum(symb .== 0)
    
    if total_meth_runs > 0
        meth_rl_mu, meth_rl_std = mean_and_std(run[symb .== 1])
        meth_min_run, meth_max_run = extrema(run[symb .== 1])
        
    else
        meth_rl_mu, meth_rl_std = 0.0, 0.0
        meth_min_run, meth_max_run = 0, 0
    end
    
    if total_nmeth_runs > 0
        nmeth_rl_mu, nmeth_rl_std = mean_and_std(run[symb .== 0])
        nmeth_min_run, nmeth_max_run = extrema(run[symb .== 0])
    else
        nmeth_rl_mu, nmeth_rl_std = 0.0, 0.0 
        nmeth_min_run, nmeth_max_run = 0, 0
    end
    
    (; total_meth_runs, meth_rl_mu, meth_rl_std, nmeth_rl_mu, nmeth_rl_std, meth_min_run, meth_max_run, nmeth_min_run, nmeth_max_run)
end


function bam_meth_stats(bamfile, genomefile, outfile=nothing)

    genomereader = open(FASTA.Reader, genomefile, index=string(genomefile, ".fai"))

    reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
    total_mapped_reads = mapreduce(d -> d[end].n_mapped, +, reader.index.index.data)
    p = Progress(total_mapped_reads)

    meth_stats = [(next!(p); bamrecord_meth_stats(r, genomereader)) for r in reader if validfrag(r)]

    update!(p, total_mapped_reads)
    close(reader)
    df = DataFrame(meth_stats)

    if !isnothing(outfile)
        CSV.write(outfile, df, delim='\t', compress=true)
    end

    df
end
