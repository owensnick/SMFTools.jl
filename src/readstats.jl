

"""
    bamrecord_meth_stats(record, genomereader, mtA=0, mtC=0)

    Calculate methylation stats for BAM `record` using `genomereader` to comapre to genome sequence
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
