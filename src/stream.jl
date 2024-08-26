
function stream_nanopolish_meth(tsvfile, outfile; T=Int32, verbose=true, compress=1)
    verbose && println("[SMF]\tStreaming nanopolish file: $tsvfile to $outfile")

    starttime = time()
    indexfile = string(outfile, ".index")
    io = open(outfile, "w")
    rio = open(indexfile, "w")
    println(io, "#bottom_header")

    pos_buff = 0
    write(io, pos_buff)

    iio = open(tsvfile) |> GzipDecompressorStream
    

    sr = streamfrags_io_compress(iio, io, rio, T=T, filelength=filesize(tsvfile), level=compress)
    verbose && println("[SMF]\tStreamed ", sum(sr.totalreads), " fragments from ", length(sr.chroms), " chromosomes")

    close(rio)

    # verbose && println("[SMF]\tWriting header.")
    read_table_pos = position(io)
    
    rio = open(indexfile)
    @time "[SMF]\tCopying header" write(io, read(rio))
    close(rio)
    rm(indexfile)

    verbose && println("[SMF]\tWriting header.")

    header_pos = position(io)

    # header 

    # println(io, "#methylation\t", modification)
    if compress > 0
        println(io, "#compressed\t", compress)
    end
    @show sr.modifications
    println(io, "#modifications\t",   join(sr.modifications, ","))
    println(io, "#totalmodifications\t",   sum(sr.totalpos))
    println(io, "#totalreads\t", sum(sr.totalreads))
    println(io, "#posmatrixfields\t", 2)
    println(io, "#posmatrixlabels\t", "postion, methprob")
    println(io, "#posmatrixencoding\t", T)
    println(io, "#readindexfields\t", 3 + length(sr.modifications)*2)
    println(io, "#readindexlabels\t", "strand_chrom_enc, leftposition, rightposition, $(length(sr.modifications))x(modstartindex, modstopindex)")
    println(io, "#readindexposition\t", read_table_pos)
    println(io, "#readindexencoding\t", Int)

    
    cumreads = cumsum(sr.totalreads)
    readind  = map((s, e) -> s:e, [0 ; cumreads[1:end-1]] .+ 1, cumreads)
    println(io, "#index\t", join(zip(sr.chroms, readind), ","))
    seekstart(io)
    # println(io, "#bottom_header")
    write(io, 0)
    write(io, header_pos)
    # write(io, pos_buff)

    close(io)
    verbose && println("[SMF]\tcomplete in ", time() - starttime, " seconds")

    ############################################

    # pos_buff = position(io)

    
    ### fields
    # start
    # 2. stop
    # 3. strand_chrom_enc
    # 4. read_index
    # 5. position
    # 6. meth prob
    # labels = ["read_index", "position", "meth_prob"]
    # println(io, "#paired\t",   false)
    # println(io, "#numregions\t",   length(chroms))
    # println(io, "#totalfrags\t", sum(totalpos))
    # println(io, "#totalecfrags\t", [sum(totalreads)])
    # println(io, "#numfields\t", 6)
    # println(io, "#fields\t[", join(["start" ; "stop" ; "strand_chrom_enc"; labels], ", "), "]")
    # println(io, "#dataencoding\t", T)
    # frags    = totalpos
    # cumfrags = cumsum(frags)
    # fragind  = map((s, e) -> s:e, [0 ; cumfrags[1:end-1]] .+ 1, cumfrags)
    # println(io, "#index\t", join(zip(chroms, fragind), ","))
    # seekstart(io)
    # println(io, "#bottom_header")
    # write(io, pos_buff)

    # close(io)
    # verbose && println("[SMF]\tcomplete in ", time() - starttime, " seconds")

end




function streambam_smf(bamfile, outfile; T=Int32, verbose=true)
    verbose && println("[SMF]\tStreaming BAM file: $bamfile to $outfile")

    starttime = time()

    io = open(outfile, "w")
    println(io, "#bottom_header")

    pos_buff = 0
    write(io, pos_buff)
    

    if isfile(string(bamfile, ".bai"))
        reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
        totalreads = totalreadsindex(reader)
        verbose && println("[SMF]\tTotal reads in BAM file: ", totalreads)
    else
        reader = open(BAM.Reader, bamfile)
        totalreads = -1
    end
    chroms, totalreads, totalpos = streamfrags(reader, io, T=T)
    verbose && println("[SMF]\tStreamed ", sum(totalreads), " fragments from ", length(chroms), " chromosomes")
    
    
    #index_table = (chrom=chroms, reads=totalreads, pos=totalpos)
    
    verbose && println("[SMF]\tWriting header.")

    pos_buff = position(io)

    
    ### fields
    # start
    # 2. stop
    # 3. strand_chrom_enc
    # 4. read_index
    # 5. position
    # 6. meth prob
    labels = ["read_index", "position", "meth_prob"]
    println(io, "#paired\t",   false)
    println(io, "#numregions\t",   length(chroms))
    println(io, "#totalfrags\t", sum(totalpos))
    println(io, "#totalecfrags\t", [sum(totalreads)])
    println(io, "#numfields\t", 6)
    println(io, "#fields\t[", join(["start" ; "stop" ; "strand_chrom_enc"; labels], ", "), "]")
    println(io, "#dataencoding\t", T)
    frags    = totalpos
    cumfrags = cumsum(frags)
    fragind  = map((s, e) -> s:e, [0 ; cumfrags[1:end-1]] .+ 1, cumfrags)
    println(io, "#index\t", join(zip(chroms, fragind), ","))
    seekstart(io)
    println(io, "#bottom_header")
    write(io, pos_buff)

    close(io)
    verbose && println("[SMF]\tcomplete in ", time() - starttime, " seconds")

end


"""
    validfrag(record)

    Function to filter BAM `record` to ensure it has an `Ml` field and is not a secondary alignment.
"""
function validfrag(record, mlt="ML", mlt_vold = "Ml")
    f = BAM.flags(record)
    
    if BAM.ismapped(record) && !iszero(BAM.refid(record)) && (haskey(record, mlt) || haskey(record, mlt_vold)) && ((f == 0) || (f == 16))
        return true
    else
        return false
    end
    
end

"""
    totalreadsindex(reader)

    Function to get the total number of reads in a BAM file reader `reader` using the index
"""
function totalreadsindex(reader)
    if isnothing(reader.index)
        return -1
    else
        return mapreduce(d -> isnothing(d[end]) ? 0 : d[end].n_mapped, +, reader.index.index.data) + reader.index.n_no_coors
    end

end


"""
    streamfrags(bamreader, io; T=Int32, filtfun=validfrag, mlt = 0.0)
    
    Function to stream a BAM file reader `bamreader` and write to io stream `io`.

    BAM entries are filtered by `filtfun` and methylation probabilities are filtered by `mlt`.

    Assumes a Nanopore BAM file generated by guppy with methylation probabilities to a binary file that can be loaded by memory mapping with `T` specifying the data type to store information.

"""
function streamfrags(bamreader, io; T=Int32, filtfun=validfrag, mlt = 0.0)

    ### get chromosome names from the bam header
    chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    chromindex = Dict(chroms[i] => i for i = 1:length(chroms))
    

    
    totalreads = zeros(Int, length(chroms))
    totalpos = zeros(Int, length(chroms))
    
    index = 0
    readid = 1

    total_reads_in_file = totalreadsindex(bamreader)
    if total_reads_in_file != -1
        p = Progress(total_reads_in_file, 1, "Streaming BAM file: ", 50)
    else
        p = nothing
    end

    for record in bamreader
        !isnothing(p) && next!(p)
        !filtfun(record) && continue
        
        try
            chrom = BAM.refname(record)
        catch err
            display(record)
            display(BAM.ismapped(record))
            rethrow(err)
        end
        chrom = BAM.refname(record)
        fragstart = BAM.leftposition(record)
        fragstop = BAM.rightposition(record)
        fragstrand = GenomeFragments.getstrand(record)
        chromind = chromindex[chrom]
        
        ### get methylation probabilities
        pos, mlp = posmlp(record)

        haswritten = false
        for (p, ml) in zip(pos, mlp)
            ### reads with zero mehtylated bases will currently not appear in the output
            if ml > mlt
                write(io, T(fragstart))         ### start
                write(io, T(fragstop))          ### stop
                write(io, GenomeFragments.set_strand_chrom_enc(fragstrand, chromind, T))          ### chrom strand enc
                write(io, T(readid)) ## read index
                write(io, T(p)) ## position
                write(io, T(ml)) ## methylation probability
                totalpos[chromind] += 1
                haswritten=true
            end 
        end
        if haswritten
            totalreads[chromind] += 1
        end
        readid += 1
    end
    chroms, totalreads, totalpos
end


function streambam_smf_pos(bamfile, outfile, modifiction="";  T=Int32, verbose=true, mlt=0, compress=1)
    verbose && println("[SMF]\tStreaming BAM file: $bamfile to $outfile")

    starttime = time()
    indexfile = string(outfile, ".index")
    io = open(outfile, "w") # main file
    iio = open(indexfile, "w")


    write(io, 0) ### organisation involving FM matrix, readtable, chromindex
    # println(io, "#bottom_header")

    header_pos = 0
    write(io, header_pos)

    read_table_pos = 0
    # total_reads = 0
    # total_read_fields = 0
    # write(io, read_table_pos)
    # write(io, total_reads)
    # write(io, total_read_fields)
    
    # chrom_index_pos = 0
    # write(io, chrom_index_pos)
    

    if isfile(string(bamfile, ".bai"))
        reader = open(BAM.Reader, bamfile, index=string(bamfile, ".bai"))
        totalreads = totalreadsindex(reader)
        verbose && println("[SMF]\tTotal reads in BAM file: ", totalreads)
    else
        reader = open(BAM.Reader, bamfile)
        totalreads = -1
    end
    # chroms, totalreads, totalpos 
    if compress > 0
        sr = streamposfrags_compress(reader, io, iio; T=T, mlt = mlt, level=compress)
    else
        sr = streamposfrags(reader, io, iio, T=T, mlt=mlt)
    end
    verbose && println("[SMF]\tStreamed ", sum(sr.totalreads), " fragments from ", length(sr.chroms), " chromosomes")
    close(iio)



    verbose && println("[SMF]\tWriting header.")
    read_table_pos = position(io)
    
    rio = open(indexfile)
    @time "[SMF]\tCopying header" write(io, read(rio))
    close(rio)
    rm(indexfile)
    header_pos = position(io)

    # header 

    # println(io, "#methylation\t", modification)
    if compress > 0
        println(io, "#compressed\t", compress)
    end
    println(io, "#modifications\t",   join(sr.modifications, ","))
    println(io, "#totalmodifications\t",   sum(sr.totalpos))
    println(io, "#totalreads\t", sum(sr.totalreads))
    println(io, "#posmatrixfields\t", 2)
    println(io, "#posmatrixlabels\t", "postion, methprob")
    println(io, "#posmatrixencoding\t", T)
    println(io, "#readindexfields\t", 3 + length(sr.modifications)*2)
    println(io, "#readindexlabels\t", "strand_chrom_enc, leftposition, rightposition, $(length(sr.modifications))x(modstartindex, modstopindex)")
    println(io, "#readindexposition\t", read_table_pos)
    println(io, "#readindexencoding\t", Int)

    
    cumreads = cumsum(sr.totalreads)
    readind  = map((s, e) -> s:e, [0 ; cumreads[1:end-1]] .+ 1, cumreads)
    println(io, "#index\t", join(zip(sr.chroms, readind), ","))
    seekstart(io)
    # println(io, "#bottom_header")
    write(io, 0)
    write(io, header_pos)
    # write(io, pos_buff)

    close(io)
    verbose && println("[SMF]\tcomplete in ", time() - starttime, " seconds")


    
    ### fields for Pos matrix
    # 1. position
    # 2. meth_prob
    

    ### fields for read index
    # 1. strand_chrom_enc
    # 2. leftposition
    # 3. rightposition
    # 4. modindex
    # 5. num mods



    # labels = ["read_index", "position", "meth_prob"]
    # println(io, "#paired\t",   false)
    # println(io, "#numregions\t",   length(chroms))
    # println(io, "#totalfrags\t", sum(totalpos))
    # println(io, "#totalecfrags\t", [sum(totalreads)])
    # println(io, "#numfields\t", 6)
    # println(io, "#fields\t[", join(["start" ; "stop" ; "strand_chrom_enc"; labels], ", "), "]")
    # println(io, "#dataencoding\t", T)
    # frags    = totalpos
    # cumfrags = cumsum(frags)
    # fragind  = map((s, e) -> s:e, [0 ; cumfrags[1:end-1]] .+ 1, cumfrags)
    # println(io, "#index\t", join(zip(chroms, fragind), ","))
    # seekstart(io)
    # println(io, "#bottom_header")
    # write(io, pos_buff)

    # close(io)
    # verbose && println("[SMF]\tcomplete in ", time() - starttime, " seconds")

end

function streamposfrags(bamreader, io, iio; T=Int32, filtfun=validfrag, mlt = 0.0)

    ### get chromosome names from the bam header
    chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    chromindex = Dict(chroms[i] => i for i = 1:length(chroms))
    

    
    totalreads = zeros(Int, length(chroms))
    totalpos = zeros(Int, length(chroms))

    ### read index
    # 1. strand_chrom_enc
    # 2. leftposition
    # 3. rightposition
    # 4. mod_start_index
    # 5. mod_stop_index
    
    # strand_chrom_enc = Int[]
    # leftposition = Int[]
    # rightposition = Int[]
    # modstartindex = Int[]
    # nummods = Int[]
    
    
    
    readid = 1
    modid = 1

    total_reads_in_file = totalreadsindex(bamreader)
    if total_reads_in_file != -1
        p = Progress(total_reads_in_file, 1, "Streaming BAM file: ", 50)
    else
        p = nothing
    end

    for record in bamreader
        !isnothing(p) && next!(p)
        !filtfun(record) && continue
        
        # try
        #     chrom = BAM.refname(record)
        # catch err
        #     display(record)
        #     display(BAM.ismapped(record))
        #     rethrow(err)
        # end
        chrom = BAM.refname(record)
        fragstart = BAM.leftposition(record)
        fragstop = BAM.rightposition(record)
        fragstrand = GenomeFragments.getstrand(record)
        chromind = chromindex[chrom]
        
        ### get methylation probabilities
        pos, mlp = posmlp(record)


        write(iio, GenomeFragments.set_strand_chrom_enc(fragstrand, chromind, Int))          ### chrom strand enc
        write(iio, fragstart)
        write(iio, fragstop)
        write(iio, modid)

        
        totalreads[chromind] += 1
        numwritten = 0
        for (p, ml) in zip(pos, mlp)
            ### reads with zero mehtylated bases will currently not appear in the output
            if ml > mlt
                # write(io, T(fragstart))         ### start
                # write(io, T(fragstop))          ### stop
                # write(io, GenomeFragments.set_strand_chrom_enc(fragstrand, chromind, T))          ### chrom strand enc
                # write(io, T(readid)) ## read index
                write(io, T(p)) ## position
                write(io, T(ml)) ## methylation probability
                numwritten += 1
                modid += 1
                totalpos[chromind] += 1
                
            end 
        end
        write(iio, modid - 1) ### TODO currently this gives redundant information with the subsequent start, we can remove this
        
        # if numwritten > 0
           
        # end
        readid += 1
    end
    (; chroms, totalreads, totalpos)
end



function streamposfrags_compress(bamreader, io, iio; T=Int32, filtfun=validfrag, mlt = 0.0, level=1)

    ### get chromosome names from the bam header
    chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    chromindex = Dict(chroms[i] => i for i = 1:length(chroms))
    

    
    totalreads = zeros(Int, length(chroms))
    totalpos = zeros(Int, length(chroms))

    ### read index
    # 1. strand_chrom_enc
    # 2. leftposition
    # 3. rightposition
    # 4. mod_start_index
    # 5. mod_stop_index
    
    # strand_chrom_enc = Int[]
    # leftposition = Int[]
    # rightposition = Int[]
    # modstartindex = Int[]
    # nummods = Int[]
    
    
    
    readid = 1
    modid = 1

    total_reads_in_file = totalreadsindex(bamreader)
    if total_reads_in_file != -1
        p = Progress(total_reads_in_file, 1, "Streaming BAM file: ", 50)
    else
        p = nothing
    end

    modifications = String[]

    for record in bamreader
        !isnothing(p) && next!(p)
        !filtfun(record) && continue
        
        # try
        #     chrom = BAM.refname(record)
        # catch err
        #     display(record)
        #     display(BAM.ismapped(record))
        #     rethrow(err)
        # end
        chrom = BAM.refname(record)
        fragstart = BAM.leftposition(record)
        fragstop = BAM.rightposition(record)
        fragstrand = GenomeFragments.getstrand(record)
        chromind = chromindex[chrom]
        
        ### get methylation probabilities
        mods = methcalls_cg_gc(record)
        if isempty(modifications)
            modifications = mods.mods
        # elseif length(modifications) < length(mods.mods)
        #     for m in modifications
        #         @assert m ∈ mods.mods
        #     end
        # else
        #     for m in mods.mods
        #         @assert m ∈ modifications
        #     end
            try
                @assert modifications == mods.mods
            catch
                display(modifications)
                display(mods.mods)
                error("")
            end
        end


        # pos, mlp = posmlp(record)
        totalreads[chromind] += 1

        write(iio, GenomeFragments.set_strand_chrom_enc(fragstrand, chromind, Int))          ### chrom strand enc
        write(iio, fragstart)
        write(iio, fragstop)

        for mi = 1:length(modifications)
            write(iio, modid)
            ind = mods.mls[mi] .> mlt
            M = T.([mods.pos[mi][ind]' ; mods.mls[mi][ind]'])
            C = compress(M, level=level)
            w = write(io, C)
            modid += w
            totalpos[chromind] += w
            write(iio, modid - 1) ### TODO currently this gives redundant information with the subsequent start, we can remove this
        end

        readid += 1

        # write(iio, modid)
        # totalreads[chromind] += 1
        # # numwritten = 0
        # ind = mlp .> mlt
        # M = T.([pos[ind]' ; mlp[ind]'])
        # C = compress(M, level=level)
        # w = write(io, C)
        # modid += w
        # totalpos[chromind] += w
        # write(iio, modid - 1) ### TODO currently this gives redundant information with the subsequent start, we can remove this

        # readid += 1
    end
    modifications

    (; chroms, totalreads, totalpos, modifications)
end

gc_lr_int_score(lr; scale = 100_000, T=Int32) = min(Int(round(lr*scale)), typemax(T))
function parseline(line)
   
    chrom = SubString{String}("")
    strand = SubString{String}("")
    start = 0
    stop = 0
    readid = SubString{String}("")
    lr = 0
    nummot = 0
    seq = SubString{String}("")
    

    for (i, s) in enumerate(eachsplit(line))
       if i == 1
            chrom = s
        elseif i == 2
            strand = s
        elseif i == 3
            start = Parsers.parse(Int, s)
        elseif i == 4
            stop = Parsers.parse(Int, s)
        elseif i == 5
            readid = s
        elseif i == 6
            lr = Parsers.parse(Float64, s) |> gc_lr_int_score
        elseif i == 10
            nummot = Parsers.parse(Int, s)
        elseif i == 11
            seq = s
        end
    end
    
    (; chrom, strand, start, stop, readid, lr, nummot, seq)
end

valid_gc_call(seq) = !occursin(r"CGC|GCG", seq)

function addgc!(pos, score, start, seq, numgc, lr)
    currentmatch = 6
    push!(pos, start + currentmatch - 5)
    push!(score, lr)
    for i = 2:numgc
        ind = findnext("GC", seq, currentmatch + 2)
        currentmatch = first(ind)
        push!(pos, start + currentmatch - 5)
        push!(score, lr)
    end
end

"""
    streamfrags_io(file, io; T=Int32, filtfun=validfrag, mlt = 0.0)
    
    Function to stream a a methylation TSV file reader `iio` and write to io stream `io`.

    Assumes nanopolish methylation TSV file format.

"""
function streamfrags_io(file, io; T=Int32, filtfun=validfrag, mlt = 0.0)

    ### get chromosome names from the bam header

    # iio = open(tsvfile) |> GzipDecompressorStream
    filelength = filesize(file)

    chroms = String[]
    chromindex = Dict{String, Int}()

    # chroms = [v["SN"] for v in findall(BAM.header(bamreader), "SQ")]
    # chromindex = Dict(chroms[i] => i for i = 1:length(chroms))

    totalreads = Dict{String, Int}()
    totalpos = Dict{String, Int}()

   
    
    p = Progress(filelength, 1, "Streaming Meth file: ", 50)

    currentchrom = ""
    ci = 1
    
    readindex = 0
    readid = ""
    readstrand = "+"
    readchrom = ""
    readpos = Int[]
    readscore = Int[]


    for line in eachline(iio)
        occursin(r"^chromosome", line) && continue
        res = parseline(line)

        if res.chrom != currentchrom
            push!(chroms, res.chrom)
            chromindex[res.chrom] = ci
            currentchrom = res.chrom
            ci += 1
            totalreads[res.chrom] = 0
            totalpos[res.chrom] = 0
        end
        !valid_gc_call(res.seq) && continue

        if readid != res.readid
            if readid != ""
                readstart = first(readpos)
                readstop = last(readpos)
                for (pos, score) in zip(readpos, readscore)
                    write(io, T(readstart))         ### start
                    write(io, T(readstop))          ### stop
                    write(io, GenomeFragments.set_strand_chrom_enc(readstrand, chromindex[readchrom], T))          ### chrom strand enc
                    write(io, T(readindex)) ## read index
                    write(io, T(pos)) ## position
                    write(io, T(score)) ## methylation probability
                    totalpos[readchrom] += 1
                end
                totalreads[readchrom] += 1
                update!(p, position(iio))
            end
            readid = res.readid
            readstrand = res.strand
            readchrom = res.chrom
            readpos = Int[]
            readscore = Int[]
            readindex += 1
        end
        addgc!(readpos, readscore, res.start, res.seq, res.nummot, res.lr)

    end

    if readid != ""
        readstart = first(readpos)
        readstop = last(readpos)
        for (pos, score) in zip(readpos, readscore)
            write(io, T(readstart))         ### start
            write(io, T(readstop))          ### stop
            write(io, GenomeFragments.set_strand_chrom_enc(readstrand, chromindex[readchrom], T))          ### chrom strand enc
            write(io, T(readindex)) ## read index
            write(io, T(pos)) ## position
            write(io, T(score)) ## methylation probability
            totalpos[readchrom] += 1
        end
        totalreads[readchrom] += 1
    end

    update!(p, filelength)



    chroms, [totalreads[c] for c in chroms], [totalpos[c] for c in chroms]
end


"""
    streamfrags_io(file, io; T=Int32, filtfun=validfrag, mlt = 0.0)
    
    Function to stream a a methylation TSV file reader `iio` and write to io stream `io`.

    Assumes nanopolish methylation TSV file format.

"""
function streamfrags_io_compress(iio, io, rio; T=Int32, filtfun=validfrag, mlt = 0.0, level=1, filelength=0)

    ### get chromosome names from the bam header

     ### read index
    # 1. strand_chrom_enc
    # 2. leftposition
    # 3. rightposition
    # 4. mod_start_index
    # 5. mod_stop_index



    chroms = String[]
    chromindex = Dict{String, Int}()

    totalreads = Dict{String, Int}()
    totalpos = Dict{String, Int}()

    if filelength > 0
        p = Progress(filelength, 1, "Streaming Meth file: ", 50)
    end

    currentchrom = ""
    ci = 1
    
    readindex = 0
    readid = ""
    readstrand = ""
    readchrom = ""
    readpos = Int[]
    readscore = Int[]

    modid = 1
    for line in eachline(iio)
        occursin(r"^chromosome", line) && continue
        res = parseline(line)

        if res.chrom != currentchrom
            push!(chroms, res.chrom)
            chromindex[res.chrom] = ci
            @show res.chrom, ci
            currentchrom = res.chrom
            ci += 1 
            totalreads[res.chrom] = 0
            totalpos[res.chrom] = 0
        end
        !valid_gc_call(res.seq) && continue

        if readid != res.readid
            if readid != ""
                readstart = first(readpos)
                readstop = last(readpos)

           

                write(rio, GenomeFragments.set_strand_chrom_enc(readstrand, chromindex[readchrom], Int))          ### chrom strand enc
                write(rio, readstart)
                write(rio, readstop)
                write(rio, modid)
                # try
                #     @show readscore, readpos
                #     ind = readscore .> mlt
                #     @show ind
                #     M = T.([readpos[ind]' ; readscore[ind]'])
                #     C = compress(M, level=level)
    
                # catch
                #     error("")
                # end
                ind = readscore .> mlt
                M = T.([readpos[ind]' ; readscore[ind]'])
                C = compress(M, level=level)

                w = write(io, C)
                modid += w
                write(rio, modid -1)
                totalpos[readchrom] += w

                totalreads[readchrom] += 1
                (filelength > 0) && update!(p, position(iio))
            end
            readid = res.readid
            readstrand = res.strand
            readchrom = res.chrom
            readpos = Int[]
            readscore = Int[]
            readindex += 1
        end
        addgc!(readpos, readscore, res.start, res.seq, res.nummot, res.lr)

    end

    if readid != ""
        readstart = first(readpos)
        readstop = last(readpos)
       
        write(rio, GenomeFragments.set_strand_chrom_enc(readstrand, chromindex[readchrom], Int))          ### chrom strand enc
        write(rio, readstart)
        write(rio, readstop)
        write(rio, modid)

        ind = readscore .> mlt
        M = T.([readpos[ind]' ; readpos[ind]'])
        C = compress(M, level=level)

        w = write(io, C)
        modid += w
        write(rio, modid -1)
        totalpos[readchrom] += w

        totalreads[readchrom] += 1
    end

    (filelength > 0) && update!(p, filelength)

    (chroms=chroms,  totalreads=[totalreads[c] for c in chroms], totalpos=[totalpos[c] for c in chroms], modifications=["G5mC"])
end

function get_ml_data(record)
    if haskey(record, "Ml")
        return record["Ml"]::Vector{UInt8}
    elseif haskey(record, "ML")
        return record["ML"]::Vector{UInt8}
    else
        error("ML not found:\n$(record)")
    end
end

function get_mm_data(record)
    if haskey(record, "Mm")
        return record["Mm"]::String
    elseif haskey(record, "MM")
        return record["MM"]::String
    else
        error("MM not found:\n$(record)")
    end
end

"""
    posmlp(record)

    Retrieve the methylation probabilities for the BAM `record` from the `Mm` and `Ml` fields

    *Warning: currently only returns 6mA methylation probabilities*
    *Warning: currently returns methylation probabilities w.r.t the read not the genome, ie some positions may not match the genome*

"""
# function posmlp(record; flipstrand=false)
#     pos = BAM.ispositivestrand(record)
#     @show pos
#     # auxdict = BAM.auxdata(record)
#     # ml = auxdict["ML"]::Vector{UInt8}
#     ml = get_ml_data(record)
#     aln = BAM.alignment(record)
#     seq = BAM.sequence(record)
#     mli = 1
    
#     ### output
#     positions = Int[]
#     # ops = Operation[]
#     mlp = Int[]

    
#     aln = BAM.alignment(record)
#     seq = BAM.sequence(record)
#     mli = 1
    
#     # positions = Int[]
#     # ops = Operation[]
#     # mlp = Int[]
  
#     mm = get_mm_data(record)
#     @assert mm[1] == 'A'
#     # @assert first(auxdict["Mm"]) == 'A'
    
#     # if pos
#     #     @assert length(ml) == (sum(s -> s == DNA_A, seq) + sum(s -> s == DNA_C, seq))
#     # else

#     #     @assert length(ml) == (sum(s -> s == DNA_T, seq) + sum(s -> s == DNA_G, seq))
#     # end
    
#     # @show "flip"
#     for i = 1:BAM.seqlength(record)
#         if (pos && (seq[i] == DNA_A)) || (!pos && (seq[i] == DNA_T))
#             # k, op = seq2ref(aln, BAM.seqlength(record) - i + 1)

#             k, op = seq2ref(aln, i)
            
#             if (op == OP_MATCH) || (op == OP_SEQ_MATCH)
#                 push!(positions, k)
#                 push!(mlp, ml[mli])
#             end
#             mli += 1

#         end
#     end

    
#     positions, mlp
# end

function posmlp(record)
    pos = BAM.ispositivestrand(record)
   
    # auxdict = BAM.auxdata(record)
    # ml = auxdict["ML"]::Vector{UInt8}
    ml = get_ml_data(record)
    mm = get_mm_data(record)
    @assert mm[1] == 'A'
    aln = BAM.alignment(record)
    seq = BAM.sequence(record)

    atind_called = BAM.ispositivestrand(record) ? (seq .== DNA_A) : (seq .== DNA_T) 
    gcind_called = BAM.ispositivestrand(record) ? (seq .== DNA_C) : (seq .== DNA_G) 

    totalat = sum(atind_called)

    ml_A = ml[1:totalat]
    ml_C = ml[(totalat + 1):end]

    if !pos
        reverse!(ml_A)
        reverse!(ml_C)
    end
    positions = Int[]
    ml_A_ind = falses(length(ml_A))
    mli_A = 1
    for i = 1:length(seq)
        if atind_called[i]
            k, op = seq2ref(aln, i)
            if op == OP_MATCH
                push!(positions, k)
                ml_A_ind[mli_A] = true
            end
            mli_A += 1
        end
    end
    
    positions, ml_A[ml_A_ind]
end

function methcalls(record)
    
    # meths = String[]
    mm = SMFTools.get_mm_data(record)
    # fields = Set{String}()
    
    
    #### obtain runs of modifications
    # modruns = Dict{String, Vector{Int}}()
    mods = String[]
    runs = Vector{Int}[]
    currentmod = ""
    
    for block in eachsplit(mm, ";", keepempty=false)
        firstfield = true
        for f in eachsplit(block, ",")
            if firstfield
                
                if f[1] == 'A'
                   @assert f[2] == '+'
                   @assert f[3] == 'a'
                   currentmod = "6mA"
                elseif f[1] == 'C'
                    @assert f[2] == '+'
                    if f[3] == 'h'
                        currentmod = "5hmC"
                    elseif f[3] == 'm'
                        currentmod = "5mC"
                    else
                        error("Modification $f not recognised")
                    end
                else
                    error("Modification $f not recognised")
                end
                firstfield = false
                push!(runs, Int[])
                push!(mods, currentmod)
            else
                push!(runs[end], Parsers.parse(Int, f))
            end
        end
    end
    
    
    #### align runs with modification calls
    
    ml = SMFTools.get_ml_data(record)
    base_indexes = BitArray{1}[]
    seq = BAM.sequence(record)
    aln = BAM.alignment(record)
    rcdict = Dict('A' => (DNA_A, DNA_T), 'C' => (DNA_C, DNA_G))
    mli = 1
    
    pos = Vector{Int}[]
    mls = Vector{UInt8}[]
    
    
    fA = Int[]
    for (k, (m, run)) in enumerate(zip(mods, runs))
        bp = rcdict[last(m)]
        # @show m, last(m), bp
        
        ind = BAM.ispositivestrand(record) ? (seq .== bp[1]) : (seq .== bp[2])
        push!(base_indexes, ind)
        f_ind = findall(ind)
        
        ml_B = ml[mli .+ (1:length(run)) .- 1]        
        mli += length(run)
        
        if !BAM.ispositivestrand(record)
            reverse!(ml_B)
        end
        mlib = 1
        fi = 0
        
        push!(pos, Int[])
        push!(mls, UInt8[])
        
        for r in run
            fi += (r + 1)
            p, op = seq2ref(aln, f_ind[fi])
            if op == OP_MATCH
                push!(pos[k], p)
                push!(mls[k], ml_B[mlib])
                push!(fA, mlib)
            end
            mlib += 1
        end
        # @show mli
    end

    (;mods, runs, base_indexes, pos, mls)
end


function methcalls_cg_gc(record)
    
    # meths = String[]
    mm = SMFTools.get_mm_data(record)
    # fields = Set{String}()
    
    
    #### obtain runs of modifications
    # modruns = Dict{String, Vector{Int}}()
    mods = String[]
    runs = Vector{Int}[]
    currentmod = ""
    
    for block in eachsplit(mm, ";", keepempty=false)
        firstfield = true
        readingh = false
        for f in eachsplit(block, ",")
            if firstfield
                
                if f[1] == 'A'
                   @assert f[2] == '+'
                   @assert f[3] == 'a'
                   currentmod = "a"
                elseif f[1] == 'C'
                    @assert f[2] == '+'
                    if f[3] == 'h'
                        currentmod = "h"
                        readingh = true
                    elseif f[3] == 'm'
                        currentmod = "m"
                    else
                        error("Modification $f not recognised")
                    end
                else
                    error("Modification $f not recognised")
                end
                firstfield = false
                if !readingh
                    push!(runs, Int[])
                    push!(mods, currentmod)
                end
            else
                !readingh && push!(runs[end], Parsers.parse(Int, f))
            end
        end
    end
    
    
    #### align runs with modification calls
    
    ml = SMFTools.get_ml_data(record)
    base_indexes = BitArray{1}[]
    seq = BAM.sequence(record)
    aln = BAM.alignment(record)
    rcdict = Dict("a" => (DNA_A, DNA_T), "m" => (DNA_C, DNA_G))
    mli = 1
    
    pos = Vector{Int}[]
    mls = Vector{UInt8}[]
    # kmers = Vector{DNAKmer{3}}[]
    
    # fA = Int[]
    k = 0
    split_m = false
    for (m, run) in zip(mods, runs)
        (m == "h") && continue
        k += 1
        
        bp = rcdict[m]
        # @show m, last(m), bp
        
        ind = BAM.ispositivestrand(record) ? (seq .== bp[1]) : (seq .== bp[2])
        push!(base_indexes, ind)
        f_ind = findall(ind)
        
        ml_B = ml[mli .+ (1:length(run)) .- 1]        
        mli += length(run)
        
        if !BAM.ispositivestrand(record)
            reverse!(ml_B)
        end
        mlib = 1
        fi = 0

        if m == "a"
            push!(pos, Int[])
            push!(mls, UInt8[])
            # push!(kmers, DNAKmer{3}[])
        elseif m == "m"
            push!(pos, Int[])
            push!(mls, UInt8[])
            # push!(kmers, DNAKmer{3}[])
            push!(pos, Int[])
            push!(mls, UInt8[])
            # push!(kmers, DNAKmer{3}[])
        end
        
        for r in run
            fi += (r + 1)
            p, op = seq2ref(aln, f_ind[fi])
            
            
            if op == OP_MATCH
                if m == "a"
                    push!(pos[k], p)
                    push!(mls[k], ml_B[mlib])
                    # kmer = DNAKmer{3}(seq[f_ind[fi] .+ (-1:1)])
                    # push!(kmers[k], kmer)
                elseif m == "m"
                    split_m = true
                    class = classify_m_mod(seq, f_ind[fi])
                    
                    if class == :gc
                        kp = k + 1
                        push!(pos[kp], p)
                        push!(mls[kp], ml_B[mlib])
                        # kmer = DNAKmer{3}(seq[f_ind[fi] .+ (-1:1)])
                        # push!(kmers[kp], kmer)
                    elseif class == :cg
                        kp = k
                        push!(pos[kp], p)
                        push!(mls[kp], ml_B[mlib])
                        # kmer = DNAKmer{3}(seq[f_ind[fi] .+ (-1:1)])
                        # push!(kmers[kp], kmer)
                    end
                    
                end
                # push!(fA, mlib)
            end
            mlib += 1
        end
        if m == "m"
            k += 1
        end
        # @show mli
    end

    if split_m
        @assert mods[end] == "m"
        push!(mods, "i")
    end

    (;mods, runs, base_indexes, pos, mls)
end

function classify_m_mod(seq, i)
    if (i == 1) || (i == length(seq))
        return :nomod
    elseif seq[i-1] == DNA_G

        if seq[i+1] == DNA_G
            return :ambig
        else
            return :gc
        end
    elseif seq[i+1] == DNA_G
        return :cg
    else
        return :nomod
    end
end