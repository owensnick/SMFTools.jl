struct PosMatrix{T, V}

    name::String
    modification::String

    nummods::Int
    modfields::Int

    numreads::Int
    readfields::Int

    chromindex::Dict{String, UnitRange{Int}}
    chromrindex::Dict{Int, String}

    readindex::V
    readivs::IntervalCollection{UnitRange{Int}}

    PM::T

end

struct PosMatrixCompressed{T, V, W}

    name::String
    modifications::Dict{String, Int}

    nummods::Int
    modfields::Int

    numreads::Int
    readfields::Int

    chromindex::Dict{String, UnitRange{Int}}
    chromrindex::Dict{Int, String}

    readindex::V
    readivs::IntervalCollection{W}

    dataencoding::DataType
    PM::T

end

function Base.show(io::IO, pm::PosMatrix)
    if get(io, :compact, false)
        ### compact
        print(io, "Pos Matrix Indexed: ", pm.name)
        print(io, "    ", pm.nummods)
        print(io, "    ", pm.numreads)
    else
        println(io, "Index PM      :   ", pm.name)
        println(io, "nummods       :   ", pm.nummods)
        println(io, "numreads      :   ", pm.numreads)
        println(io, "numfields     :   ", pm.modfields)
        # println(io, "fields        :   ", fm.fields)
        # println(io, "totalecfrags  :   ", fm.totalecfrags)
        # println(io, "index         :   ", fm.index)
        # println(io, "chromindex    :   ", fm.chromindex)
        println(io, "Read  Index   :   ", summary(pm.readindex))
        println(io, "PM            :   ", summary(pm.PM))
    end
end


function load_pos_matrix(file)

    io = open(file)

    modifications, nummods, numreads, numpmfields, numrifields, pmfields, rifields, read_table_pos, index, chromindex, dataencoding, compressed = pm_read_header_footer(io)


    if compressed
        a = Mmap.mmap(io, Vector{UInt8}, sizeof(UInt8)*nummods, position(io))
        P = UnalignedVector{UInt8}(a)
    else
        a = Mmap.mmap(io, Vector{UInt8}, sizeof(dataencoding)*numpmfields*nummods, position(io))
        ua = UnalignedVector{dataencoding}(a)
        P = reshape(ua, (numpmfields, nummods))
    end

    i = Mmap.mmap(io, Vector{UInt8}, sizeof(Int)*numrifields*numreads, read_table_pos)
    iua = UnalignedVector{Int64}(i)
    RI = reshape(iua, (numrifields, numreads))
    close(io)

    ivs = Vector{Interval{Int}}(undef, size(RI, 2))

    for i = 1:size(RI, 2)
        strand, ci = GenomeFragments.get_strand_chrom_enc(RI[1, i])
        ivs[i] = Interval(chromindex[ci], RI[2, i], RI[3, i], strand, i) #RI[4, i]:RI[5, i])
    end
    sort!(ivs)

    ic = IntervalCollection(ivs)

    if compressed
        moddict = Dict(m => i for (i, m) in enumerate(modifications))
        PosMatrixCompressed{typeof(P), typeof(RI), Int}(file, moddict, nummods, numpmfields, numreads, numrifields, index, chromindex, RI, ic, dataencoding, P)
    else
        PosMatrix{typeof(P), typeof(RI)}(file, modification, nummods, numpmfields, numreads, numrifields, index, chromindex, RI, ic, P)
    end

end



function pm_read_header_footer(io,)

    fileversion = read(io, Int64)
    @assert fileversion == 0


    modifications = String[]
    nummods = -1
    numreads = -1

    numpmfields = -1
    numrifields = -1
    pmfields = String[]
    rifields = String[]

    read_table_pos = 0 

    index = Dict{String, UnitRange{Int}}()
    chromindex = Dict{Int,  String}()
    compressed = false

    dataencoding = Int32

    headerpos = read(io, Int64)
    seek(io, headerpos)


    while true
        fields = split(readline(io), '\t')
       (fields[1] == "#modifications")     && (modifications     = split(fields[2], ","))
       (fields[1] == "#compressed")     && (compressed = true)
       (fields[1] == "#totalmodifications")   && (nummods   = parse(Int,  fields[2]))
       (fields[1] == "#totalreads")   && (numreads   = parse(Int,  fields[2]))
       (fields[1] == "#posmatrixfields")   && (numpmfields   = parse(Int,  fields[2]))
       (fields[1] == "#posmatrixlabels")   && (pmfields   = split(fields[2], r"[, ]", keepempty=false))
       (fields[1] == "#readindexfields")   && (numrifields   = parse(Int,  fields[2]))
       (fields[1] == "#readindexlabels")   && (rifields   = split(fields[2], r"[, ]", keepempty=false))
       (fields[1] == "#readindexposition")   && (read_table_pos   = parse(Int,  fields[2]))

       (fields[1] == "#posmatrixencoding")   && (dataencoding   = GenomeFragments.get_datatype(fields[2]))

   





       if fields[1] == "#index"
           index_str = split(replace(chomp(fields[2]), r"\(|\"" => ""), r"\)[,]*", keepempty=false)
           index_vec = GenomeFragments.parse_chrom_range.(index_str)

           index = Dict(index_vec)
           chromindex = Dict(i => c[1] for (i, c) in enumerate(index_vec))
           break
       end
   end

   seekstart(io)
   read(io, Int64)
   read(io, Int64)


   modifications, nummods, numreads, numpmfields, numrifields, pmfields, rifields, read_table_pos, index, chromindex, dataencoding, compressed
end


function getreadindex(chrom, location, pm::PosMatrix)
    !haskey(pm.chromindex, chrom) && return 1:0
    
    ind = pm.chromindex[chrom]
    starts = view(pm.readindex, 2, ind)
    stops  = view(pm.readindex, 3, ind)

    i = searchsortedfirst(starts, location.start)
    @show i
    while (i > 1) && !isempty(intersect(starts[i-1]:stops[i-1], location))
        i -=1
        @show i
    end
    
    (i > length(starts)) && return 1:0

    j = searchsortedlast(stops, location.stop)
    @show j
    while (j < length(stops)) &&  (starts[j + 1] <= location.stop <= stops[j+1])  
        j += 1
    end
    
    return i:j
end

function getposview(chrom, location, pm::PosMatrix)

    ind = getreadindex(chrom, location, pm)
    @show ind
    if isempty(ind)
        return view(pm.PM, :, 1:0)
    end
    pind = pm.readindex[4, first(ind)]:pm.readindex[5, last(ind)]
    @show pind
    view(pm.PM, :, pind)
end