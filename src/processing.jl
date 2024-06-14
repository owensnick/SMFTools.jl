

function modreadsintersecting(chroms, locations, PM)
    ic = IntervalCollection(Interval.(chroms, locations, Ref(STRAND_NA), 1:length(chroms)))
    modreadsintersecting(ic, PM)
end

function modreadsintersecting(ic::IntervalCollection{T}, pm) where{T}
    total = zeros(Int, length(ic))

    for (ii, ri) in eachoverlap(ic, pm.readivs)
        k = GenomicFeatures.metadata(ii)
        total[k] += 1
    end
    total
end

# function modmeta(chroms, locations, strands, PM; mlt=-1, fn=identity, spos=true, sneg=true)

#     n = length(locations[1])
#     pile = zeros(n)
#     bgpile = zeros(n)

#     for (k, (c, l, s)) in enumerate(zip(chroms, locations, strands))
#         rind = getreadindex(c, l, PM)
        
#         for ri in rind
#             rstrand, _ = GenomeFragments.get_strand_chrom_enc(PM.readindex[1, ri])
#             !spos && (rstrand == STRAND_POS) && continue
#             !sneg && (rstrand == STRAND_NEG) && continue
#             pind = PM.readindex[4, ri]:PM.readindex[5, ri]

#             for pi in pind
#                 pos = PM.PM[1, pi]
#                 ml  = PM.PM[2, pi]
                
#                 if s == "+"
#                     wpos = pos - l.start .+ 1
#                 else
#                     wpos = l.stop - pos + 1
#                 end
#                 # § wpos = pos - l.start .+ 1
#                 if (1 <= wpos <= length(pile)) && (ml > mlt)
#                     pile[wpos] += fn(ml)
#                     bgpile[wpos] += 1
#                 end
#             end


#         end

#     end

#     pile, bgpile
# end

function modmeta(chroms, locations, strands, PM, modification=first(keys(PM.modifications)); mlt=-1, fn=identity, spos=true, sneg=true)
    ic = IntervalCollection(Interval.(chroms, locations, first.(strands), 1:length(chroms)))
    modmeta(ic, PM, modifciation, mlt=mlt, fn=fn, spos=spos, sneg=sneg)
end

function modmeta(ic::IntervalCollection{T}, PM, modification=first(keys(PM.modifications)); mlt=-1, fn=identity, spos=true, sneg=true) where {T}

    n = GenomicFeatures.span(first(ic))
    pile = zeros(n)
    bgpile = zeros(n)

    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ (1:2)
    
    for (ii, ri) in eachoverlap(ic, PM.readivs)
        rii = GenomicFeatures.metadata(ri)
        pind = PM.readindex[rind[1], rii]:PM.readindex[rind[2], rii]

        !spos && (strand(ri) == STRAND_POS) && continue
        !sneg && (strand(ri) == STRAND_NEG) && continue

        for pi in pind
            pos = PM.PM[1, pi]
            ml  = PM.PM[2, pi]
                
            if strand(ii) == STRAND_POS
                wpos = pos - leftposition(ii) .+ 1
            else
                wpos = rightposition(ii) - pos + 1
            end
            # § wpos = pos - l.start .+ 1
            # wpos = pos - leftposition(ii) .+ 1

            if (1 <= wpos <= length(pile)) && (ml > mlt)
                pile[wpos] += fn(ml)
                bgpile[wpos] += 1
            end
        end
    end

    pile, bgpile
end


function readindexoffset(PM, modification)
    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ 3 .+ (1:2)

    rind
end

function modmeta(ic::IntervalCollection{T}, PM::PosMatrixCompressed, modification=first(keys(PM.modifications)); mlt=-1, norm=:none, fn=identity, spos=true, sneg=true) where {T}

    n = GenomicFeatures.span(first(ic))
    pile = zeros(n)
    # bgpile = zeros(n)

    data = Int32[] ### need to make this a type parameter in PM
    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ 3 .+ (1:2)

    totalspan = 0

    for (ii, ri) in eachoverlap(ic, PM.readivs)


        rii = GenomicFeatures.metadata(ri)
        totalspan += length(intersect(leftposition(ii):rightposition(ii), PM.readindex[2, rii]:PM.readindex[3, rii]))
        pind = PM.readindex[rind[1], rii]:PM.readindex[rind[2], rii]

        # pind = GenomicFeatures.metadata(ri)
        !spos && (strand(ri) == STRAND_POS) && continue
        !sneg && (strand(ri) == STRAND_NEG) && continue

        decompress!(data, PM.PM[pind])
        # posmeth = reshape(data, PM.modfields, div(length(data), PM.modfields))

        for pi = 1:PM.modfields:length(data)
            # pos = posmeth[1, pi]
            # ml  = posmeth[2, pi]
            pos = data[pi]
            ml = data[pi + 1] #### dirty
                
            if strand(ii) == STRAND_POS
                wpos = pos - leftposition(ii) .+ 1
            else
                wpos = rightposition(ii) - pos + 1
            end
            # § wpos = pos - l.start .+ 1
            # wpos = pos - leftposition(ii) .+ 1

            if (1 <= wpos <= length(pile)) && (ml > mlt)
                pile[wpos] += fn(ml)
                # bgpile[wpos] += 1
            end
        end
    end

    if norm == :totalspan
        pile ./= totalspan
    end

    # pile, bgpile

    pile, totalspan
end

function modmetaqc(ic::IntervalCollection{T}, PM::PosMatrixCompressed, modification=first(keys(PM.modifications)); mlt=-1, norm=:none, fn=identity, spos=true, sneg=true) where {T}

    n = GenomicFeatures.span(first(ic))
    pile = zeros(n, 256)
    # bgpile = zeros(n)

    data = Int32[] ### need to make this a type parameter in PM
    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ 3 .+ (1:2)

    totalspan = 0

    for (ii, ri) in eachoverlap(ic, PM.readivs)


        rii = GenomicFeatures.metadata(ri)
        totalspan += length(intersect(leftposition(ii):rightposition(ii), PM.readindex[2, rii]:PM.readindex[3, rii]))
        pind = PM.readindex[rind[1], rii]:PM.readindex[rind[2], rii]

        # pind = GenomicFeatures.metadata(ri)
        !spos && (strand(ri) == STRAND_POS) && continue
        !sneg && (strand(ri) == STRAND_NEG) && continue

        decompress!(data, PM.PM[pind])
        # posmeth = reshape(data, PM.modfields, div(length(data), PM.modfields))

        for pi = 1:PM.modfields:length(data)
            # pos = posmeth[1, pi]
            # ml  = posmeth[2, pi]
            pos = data[pi]
            ml = data[pi + 1] #### dirty
                
            if strand(ii) == STRAND_POS
                wpos = pos - leftposition(ii) .+ 1
            else
                wpos = rightposition(ii) - pos + 1
            end
            # § wpos = pos - l.start .+ 1
            # wpos = pos - leftposition(ii) .+ 1

            if (1 <= wpos <= n) && (ml > mlt)
                pile[wpos, ml] += 1
                # bgpile[wpos] += 1
            end
        end
    end

    if norm == :totalspan
        pile ./= totalspan
    end

    # pile, bgpile

    pile, totalspan
end

function modheat(chroms, locations, strands, PM; mlt=-1, fn=identity, spos=true, sneg=true)
    ic = IntervalCollection(Interval.(chroms, locations, first.(strands), 1:length(chroms)))
    modheat(ic, PM, mlt=mlt, fn=fn, spos=spos, sneg=sneg)
end

function modheat(ic::IntervalCollection{T}, PM; mlt=-1, fn=identity, spos=true, sneg=true) where {T}

    n = GenomicFeatures.span(first(ic))
    H = zeros(n, length(ic))
    B = zeros(n, length(ic))

    for (ii, ri) in eachoverlap(ic, PM.readivs)
        pind = GenomicFeatures.metadata(ri)
        !spos && (strand(ri) == STRAND_POS) && continue
        !sneg && (strand(ri) == STRAND_NEG) && continue

        for pi in pind
            pos = PM.PM[1, pi]
            ml  = PM.PM[2, pi]
                
            if strand(ii) == STRAND_POS
                wpos = pos - leftposition(ii) .+ 1
            else
                wpos = rightposition(ii) - pos + 1
            end
            # wpos = pos - leftposition(ii) .+ 1
            if (1 <= wpos <= size(H, 1)) && (ml > mlt)
                H[wpos, GenomicFeatures.metadata(ii)] += fn(ml)
                B[wpos, GenomicFeatures.metadata(ii)] += 1
            end
        end
    end

    H, B
end



function modheat(ic::IntervalCollection{T}, PM::PosMatrixCompressed, modification=first(keys(PM.modifications)); mlt=-1, norm=:none, fn=identity, spos=true, sneg=true) where {T}

    n = GenomicFeatures.span(first(ic))
    pile = zeros(n, length(ic))
    # bgpile = zeros(n)

    data = Int32[] ### need to make this a type parameter in PM
    mi = PM.modifications[modification]
    off = (mi - 1)*2
    rind = off .+ 3 .+ (1:2)

    totalspan = 0

    for (ii, ri) in eachoverlap(ic, PM.readivs)


        rii = GenomicFeatures.metadata(ri)
        iii = GenomicFeatures.metadata(ii)
        totalspan += length(intersect(leftposition(ii):rightposition(ii), PM.readindex[2, rii]:PM.readindex[3, rii]))
        pind = PM.readindex[rind[1], rii]:PM.readindex[rind[2], rii]

        # pind = GenomicFeatures.metadata(ri)
        !spos && (strand(ri) == STRAND_POS) && continue
        !sneg && (strand(ri) == STRAND_NEG) && continue

        decompress!(data, PM.PM[pind])

        for pi = 1:PM.modfields:length(data)

            pos = data[pi]
            ml = data[pi + 1] #### dirty
                
            if strand(ii) == STRAND_POS
                wpos = pos - leftposition(ii) .+ 1
            else
                wpos = rightposition(ii) - pos + 1
            end

            if (1 <= wpos <= n) && (ml > mlt)
                pile[wpos, iii] += fn(ml)
                # bgpile[wpos] += 1
            end
        end
    end

    if norm == :totalspan
        pile ./= totalspan
    end

    # pile, bgpile

    pile, totalspan
end