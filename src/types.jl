## types.jl
## (c) 2014, 2018 David A. van Leeuwen
## MIT license

const sr = 100

## convert an RTTM to a representation with only intervals
function intervals(b::BitMatrix)
    res = Array{Matrix{Int}}(0)
    for j in 1:size(b, 2)
        d = diff(vcat([0], map(Int, b[:, j]), [0]))
        starts = (find(d .== 1) .- 1)
        stops = (find(d .== -1) .- 1)
        l = length(starts)
        @assert l == length(stops)
        push!(res, hcat(starts, stops))
    end
    return res
end

type RTTM
    b::BitArray{2}
    names::Vector
    file::String
    intervals::Vector{Matrix{Int}}
    function RTTM(b::BitMatrix, names::Vector, file::AbstractString)
        new(b, names, file, intervals(b))
    end
end
