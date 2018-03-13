## types.jl
## (c) 2014, 2018 David A. van Leeuwen
## MIT license

const sr = 100

type RTTM
    b::BitArray{2}
    names::Vector
    file::String
end
