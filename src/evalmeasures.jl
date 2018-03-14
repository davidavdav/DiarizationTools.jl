## evalmeasures.jj
## (c) 2018 David A. van Leeuwen
## MIT License

using DiarizationTools

## sad(x): compute speech activity detection,
## or a mask of non-null classes.
sad(x::BitMatrix) = any(x, 2)
sad(x::BitVector) = x

## impurity of bitarray hyp w.r.t. bitarray ref.
## rows represent time or data points
## For ref, columns specify different classes,
## a row is a true-hot encoding of the class id.
## For hyp, the bits indicate time or data point
## belonging to the same class
## FA time is counted as error, missed time is ignored
function impurity(ref::BitMatrix, hyp::BitVector)
    overlap = ref' * hyp
    return 1 - maximum(overlap)/sum(hyp)
end

## for a number of reference vectors, the impurity is the
## average over al hypothesized classes
## Note: this function scales very badly with the sizes of ref and hyp!
function impurity(ref::BitMatrix, hyp::BitMatrix)
    overlap = ref' * hyp
    return 1 - mean(maximum(overlap, 1) ./ sum(hyp, 1))
end

function famiss(ref::BitMatrix, hyp::BitMatrix)
    refsad = sad(ref)
    hypsad = sad(hyp)
    fa = sum(hypsad .& .!refsad) / sum(.!refsad)
    miss = sum(.!hypsad .& refsad) / sum(refsad)
    return fa, miss
end
famiss(ref::RTTM, hyp::RTTM) = famiss(ref.b, hyp.b)

spclimpurity(ref::BitArray, hyp::BitArray) = impurity(hyp, ref), impurity(ref, hyp)
# spclimpurity(ref::RTTM, hyp::RTTM) = spclimpurity(ref.b, hyp.b)

## using intervals, is this faster?

## Helper from rocanalysis
## Returns the index to the largest value in the sorted array `a` ≤ `x` if lower==false
## If lower==true, the value must be strictly < `x`
function binsearch(x::Real, a::Vector{T}; lower=false, check=false) where T<:Real
    check && (issorted(a) || error("Array needs to be sorted"))
    mi = 1
    ma = length(a)
    if x < a[mi] || lower && x == a[mi]
        return 0
    elseif x > a[ma] || !lower && x == a[ma]
        return ma
    end
    while ma - mi > 1
        h = mi + (ma - mi) ÷ 2
        if x > a[h] || !lower && x == a[h]
            mi = h
        else
            ma = h
        end
    end
    return mi
end

## in1 and in2 are matrices of fame start-stop times as rows
## they do not need to have the same length

## complexity nseg(in1) × nseg(in2)
function overlap1(in1::Matrix{T}, in2::Matrix{T}) where T
    tot_overlap = zero(T)
    for i in 1:size(in1, 1), j in 1:size(in2, 1)
        tot_overlap += max(0, min(in1[i, 2], in2[j, 2]) - max(in1[i, 1], in2[j, 1]))
    end
    return tot_overlap
end

## complexity nseg(in1) + nseg(in2)
function overlap(in1::Matrix{T}, in2::Matrix{T}) where T
    tot_overlap = zero(T)
    ni, nj = size(in1, 1), size(in2, 1)
    i = j = 1
    while i ≤ ni && j ≤ nj
        tot_overlap += max(0, min(in1[i, 2], in2[j, 2]) - max(in1[i, 1], in2[j, 1]))
        if in1[i, 2] < in2[j, 2]
            i += 1
        else
            j += 1
        end
    end
    return tot_overlap
end

## impurity between list-of-intervals and intervals
function impurity(ref::Vector{Matrix{T}}, h::Matrix{T}) where T
    overlaps = [overlap(r, h) for r in ref]
    return 1.0 - maximum(overlaps) / sum(diff(h, 2))
end

impurity(ref::Vector{Matrix{T}}, hyp::Vector{Matrix{T}}) where T = mean(impurity(ref, h) for h in hyp)

spclimpurity(ref::Vector{Matrix{T}}, hyp::Vector{Matrix{T}}) where T = impurity(hyp, ref), impurity(ref, hyp)
spclimpurity(ref::RTTM, hyp::RTTM) = spclimpurity(intervals(ref), intervals(hyp))
