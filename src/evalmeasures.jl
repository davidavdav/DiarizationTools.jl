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

## complexity nseg(ref) × nseg(hyp)
function overlap1(r::Matrix{T}, h::Matrix{T}) where T
    tot_overlap = zero(T)
    for i in 1:size(r, 1), j in 1:size(h, 1)
        tot_overlap += max(0, min(r[i, 2], h[j, 2]) - max(r[i, 1], h[j, 1]))
    end
    return tot_overlap
end

## complexity nseg(ref) + nseg(hyp)
function overlap(r::Matrix{T}, h::Matrix{T}) where T
    tot_overlap = zero(T)
    ni, nj = size(r, 1), size(h, 1)
    i = j = 1
    while true
        start = max(r[i, 1], h[j, 1])
        if r[i, 2] < h[j, 2]
            tot_overlap += max(0, r[i, 2] - start)
            i += 1
            i ≤ ni || break
        else
            tot_overlap += max(0, h[j, 2] - start)
            j += 1
            j ≤ nj || break
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

function sad(intervals::Vector{Matrix{T}}) where T
    ni = length(intervals)
    indices = repmat([1], ni)
    limits = [size(i, 1) for i in intervals]
    res = Array{Vector{T}}(0)
    start = minimum(intervals[i][indices[i], 1] for i in 1:ni if indices[i] ≤ limits[i])
    while True
        stop = maximum(intervals[i][indices[i], 2] for i in 1:ni if indices[i] ≤ limits[i])
        for i in 1:ni
            while indices[i] ≤ limit[i] && intervals[i][indices[i], 2] ≤ stop
                indices[i] += 1
            end
        end
        any(indices .≤ limits) || break
        newstart = minimum(intervals[i][indices[i], 1] for i in 1:ni if indices[i] ≤ limits[i])
        if newstart > stop
            res.push!([start, stop])
            start = newstart
        end
    end

end

function famiss(r::Vector{Matrix{T}}, h::Vector{Matrix{T}}) where T
end

spclimpurity(ref::Vector{Matrix{T}}, hyp::Vector{Matrix{T}}) where T = impurity(hyp, ref), impurity(ref, hyp)
spclimpurity(ref::RTTM, hyp::RTTM) = spclimpurity(intervals(ref), intervals(hyp))
