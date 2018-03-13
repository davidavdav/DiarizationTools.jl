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

## in1 and in2 are matrices of fame start-stop times as rows
## they do not need to have the same length
function overlap(in1::Matrix{T}, in2::Matrix{T}) where T
    o = zero(T)
    for i in 1:size(in1, 1), j in 1:size(in2, 1)
        o += max(0, min(in1[i,2], in2[j,2]) - max(in1[i,1], in2[j,1]))
    end
    return o
end

## impurity between list-of-intervals and intervals
function impurity(ref::Vector{Matrix{T}}, h::Matrix{T}) where T
    overlaps = [overlap(r, h) for r in ref]
    return 1.0 - maximum(overlaps) / sum(diff(h, 2))
end

impurity(ref::Vector{Matrix{T}}, hyp::Vector{Matrix{T}}) where T = mean(impurity(ref, h) for h in hyp)

spclimpurity(ref::Vector{Matrix{T}}, hyp::Vector{Matrix{T}}) where T = impurity(hyp, ref), impurity(ref, hyp)
spclimpurity(ref::RTTM, hyp::RTTM) = spclimpurity(intervals(ref), intervals(hyp))
