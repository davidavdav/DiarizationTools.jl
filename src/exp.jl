## exp.jl
## analyze the experiment

using DiarizationTools
import Logging

## merge RTTMs, keeping track of the same names, and simply concatenating all speaking time
function Base.merge(rttms::RTTM...)
    allspeakers = union([speakers(rttm) for rttm in rttms]...)
    index = Dict(s => i for (i, s) in enumerate(allspeakers))
    N = sum(size(rttm.b, 1) for rttm in rttms)
    res = falses(N, length(allspeakers))
    offset = 0
    for rttm in rttms
        b = rttm.b
        n = size(b, 1)
        for (i, s) in enumerate(speakers(rttm))
            res[offset+(1:n), index[s]] = b[:, i]
        end
        offset += n
    end
    return RTTM(res, allspeakers, "merged")
end

## fixspeaker makes speakers named "spreker-$n" unique tot the file,
## because otherwise they can be merged with similarly named speakers in other rttms
function fixspeaker!(r::RTTM)
    s = speakers(r) ## array of String
    for i in 1:length(s)
        if startswith(s[i], "spreker-")
            s[i] = r.file * s[i]
        end
    end
    return r
end

## just list the reference file names, ref and hyp need to be stacked in order
function getbasenames(dir)
    names = Vector{String}(0)
    for (root, dirs, files) in Base.Filesystem.walkdir(dir)
        for file in files
            if endswith(file, ".rttm")
                push!(names, file)
            end
        end
    end
    return names
end

function readrttms(dir, files)
    rttms = Vector{RTTM}(0)
    for file in files
        rttm = RTTM(readrttm(joinpath(dir, file)))
        fixspeaker!(rttm)
        push!(rttms, rttm)
    end
    return rttms
end

## make sure that ref and hyp are the same length, because the last segment may end at a different time
function ensuresamelength!!(ref::RTTM, hyp::RTTM)
    nref, nhyp = length(ref), length(hyp)
    n = max(nref, nhyp)
    if n > nref
        ref.b = vcat(ref.b, falses(n - nref, size(ref.b, 2)))
    elseif n > nhyp
        hyp.b = vcat(hyp.b, falses(n - nhyp, size(hyp.b, 2)))
    end
end

function impurities(;sys="ivectors", nfiles=5)
    out = open("performance.txt", "a")
    basenames = getbasenames(joinpath("data", "ground_truth"))[1:nfiles]
    root = joinpath("data", "linking_output", sys)
    for dir in ("thres0.0", "thres0.75", "thres1.1542")
            m = match(r"thres([\d.]+)", dir)
            if m != nothing
                thres = parse(Float64, m.captures[1])
                Logging.info(@sprintf("Running threshold %f", thres))
                refs = readrttms(joinpath("data", "ground_truth"), basenames) ## re-init because we run ensuresamelength!!
                hyps = readrttms(joinpath(root, dir), basenames)
                for (ref, hyp) in zip(refs, hyps)
                    ensuresamelength!!(ref, hyp)
                end
                Logging.info("Merging")
                ref = merge(refs...)
                hyp = merge(hyps...)
                Logging.info("Computing impurities")
                spi, cli = spclimpurity(ref, hyp)
                fa, miss = famiss(ref, hyp)
                line = @sprintf("%s %d %6.4f %6.4f %6.4f %6.4f %6.4f", sys, nfiles, thres, spi, cli, fa, miss)
                println(line)
                println(out, line)
        end
    end
    close(out)
end
