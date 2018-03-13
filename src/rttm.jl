## rttm.jl
## (c) 2014, 2018 David A. van Leeuwen
## MIT license

## Continued from git repository vrt, modernized for Julia-0.6.

module DiarizationTools

export RTTM, readrttm, writerttm, files, speakers, intervals, perfectsad, rmspeaker, mosttalkative, rmsadfromspeaker, insertspeaker, splitsad, readuem, writeuem

using DataFrames
using CSV
using Logging

include("types.jl")

Base.copy(r::RTTM) = RTTM(copy(r.b), copy(r.names), copy(r.file))
Base.length(r::RTTM) = size(r.b, 1)
speakers(r::RTTM) = r.names
intervals(r::RTTM) = r.intervals

## constructor from dataframe
function RTTM(df::DataFrame)
    df = df[df[:kind] .== "SPEAKER",:]
    files = unique([s for s=df[:file]])
    length(files)==1 || error("RTTM type is for a single file")
    speakers = convert(Vector{String}, sort(unique(df[:name])))
    nsp = length(speakers)
    dur = maximum(df[:start] + df[:dur])
    b = falses(round(Int, sr*dur), nsp)
    d = Dict(s => i for (i, s) in enumerate(speakers))
    for row in eachrow(df)
        if  row[:kind] == "SPEAKER"
            start = round(Int, row[:start] * sr) + 1
            stop = start + round(Int, row[:dur] * sr) - 1
            b[start:stop, d[row[:name]]] = true
        end
    end
    RTTM(b, speakers, files[1])
end
## from IO, string, whatever CSV.read accepts:
RTTM(file) = RTTM(readrttm(file))

## export to dataframe
import DataFrames.DataFrame
function DataFrame(r::RTTM)
    df = Array{DataFrame}(0)
    for (j, sp) in enumerate(speakers(r))
        d = diff(vcat([0], map(Int, r.b[:, j]), [0]))
        starts = (find(d .== 1) .- 1)/sr
        stops = (find(d .== -1) .- 1)/sr
        l = length(starts)
        @assert l == length(stops)
        info(l)
        nas = repmat([missing], l)
        push!(df, DataFrame(t=repmat(["SPEAKER"], l), file=repmat([r.file], l),
                          chan=ones(Int, l),
                          start=starts, dur=stops-starts, ortho=nas,
                          stype=nas, name=repmat([sp], l), conf=nas, slat=nas))
#        for (start,stop) in zip(starts,stops)
#            println(@sprintf("SPEAKER %s 1 %4.2f %4.2f <NA> <NA> %s <NA> <NA>", r.file, start, start+stop, sp[j]))
#        end
    end
    vcat(df...)
end


## create RTTM object from one file in overall DataFrame
RTTM(df::DataFrame, file::String) = RTTM(df[:(file .== $file),:])

rttmheader = [:kind, :file, :chan, :start, :dur, :ortho, :stype, :name, :conf, :slat]

## reads a rttm
function readrttm(file, skipinfo=true)
    rttm = CSV.read(file, delim=' ', header=false, null="<NA>")
    names!(rttm, rttmheader[1:ncol(rttm)])
    if skipinfo
        rttm = rttm[rttm[:kind] .== "SPEAKER",:]
    end
    rttm[:file] = map(string, rttm[:file])
    categorical!(rttm, :file)         # as.factor()
    rttm
end

## using writetable is not flexible enoough and shout_cluster is picky
function writerttm(file::String, df::DataFrame; mindur=0.2)
    fd = open(file, "w")
    for i in 1:size(df,1) if df[i,:dur] > mindur
        println(fd, @sprintf("%s %s %d %4.2f %4.2f <NA> <NA> %s <NA>",
                             df[i,1], df[i,2],
                             df[i,3], df[i,4], df[i,5], df[i,8]))
    end end
    close(fd)
end
writerttm(file::String, r::RTTM; mindur=0.2) = writerttm(file, DataFrame(r))

files(rttm::DataFrame) = unique(rttm[:file])

perfectsad(r::RTTM) = RTTM(any(r.b,2), ["SPEECH"], r.file)

function rmspeaker(r::RTTM, speaker::String)
    sp = speakers(r)
    d = Dict(sp, 1:length(sp))
    cols = setdiff(1:length(sp), d[speaker])
    RTTM(r.b[:,cols], r.names[cols], r.file)
end
rmspeaker(r::RTTM, i::Int) = rmspeaker(r, r.names[i])

function rmspeaker(rttm::DataFrame, speaker::String)
    rttm[:(name .!= $speaker),:]
end

function mosttalkative(r::RTTM)
    time = sum(r.b, 1)
    mi = indmax(time)
    info(@sprintf("Most talkative speaker %s, %2.0f%% of the time",
                         r.names[mi],100time[mi]/sum(time)))
    mi
end

function rmsadfromspeaker(sad::RTTM, id::Int, ref::RTTM)
    new = copy(sad)
    common = 1:minimum(map(length, (new,ref)))
    new.b[common, 1] &= !ref.b[common, id]
    new
end


function insertspeaker(hyp::RTTM, id::Int, ref::RTTM)
    common = 1:minimum(map(length, (hyp, ref)))
    b = hcat(hyp.b, falses(size(hyp.b, 1)))
    b[common, size(b,2)] = ref.b[common, id]
    RTTM(b, vcat(hyp.names, ref.names[id]), hyp.file)
end

## take a (SAD) RTTM, and split it in partitions less than 20 minutes.
## returns the split times as a Data frame of N begin and end points, similar to an UEM file
function splitsad(r::RTTM, dur=20*60*sr)
    s = perfectsad(r)
    nseg = iceil(sum(s.b) / dur) # number of segments
    segs = Array(Float64, nseg, 2)
    cs = cumsum(s.b)
    dt = cs[end]/nseg
    i=0
    for j = 1:length(cs)
        if cs[j] > i*dt
            i += 1
            segs[i,1] = j / sr
            if i > 1
                segs[i-1,2] = (j-1) / sr
            end
        end
    end
    segs[end,2] = length(cs) / sr
    DataFrame(file=rep(r.file, nseg), chan=1, start=segs[:,1], stop=segs[:,2])
end

## filter out areas not indicated by the uem
function Base.filter(r::RTTM, uem::DataFrame)
    mask = falses(size(r.b,1),1)
    for i in 1:size(uem,1) if uem[i,1] == r.file
        start = round(uem[i,"start"]*sr) + 1
        stop = round(uem[i,"stop"]*sr) - 1
        mask[start:stop] = true
    end end
    res = copy(r)
    for i in 1:size(res.b,2)
        res.b[:,i] &= mask
    end
    res
end

function writeuem(file::String, uem::DataFrame)
    fd = open(file, "w")
    for i= 1:size(uem,1)
        println(fd, @sprintf("%s %d %8.2f %8.2f", uem[i,1], uem[i,2], uem[i,3], uem[i,4]))
    end
    close(fd)
end

function writeuem(file::String, uem::DataFrame, split::Bool)
    if split
        for i=1:size(uem,1)
            writeuem(string(file, ".s", i), uem[i,:])
        end
    else
        writeuem(file, uem)
    end
end

function readuem(file::String)
    uem = CSV.read(file, delim=' ', header=false)
    names!(uem, ["file", "chan", "start", "stop"])
    uem
end

end
