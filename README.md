# DiarizationTools.jl

Tools for reading / writing "Rich Transcription" RTTM files and computing performances.

We started out representing the information in an RTTM (which speaker speaks when) as a largeish `BitArray`, rows indicating time frames and columns speakers.  For this representation, determining when anybody is speaking is easy, it is just
```julia
sad(x::BitArray) = any(x, 2)
```

However, it turns out that for speaker/cluster impurity evaluation, we basically have to do a big matrix-matrix multiplication, which is nice and simple to write down, but inefficient to compute.

So we added a list-of-segments-representation, which is basically the same as the original RTTM text format.

This all makes the code a bit hybrid.  Sorry about that.

## Install

If this were a proper package, you could do something like

```julia
Pkg.clone("https://github.com/davidavdav/DiarizationTools.jl.git")
```
but it isn't---you'll only get DiarizationTools

So install locally, work from the directory where you see `src/`, and make a softlink `data` pointing to the data.

## Usage

```julia
Pkg.add("DataFrames") ## these three lines are only necessary once
Pkg.add("CSV")
Pkg.add("Logging")

include("src/rttm.jl")
include("src/evalmeasures.jl")
include("src/exp.jl")
impurities(sys="xvectors", nfiles=10)
```

`nfiles` should really be 83 to do the entire analysis.

The data printed to the screen is

- system type (ivectors or xvectors)
- number of files
- threshold
- speaker impurity
- cluster impurity
- false alarm fraction (fraction of silence recognized as speech)
- missed fraction (fraction of speech recognized as silence)

The current code only evaluates on three built-in thresholds (0.0, 0.75 and 0.1542) because it is so slow.
