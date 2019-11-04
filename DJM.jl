module DJM

using DataFrames, Statistics, StatsBase, DelimitedFiles
using Plots, FITSIO, CSV, FFTW

include("Utility.jl")
include("Timing.jl")
#include("Data.jl")
include("HEASARC.jl")
include("CXC.jl")
#include("XMM.jl")
include("Swift.jl")

end
