module DJM

using DataFrames, Statistics, StatsBase, DelimitedFiles
using Plots, FITSIO, CSV, FFTW
using Formatting #Kirk added

include("Utility.jl")
include("Timing.jl")
#include("Data.jl")
include("HEASARC.jl")
include("CXC.jl")
#include("XMM.jl")
include("Swift.jl")

#KIRK ADDED
include("Probability.jl")
include("MatrixGen.jl")

end
