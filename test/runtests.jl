include("../src/ScanSequenceExact.jl")
#include("../src/ScanSequenceFuzzy.jl")

## To run test locally, type "] test" in julia REPL

using .ScanSequenceExact
using Test

@testset "MotifScanner.jl" begin

    @testset "ScanSequenceExact.jl" begin
        include("test_ScanSequenceExact.jl")
    end

end
