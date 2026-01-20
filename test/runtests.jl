using ScanSequenceExact
using Test

# ScanSequenceExact.jl
## matchToDataFrame() only accepts a very specific type Array{FASTX.FASTA.Record, 1})
test_vector = []
@test_throws ArgumentError matchToDataFrame("aaa")
@test_throws ArgumentError matchToDataFrame(5.0)
@test_throws ArgumentError matchToDataFrame(114)
@test_throws ArgumentError matchToDataFrame(test_vector)

@test  



@testset "ScanSequenceExact.jl" begin

end
