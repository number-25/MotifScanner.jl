include("../src/ScanSequenceExact.jl")
#include("../src/ScanSequenceFuzzy.jl")

## To run test locally, type "] test" in julia REPL

using .ScanSequenceExact
using Test

@testset "ScanSequenceExact.jl" begin
    fasta_io = FASTAReader((open("./test_data/sample_lncRNAs_20.fasta")))
    fasta_sequence_records = collect(fasta_io) ; close(fasta_io)
    motif_query = ExactSearchQuery(LongDNA{4}("ATTAA"))
    output_directory = "."
    fasta_sequence_filename = "test_fasta"
    transcript_end_range = 3000
    @test matchToDataFrame(fasta_sequence_records) isa String
    matches_vector = reduce(vcat, match_dataframe[:,:motif_loci])
    transcript_length_vector = match_dataframe[:,:length]
    @test bpBinDict(matches_vector) isa Dict
    @test bpMotifDensity(matches_vector) isa Dict
    fasta_sequence_reads_mean = 2000
    fasta_sequence_reads_stdev = 2500
    @test myLogNormal(fasta_sequence_reads_mean, fasta_sequence_reads_stdev) isa LogNormal
    simulated_match_dataframe = DataFrame(record = String[], length = Int64[], gc_content = Float64[], motif_loci = Vector[], count = Int64[])
    @test simulateTranscripts(transcript_length_vector) isa Vector 
    @test simulatedMatchToDataFrame(simulateTranscripts(transcript_length_vector)) isa String  
    simulated_matches_vector = reduce(vcat, simulated_match_dataframe[:,:motif_loci])
    simulated_transcript_length_vector = simulated_match_dataframe[:,:length]
    @test simulatedbpBinDict(simulated_matches_vector) isa Dict 
    # Last function I believe?
    @test simulatedbpMotifDensity(simulated_transcript_length_vector) isa Dict 
end
