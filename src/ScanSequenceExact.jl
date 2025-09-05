# ScanSequenceExact.jl 
# Only modularize it if it'll be reused 
#module ScanSequenceExact

# Add more info so that these packages are installed if they are not found on
# the users system
#

using Pkg, CodecZlib, BioSequences, FASTX, ArgParse, CSV, DataFrames

#export ScanSequence 

# Accept and manage argument using the ArgParse package https://argparsejl.readthedocs.io/en/latest/argparse.html 
# ArgParse 

function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table settings begin
        "fasta_file" 
            help = "fasta file containing RNA/DNA sequences of interest" 
            required = true
        "motif_sequence"
            help = "string of the motif that is being searched for"
            required = true
            arg_type = String
#TODO - add in another ARG for output directory

#        "--opt1"
#            help = "an option with an argument"
#        "--opt2", "-o"
#            help = "another option with an argument"
#            arg_type = Int
#            default = 0
#        "--flag1"
#            help = "an option without argument, i.e. a flag"
#            action = :store_true
#        "arg1"
#            help = "a positional argument"
#            required = true
    end

    return parse_args(settings)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end
end

main()

# Set ARGS - arguments
# https://docs.julialang.org/en/v1/manual/command-line-interface/ 

fasta_sequence = ARGS[1]
fasta_sequence_filename = basename(fasta_sequence)

if typeof(ARGS[2]) != String
    throw(ArgumentError("The motif_sequence positional argument is not a string of letters (e.g. AGTA), please verify the motif"))
else
    motif_sequence = ARGS[2]
end 

# Validate the input arguments - throw errors if input files are not FASTA
# formatted, or if they have an incorrect extension

if endswith(fasta_sequence, ".gz")
    validate_fasta(GzipDecompressorStream(open(fasta_sequence))) !== nothing ? throw(ArgumentError("The input gzipped FASTA file (first argument) is not correctly formatted in FASTA format, please quality check the file")) : nothing
elseif endswith(fasta_sequence, r".fasta|.fa")
    validate_fasta(open(fasta_sequence)) !== nothing ? throw(ArgumentError("The input file FASTA file (first argument) is not correctly formatted in FASTA format, please quality check the file")) : nothing
elseif validate_fasta(open(fasta_sequence)) === nothing  
    throw(ArgumentError("The input file file does not have a FASTA file extension, but it appears to be a correctly formatted FASTA file, please add an .fa or .fasta file extension to avoid future confusion"))
else 
    throw(ArgumentError("The input file file does not have a FASTA file extension, and it doesn't appear to be a correctly formatted FASTA file either, please look into the file you are providing"))
end 

# TODO Convert the motif sequence into a BioSequence query type - LongDNA/RNA{4} - BUT both sequence AND query HAVE TO BE the same SequenceType

if 'U' ∈ motif_sequence 
    motif_biosequence = LongRNA{4}(motif_sequence)
else 
    motif_biosequence = LongDNA{4}(motif_sequence)
end 

## Create the BioSequence query struct 

motif_query = ExactSearchQuery(motif_biosequence)

## Load the fasta_sequence file into memory 

if endswith(fasta_sequence, r".gz")
    fasta_io = FASTAReader(GzipDecompressorStream(open(fasta_sequence)))
    fasta_sequence_records = collect(fasta_io) ; close(fasta_io)
else
    fasta_io = FASTAReader((open(fasta_sequence)))
    fasta_sequence_records = collect(fasta_io) ; close(fasta_io)
end 

#//TODO
# Verify that the identifiers for each fasta record are unique, otherwise throw
# an error that they are duplicated and will cause double ups when counting
# motifs occurance 

#record_identifier_vector = identifier.(fasta_sequence_records)

#for id in record_identifier_vector

#length(identifier.(fasta_sequence_records)) != length(unique(identifier.(fasta_sequence_records))) : throw(ErrorException("The input FASTA file contains duplicate FASTA identifiers/headers, please investigate the FASTA file"))

## sequence(fasta_sequence_records) needs to be in LongDNA{} type in order to
## perform search 

# Loop through the sequences a perform a motif search, get length of sequence
# gc content and the location of the motif match, as well as a motif match count

# Create a DataFrame to store matches, count of matches, and record length, gc_content?

match_dataframe = DataFrame(record = String[], length = Int64[], gc_content = Float64[], motif_loci = Vector[], count = Int64[])

for record in fasta_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    record_length = length(sequence(record))
    gc_content = round(BioSequences.gc_content(record_sequence), sigdigits = 3)
    # Search the motif against the sequence)
    motif_search = findall(motif_query, record_sequence)
    if !isempty(motif_search)
        start_range_vector = []
        for range in motif_search
            push!(start_range_vector, range.start)
        end
        match_count = length(start_range_vector)
        push!(match_dataframe, [record_id, record_length, gc_content, start_range_vector, match_count])
    end
end 

CSV.write("$(fasta_sequence_filename).csv", match_dataframe)


## Store the location of the first base of a match, pasting into a single column
## vector

matches_vector = []

for record in fasta_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    # Search the motif against the sequence)
    motif_search = findall(motif_query, record_sequence)
    if !isempty(motif_search)
        for range in motif_search
            push!(matches_vector, range.start)
        end
    end
end 


## TODO 
## Average motif density per base-pair - in 100bp bins?


### Remove the bins with zero values, they distort the plot too much 

#filter(p -> !iszero(p.second), bin_dict)

#plot(filter(p -> !iszero(p.second), bin_dict), size=(700,200)) 


#function ScanSequenceExact(fasta, motif::AbstractString)

#end 
