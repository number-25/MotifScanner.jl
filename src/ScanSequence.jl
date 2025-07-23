# ScanSequence.jl 
# Only modularize it if it'll be reused 
#module ScanSequence

# Add more info so that these packages are installed if they are not found on
# the users system
#
using Pkg, CodecZlib, BioSequences, FASTX, ArgParse

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
#
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
if typeof(ARGS[2]) != String
    throw(ArgumentError("The motif_sequence positional argument is not a string of letters (e.g. AGTA), please verify the motif"))
else
    motif_sequence = ARGS[2]
end 

# Validate the input arguments - throw errors if input files are not FASTA
# formatted, or if they have an incorrect extension

if endswith(fasta_sequence, r".gz")
    validate_fasta(GzipDecompressorStream(open(fasta_sequence))) !== nothing ? throw(ArgumentError("The input gzipped FASTA file (first argument) is not correctly formatted in FASTA format, please quality check the file")) : nothing
elseif endswith(fasta_sequence, r".fasta|.fa")
    validate_fasta(open(fasta_sequence)) !== nothing ? throw(ArgumentError("The input file FASTA file (first argument) is not correctly formatted in FASTA format, please quality check the file")) : nothing
elseif validate_fasta(open(fasta_sequence)) === nothing  
    throw(ArgumentError("The input file file does not have a FASTA file extension, but it appears to be a correctly formatted FASTA file, please add an .fa or .fasta file extension to avoid future confusion"))
else 
    throw(ArgumentError("The input file file does not have a FASTA file extension, and it doesn't appear to be a correctly formatted FASTA file either, please look into the file you are providing"))
end 

# Convert the motif sequence into a BioSequence query type - LongDNA/RNA{4}

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
    fasta_sequence_records = collect(fasta_io) ; close(io)
else
    fasta_io = FASTAReader((open(fasta_sequence))
    fasta_sequence_records = collect(io) ; close(io)
end 

## sequence(fasta_sequence_records) needs to be in LongDNA{} type in order to
# perform search 

for record in fasta_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    motif_search = findall(motif_query, 






#function ScanSequence(fasta, motif::AbstractString)

#end 
