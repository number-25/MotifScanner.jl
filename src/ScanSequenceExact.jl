# ScanSequenceExact.jl 
# Only modularize it if it'll be reused - unless we need it for testing purposes? 
#module ScanSequenceExact

# Add more info so that these packages are installed if they are not found on
# the users system
#

using Pkg, CodecZlib, BioSequences, FASTX, ArgParse, CSV, DataFrames, Plots, Distributions

#export ScanSequence 

# Accept and manage argument using the ArgParse package https://argparsejl.readthedocs.io/en/latest/argparse.html 
# ArgParse 

function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "fasta_file" 
            help = "fasta file containing RNA/DNA sequences of interest" 
            required = true
        "motif_sequence"
            help = "string of the motif that is being searched for"
            required = true
            arg_type = String
        "--output_directory"
            help = "path to output directory, will be created if it doesn't exist"
            required = true
            arg_type = String
        "--plot"
            help = "produce a plot of motif frequency across transcript"
            action = :store_true
        "--plot_end_range"
            help = "end range of plot (3000bp by default)"
            arg_type = Int64
            default = 3000
            required = false

# montecarlo option to run motif detection on simulated sequences and
# plot results alongside motif frequency
        "--random_reference"
            help = "run motif detection on equal number of simulated sequences of similar length"
            required = false
            action = :store_true

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

args = parse_commandline()


# Set ARGS - arguments
# https://docs.julialang.org/en/v1/manual/command-line-interface/ 

output_directory = args["output_directory"]

# check if the directory exists, and create it if it doesn't

isdir(output_directory) ? println("Output is in $(output_directory)") : mkdir(output_directory) 

fasta_sequence = args["fasta_file"]

fasta_sequence_filename = basename(fasta_sequence)

if typeof(args["motif_sequence"]) != String
    throw(ArgumentError("The motif_sequence positional argument is not a string of letters (e.g. AGTA), please verify the motif"))
else
    motif_sequence = args["motif_sequence"]
end 

transcript_end_range = args["plot_end_range"]


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

## Ensure that all biosequences are of DNA alphabet

if 'U' ∈ motif_sequence 
    motif_biosequence = convert(LongDNA{4}, LongRNA{4}(motif_sequence))
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

function matchToDataFrame(sequence_records)
    for record in sequence_records
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
    CSV.write("$(output_directory)/$(fasta_sequence_filename).csv", match_dataframe)
end 

export matchToDataFrame

matchToDataFrame(fasta_sequence_records)

## Store the location of the first base of a match, pasting into a single column
## vector
matches_vector = reduce(vcat, match_dataframe[:,:motif_loci])
transcript_length_vector = match_dataframe[:,:length]

# Store the number of matches at each base pair inside a dict,
# keys are base pair locus, values are number of matches at each base pair
# plot_end_range default is 3000

function bpBinDict(vector)
    bp_bin_dict = Dict(map(x -> x => 0, collect(1:transcript_end_range)))
    for match in vector
        for value in collect(1:transcript_end_range)
            if match == value
                bp_bin_dict[value] += 1
                break
            end
        end
    end
    return bp_bin_dict
end 

export bpBinDict

bp_matches = bpBinDict(matches_vector)

# Count the number of transcripts in each bp locus e.g. 1 all the way up to end range. This should provide a figure which we can use to normalize the motif count

function bpMotifDensity(vector)
    bp_length_dict = Dict(map(x -> x => 0, collect(1:transcript_end_range)))
    for rna in vector
        for num in 1:transcript_end_range
            if num ∈ 0:rna
                bp_length_dict[num] += 1
            end
        end
    end
    return bp_length_dict
end    

export bpMotifDensity

bp_density = bpMotifDensity(transcript_length_vector)

# Now work through the two dicts and create a third dict containing the normalized values

normalized_transcript_dict = Dict()

for (key,value) in bp_matches   
    normalized_transcript_dict[key] = (value / bp_density[key]) * 10^4
end

# Write out motif density at each bp in CSV format
CSV.write(normalized_transcript_dict, "", headers = []) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

## Motif density on simulated reads - this is executed if the random_reference
# option is provided

## simulate reads - same number of reads will be simualted as those in the main
## dataset 

## Generate a length distribution in the simulated reads which closely follows
## the true reads 

### Using logNormal distribution --- find the mean and std-dev using this
#function

function myLogNormal(m, std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))

    return LogNormal(μ,σ)
end

export myLogNormal

function simulateTranscripts(length_vector)
    isempty(length_vector) ? throw(ArgumentError("please check the provided argument, it should be a vector of Floats")) : nothing
    ## Get the mean and std-dev of the real sequence reads
    fasta_sequence_reads_mean = Int64(round(mean(length_vector), digits = 0))
    fasta_sequence_reads_stdev = Int64(round(std(length_vector), digits = 0))
    # Seed the log normal distribution
    log_normal_distribution = myLogNormal(fasta_sequence_reads_mean, fasta_sequence_reads_stdev)
    # Vector of simulated transcript lengths
    simulated_transcript_length_vector = []
    for len in rand(log_normal_distribution, length(length_vector))
        push!(simulated_transcript_length_vector, Int64(round(len, digits = 0)))
    end
    ## Vector to store reads 
    simulated_transcripts = []
    ## Simulate reads 
    for transcript in simulated_transcript_length_vector
        push!(simulated_transcripts, randdnaseq(transcript))
    end
    return simulated_transcripts
end

export simulateTranscripts

simulateTranscripts(transcript_length_vector)


# Create a dataframe for the simulated transcript matches
simulated_match_dataframe = DataFrame(record = String[], length = Int64[], gc_content = Float64[], motif_loci = Vector[], count = Int64[])

function simulatedMatchToDataFrame(simulated_transcripts)
    counter = 1
    for transcript in simulated_transcripts
        transcript_name = "MTFSC_" * string(counter)
        transcript_length = length(transcript)
        gc_content = round(BioSequences.gc_content(transcript), sigdigits = 3)
        # Search the motif against the sequence)
        motif_search = findall(motif_query, transcript)
        if !isempty(motif_search)
            start_range_vector = []
            for range in motif_search
                push!(start_range_vector, range.start)
            end
            match_count = length(start_range_vector)
            push!(simulated_match_dataframe, [transcript_name, transcript_length, gc_content, start_range_vector, match_count])
            counter += 1
        end
        CSV.write("$(output_directory)/$(fasta_sequence_filename)_sim_standard.csv", simulated_match_dataframe)
    end
end

export simulatedMatchToDataFrame

simulatedMatchToDataFrame(simulateTranscripts(transcript_length_vector))

#### MORE TESTING - TODO

## Store the location of the first base of a match, pasting into a single column
## vector
simulated_matches_vector = reduce(vcat, simulated_match_dataframe[:,:motif_loci])
simulated_transcript_length_vector = simulated_match_dataframe[:,:length]

# Store the number of matches at each base pair inside a dict,
# keys are base pair locus, values are number of matches at each base pair
# plot_end_range default is 3000

function simulatedbpBinDict(vector)
    bp_bin_dict = Dict(map(x -> x => 0, collect(1:transcript_end_range)))
    for match in vector
        for value in collect(1:transcript_end_range)
            if match == value
                bp_bin_dict[value] += 1
                break
            end
        end
    end
    return bp_bin_dict
end 

export simulatedbpBinDict

simulated_bp_matches = simulatedbpBinDict(simulated_matches_vector)

# Count the number of transcripts in each bp locus e.g. 1 all the way up to end range. This should provide a figure which we can use to normalize the motif count

function simulatedbpMotifDensity(vector)
    bp_length_dict = Dict(map(x -> x => 0, collect(1:transcript_end_range)))
    for rna in vector
        for num in 1:transcript_end_range
            if num ∈ 0:rna
                bp_length_dict[num] += 1
            end
        end
    end
    return bp_length_dict
end    

export simulatedbpMotifDensity

simulated_bp_density = simulatedbpMotifDensity(simulated_transcript_length_vector)

# Now work through the two dicts and create a third dict containing the normalized values

simulated_normalized_transcript_dict = Dict()

for (key,value) in simulated_bp_matches   
    simulated_normalized_transcript_dict[key] = (value / simulated_bp_density[key]) * 10^4
end

## TODO
## Average motif density per base-pair - in 100bp bins?

### Remove the bins with zero values, they distort the plot too much?? Left
# code?

# If the --plot option is used, the following code will be executed, producing
# plot(s) of the motif density across the transcript(s)

if args["plot"] == true
    normalized_transcript_dict_plot = plot(normalized_transcript_dict, xticks = 0:300:3000, plot_title = "Motif density per base pair", seriescolor = :green, seriesalpha = 0.5, grid = false, label=false, xlabel = "Position (bp)", ylabel = "Motif density (10⁴)")
    png(normalized_transcript_dict_plot, "$(output_directory)/$(fasta_sequence_filename).png") 
end

# If the --random_reference option is used, the simulated sequences will be
# plotted on their own, and also overlayed against the real sequence data

## Plot on their own 

if args["plot"] == true && args["random_reference"] == true
    simulated_normalized_transcript_dict_plot = plot(simulated_normalized_transcript_dict, xticks = 0:300:3000, plot_title = "Motif density per base pair of simulated sequences", seriescolor = :orange, seriesalpha = 0.5, grid = false, label=false, xlabel = "Position (bp)", ylabel = "Motif density (10⁴)")
    png(simulated_normalized_transcript_dict_plot, "$(output_directory)/$(fasta_sequence_filename)_sim_standard.png") 
## Plot overlayed with 
p = plot(simulated_normalized_transcript_dict, xticks = 0:300:3000, plot_title = "Motif density per base pair of simulated sequences", seriescolor = :orange, seriesalpha = 0.5, grid = false, label=false, xlabel = "Position (bp)", ylabel = "Motif density (10⁴)")
plot!(p, normalized_transcript_dict, seriescolor = :green, seriesalpha = 0.5, label = "simulated")
#end
#function ScanSequenceExact(fasta, motif::AbstractString)

