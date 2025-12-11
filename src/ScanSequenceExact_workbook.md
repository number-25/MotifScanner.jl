# Workbook for ScanSequenceExact.jl

Scan a sequence for an exact motif match - ScanSequenceExact.jl

TODO
* Allow the input of a list of motifs? 

## Handle arguments

This functionality will be handled with the `ArgParse.jl` [package](https://argparsejl.readthedocs.io/en/latest/argparse.html).

**Safechecks**
- [x] Is the input file a FASTA file? Does it have the correct extension?
- [ ] Ensure that the identifier for each sequence is unique, otherwise clashes will take place, and results will be combined for sequences with the same id.
- Warn the user if the identifiers are not unique, and print which ones they are so that they can be changed manually.

## Basic FASTX and BioSequences functionality

```julia

# Load necessary packages
using Pkg, CodecZlib, BioSequences, FASTX, ArgParse, CSV, DataFrames, Plots, Random

# Create a motif sequence in BioSequence format
motif_sequence = LongDNA{4}("AGTC")

# Polyadenylation sequence motif
motif_sequence = LongDNA{4}("AATAAA")

# Create an exact query
motif_query = ExactSearchQuery(motif_sequence)

# Perform a quick match to demonstrate functionality
tmp_search = findall(motif_query, LongDNA{4}(sequence(fasta_sequence_records[2])))

```

## Create a DataFrame to store matches, count of matches, and record length, gc_content

```julia

# Loop through the sequences a perform a motif search, get length of sequence
# gc content and the location of the motif match, as well as a motif match count
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
    CSV.write("test.csv", match_dataframe)
end 

matchToDataFrame(fasta_sequence_records)

```

---

## Collect all exact matches and transcript lengths into single vector

```julia

matches_vector = reduce(vcat, match_dataframe[:,:motif_loci])
transcript_length_vector = match_dataframe[:,:length]

```


## Store the number of matches at each base pair inside a dict, keys are base pair locus, values are number of matches at each base pair

```julia

transcript_end_range = 3000 

function bpBinDict(vector)
    bp_bin_dict = Dict(map(x -> x => 0, collect(1:trnascript_end_range)))
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

bp_matches = bpBinDict(matches_vector)

```
   
## Normalization of motif matches by bp frequency at each locus

Normalise the number of matches at each base pair by the number of total base
pairs at that position for all transcripts. This will ensure that those loci with
many matches are not inflated simply because there are more transcripts of that
size in them, since we're looking at overall motif density.

Attempt to count the number of transcripts in each bp locus e.g. 1 all the way up to an end range (3000). This should provide a figure which we can use to normalize the motif count

```julia

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

bp_density = bpMotifDensity(transcript_length_vector)

```

Now work through the two dicts and create a third dict containing the normalized values

```julia

normalized_transcript_dict = Dict()

for (key,value) in bp_matches   
    normalized_transcript_dict[key] = (value / bp_density[key]) * 10^4
end

# Write out motif density at each bp in CSV format
CSV.write(normalized_transcript_dict, "", headers = []) 

```

```julia

# Random, unrelated filter function
filtered_bin_dict = filter(p -> !iszero(p.second), bp_bin_dict)

```

## Do a basic plot of the normalized dict
```julia

plot(normalized_transcript_dict, xticks = 0:300:3000, plot_title = "Motif density per base pair", seriescolor = :green, seriesalpha = 0.5, grid = false, label=false, xlabel = "Position (bp)", ylabel = "Motif density (10⁴)")

```
 
## Add the option of performing motif detection on simulated sequences 
Motif density on simulated reads - this is executed if the random_reference option is provided. 

* Simulate reads - same number of reads will be simulated as those in the main dataset/fasta
* Generate a length distribution in the simulated reads which closely follows the true reads 

### Using logNormal distribution --- find the mean and std-dev using this function

```julia

function myLogNormal(m,std)
    γ = 1+std^2/m^2
    μ = log(m/sqrt(γ))
    σ = sqrt(log(γ))

    return LogNormal(μ,σ)
end

```

```julia

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

simulateTranscripts(transcript_length_vector)

```

### Create a dataframe for the simulated transcript matches

```julia

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

simulatedMatchToDataFrame(simulateTranscripts(transcript_length_vector))

```

### Store the location of the first base of a match, pasting into a single column vector

```julia

simulated_matches_vector = reduce(vcat, simulated_match_dataframe[:,:motif_loci])
simulated_transcript_length_vector = simulated_match_dataframe[:,:length]

```

### Store the number of matches at each base pair inside a dict, keys are base pair locus, values are number of matches at each base pair plot_end_range default is 3000

```julia

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

simulated_bp_matches = simulatedbpBinDict(simulated_matches_vector)

```

### Count the number of transcripts in each bp locus e.g. 1 all the way up to end range. This should provide a figure which we can use to normalize the motif count

```julia

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

simulated_bp_density = simulatedbpMotifDensity(simulated_transcript_length_vector)

```
   
### Now work through the two dicts and create a third dict containing the normalized values

```julia

simulated_normalized_transcript_dict = Dict()

for (key,value) in simulated_bp_matches   
    simulated_normalized_transcript_dict[key] = (value / simulated_bp_density[key]) * 10^4
end

```

