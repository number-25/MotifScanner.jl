# Workbook for ScanSequenceExact.jl

Scan a sequence for an exact motif match - ScanSequenceExact.jl

TODO
* Rewrite code below to reflect functonalisation of the src code
* Allow the input of a list of motifs? 

## Handle arguments

This functionality will be handled with the `ArgParse.jl` [package](https://argparsejl.readthedocs.io/en/latest/argparse.html).

**Safechecks**

- Is the input file a FASTA file? Does it have the correct extension? [done]
- Ensure that the identifier for each sequence is unique, otherwise clashes will take place, and results will be combined for sequences with the same id.
  - Warn the user if the identifiers are not unique, and print which ones they are so that they can be changed manually.

## Basic FASTX and BioSequences functionality

```julitranscripta

# Create a motif sequence in BioSequence format
motif_sequence = LongDNA{4}("AGTC")

# Polyadenylation sequence motif
PAS_sequence = LongDNA{4}("AATAAA")

# Create an exact query
motif_query = ExactSearchQuery(motif_sequence)

PAS_query = ExactSearchQuery(PAS_sequence)

# Perform a quick match to demonstrate functionality
tmp_search = findall(motif_query, LongDNA{4}(sequence(fasta_sequence_records[2])))
```

## Create a DataFrame to store matches, count of matches, and record length, gc_content

```julia

match_dataframe = DataFrame(record = String[], length = Int64[], gc_content = Float64[], motif_loci = Vector[], count = Int64[])

# Loop through the sequences a perform a motif search, get length of sequence
# gc content and the location of the motif match, as well as a motif match count

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

CSV.write("test.csv", match_dataframe)
```

---

## Collect all exact matches and transcript lengths into single vector

```julia

matches_vector = reduce(vcat, match_dataframe[:,:motif_loci])
transcript_length_vector = match_dataframe[:,:length]

```


## Store the number of matches at each base pair inside a dict, keys are base pair locus, values are number of matches at each base pair

```julia

bp_bin_dict = Dict(map(x -> x => 0, collect(1:3000)))

for match in matches_vector
    for value in collect(1:3000)
        if match == value
            bp_bin_dict[value] += 1
            break
        end
    end
end

```

## Normalization of motif matches by bp frequency at each locus

Normalise the number of matches at each base pair by the number of total base
pairs at that position for all transcripts. This will ensure that those loci with
many matches are not inflated simply because there are more transcripts of that
size in them, since we're looking at overall motif density.

```julia

# Random, unrelated filter function
filtered_bin_dict = filter(p -> !iszero(p.second), bin_dict)

# Attempt to count the number of transcripts in each bp locus e.g. 1 all the way up to 3000. This should provide a figure which we can use to normalize the motif count

bp_length_dict = Dict(map(x -> x => 0, collect(1:3000)))

for rna in transcript_length_vector
    for num in 1:3000
        if num ∈ 0:rna
            bp_length_dict[num] += 1
        end
    end
end

# Now work through the two dicts and create a third dict containing the normalized values?
# This appears to be functioning correctly now.

normalized_transcript_dict = Dict()

for (key,value) in bp_bin_dict
    normalized_transcript_dict[key] = (value / bp_length_dict[key]) * 10^4
end

```

## Do a basic plot of the normalized dict
```julia

plot(normalized_transcript_dict, xticks = 0:300:3000, plot_title = "Motif density per base pair", seriescolor = :green, seriesalpha = 0.5, grid = false, label=false, xlabel = "Position (bp)", ylabel = "Motif density (10⁴)")

```
