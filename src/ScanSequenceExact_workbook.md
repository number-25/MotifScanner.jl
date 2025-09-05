# Workbook for ScanSequenceExact.jl

Scan a sequence for an exact motif match - ScanSequenceExact.jl

## Handle arguments

This functionality will be handled with the `ArgParse.jl` [package](https://argparsejl.readthedocs.io/en/latest/argparse.html).

**Safechecks**

- Is the input file a FASTA file? Does it have the correct extension? [done]
- Ensure that the identifier for each sequence is unique, otherwise clashes will take place, and results will be combined for sequences with the same id.
  - Warn the user if the identifiers are not unique, and print which ones they are so that they can be changed manually.

## Basic FASTX and BioSequences functionality

```julia

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

## Create a DataFrame to store matches, count of matches, and record length, gc_content?

```julia

match_dataframe = DataFrame(record = String[], length = Int64[], gc_content = Float64[], motif_loci = Vector[], count = Int64[])

# Loop through the sequences a perform a motif search, get length of sequence
# gc content and the location of the motif match, as well as a motif match count

for record in fasta_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    record_length = length(sequence(record))
    gc_content = round(gc_content(LongDNA{4}(sequence(record_sequence))), sigdigits = 3)
    # Search the motif against the sequence)
    motif_search = findall(motif_query, record_sequence)
    if !isempty(motif_search)
        start_range_vector = []
        for range in motif_search
            push!(start_range_vector, range.start)
        end
        match_count = length(start_range_vector)
        push!(match_dataframe, [record_id, length, gc_content, start_range_vector, match_count])
    end
end

CSV.write("test.csv", match_dataframe)
```

---

Below hasn't been incorporated into the src yet, still testing

## Collect all exact matches into single vector

```julia

lncRNA_matches_vector = []

for record in lncRNA_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    # Search the motif against the sequence)
    motif_search = findall(PAS_query, record_sequence)
    if !isempty(motif_search)
      for range in motif_search
          push!(lncRNA_matches_vector, range.start)
      end
    end
end

# Store the number of matches at each base pair inside a dict,
# keys are base pair locus, values are number of matches at each base pair

bp_bin_dict = Dict(map(x -> x => 0, collect(1:3000)))

for match in lncRNA_matches_vector
    for value in collect(1:3000)
        if match == value
            bp_bin_dict[value] += 1
            break
        end
    end
end
```

### Normalization

Normalise the number of matches at each base pair by the number of total base
pairs at that position for all lncRNAs. This will ensure that those loci with
many matches are not inflated simply because there are more transcripts of that
size in them, since we're looking at overall motif density.

```julia

# Random, unrelated filter function
filtered_bin_dict = filter(p -> !iszero(p.second), bin_dict)

lncRNA_length_vector = []

# Create a vector of the transcript lengths
for records in lncRNA_sequence_records
    push!(lncRNA_length_vector, length(sequence(records)))
end

# Attempt to count the number of transcripts in each bp locus e.g. 1 all the way up to 3000. This should provide a figure which we can use to normalize the motif count

lncRNA_bp_length_dict = Dict(map(x -> x => 0, collect(1:3000)))

for rna in lncRNA_length_vector
    for num in 1:3000
        if num ∈ 0:rna
            lncRNA_bp_length_dict[num] += 1
        end
    end
end

# Now work through the two dicts and create a third dict containing the normalized values?
# This appears to be functioning correctly now.

normalized_lncRNA_dict = Dict()

for (key,value) in bp_bin_dict
    normalized_lncRNA_dict[key] = (value / lncRNA_bp_length_dict[key]) * 10^4
end
```

### Do a basic plot of the normalized dict

`plot(normalized_lncRNA_dict)`
