# Workbook for MotifScanner.jl 

Notes, explanations, summaries and general information on the development of this package. 

## Primary Goal
I want to be able to input the fasta sequence of an RNA, and scan the RNA
string to search for the presence of specific motifs. In particular, I am
foremost interested in the U1 SNRP and poly-adenylation signal motifs.  

For every RNA input, the position of these motifs within the RNA is tabulated
(dict?), along with the length of the RNA string, and perhaps a length
normalization quotient. From here, it would be nice to plot the presence of the
motif relative to the length, allowing for the comparison of all the RNAs based
on a common scale of [start:end].   

### Download test fasta data 
Got gencode lncRNA fasta sequences from Release 47 (GRCh38.p14), subsample 20 sequences and formatted them in multifasta format

```bash
seqtk sample gencode.v47.lncRNA_transcripts.fa 20 | seqkit seq -w 60 - > ~/local_analysis/MotifScanner/test/test_data/sample_lncRNAs_20.fasta

# create a gzipped form also, as data should be able to be provided in gzip format
gzip -k sample_lncRNAs_20.fasta
```

### Understand how BioSequences, FASTX IO functions 

#### Store FASTA data inside vector
The struct which the FASTX package employs is called a Record - each Record
represents a fasta sequence; it's *identifier* which is the first word after
'>' e.g. ">chr1", the *description*, which is the remainder of the 'header'
after the identifier e.g. ">chr1 | GHH18.ENSEM09912", and the *sequence* on the
next line, which is the actual sequence content of the record. 

The vector will thus be one consisting of FASTX Records which can be accessed in a indexed manner.   

Using a GZIP decoder, the gzipped fasta file can be validated, and accessed. 

TODO
* What if the fasta file given is massive e.g. over a gig? That's pretty damn cooked though.... all of the transcripts in ensemble 112 are only 80M...  

```julia

using CodecZlib 

validate_fasta(open("test/test_data/sample_lncRNAs_20.fasta")) === nothing
validate_fasta(GzipDecompressorStream(open("test/test_data/sample_lncRNAs_20.fasta.gz"))) === nothing

io = FASTAReader(GzipDecompressorStream(open("test/test_data/sample_lncRNAs_20.fasta.gz")))

# On QIMR rig
io = FASTAReader(GzipDecompressorStream(open("/mnt/hdd1/references/annotations/transcriptome/human/gh38/ensemble/gencode.v46_112.transcripts.fa.gz")))

fasta_sequence_records = collect(io) ; close(io)

lncRNA_sequence_records = collect(io) ; close(io)

sequence(fasta_sequence_records[1])

identifier(fasta_sequence_records[1]) 

sequence(LongDNA{2},fasta_sequence_records[1])
```

## Create src 

### Scan a sequence for an exact motif match - ScanSequenceExact.jl

#### Handle arguments 

This functionality will be handled with the `ArgParse.jl` [package](https://argparsejl.readthedocs.io/en/latest/argparse.html). 

**Safechecks** 
* Is the input file a FASTA file? Does it have the correct extension? [done] 
* Ensure that the identifier for each sequence is unique, otherwise clashes will take place, and results will be combined for sequences with the same id. 
  * Warn the user if the identifiers are not unique, and print which ones they are so that they can be changed manually. 

```julia

# Create a motif sequence in BioSequence format
motif_sequence = LongDNA{4}("AGTC")

# Create an exact query 
motif_query = ExactSearchQuery(motif_sequence)

# Perform a quick match
tmp_search = findall(motif_query, LongDNA{4}(sequence(fasta_sequence_records[2])))

# Create a DataFrame to store matches, count of matches, and record length, gc_content?

match_dataframe = DataFrame(record = "", length = [Int64], gc_content = [Float64], motif_loci = [], count = [Int64])

# Loop through the sequences a perform a motif search, get length of sequence
# gc content and the location of the motif match, as well as a motif match count

for record in fasta_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    record_length = length(sequence(record))
    gc_content = round(gc_content(LongDNA{2}(sequence(record_sequence))), sigdigits = 3)
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

# Perform this on all gencode lncRNA transcripts
```julia

# Polyadenylation sequence motif
PAS_sequence = LongDNA{4}("AATAAA")
PAS_query = ExactSearchQuery(PAS_sequence)

# Collect all exact matches into single vector  

lncRNA_matches_vector = []

for record in lncRNA_sequence_records
    record_sequence = LongDNA{4}(sequence(record))
    record_id = identifier(record)
    #record_length = length(sequence(record))
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

Normalise the number of matches at each base pair by the number of total base
pairs at that position for all lncRNAs. This will ensure that those loci with
many matches are not inflated simply because there are more transcripts of that
size in them, since we're looking at overall motif density.

```julia

# Filter function
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

Do a basic plot.

`plot(normalized_lncRNA_dict)`   





//TODO
### Perform unbiased enrichment of motif - MotifEnrichment.jl
* PWM/enrichment using a vector of motifs of size k
* Unbiased PWM/enrichment using all possible nmers of size k 
* Normalize for nucleotide content e.g. A/T rich?



### Montecarlo of motif on randomly generated sequence of size k

Original idea taken from [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html)


