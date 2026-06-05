# Workbook for MotifScanner.jl 

Notes, explanations, summaries and general information on the development of this package. 

* https://julialang.org/blog/2021/01/precompile_tutorial/ 
* https://github.com/GunnarFarneback/PackageCompatUI.jl

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

## Handle arguments 

This functionality will be handled with the `ArgParse.jl` [package](https://argparsejl.readthedocs.io/en/latest/argparse.html). 


# TODO
* Write tests and test set for ScanSequenceExact.jl


### Perform unbiased enrichment of motif - MotifEnrichment.jl
* PWM/enrichment using a vector of motifs of size k

sample function for unbiased enrichment

```julia
function seqSlider(se, window::Int64)
    motif_vector = []
    end_index = lastindex(se)
    counter = 1
    rev_se = reverse(se)
    while (counter + window) <= end_index
      push!(motif_vector, se[counter:counter+window-1])
      push!(motif_vector, rev_se[counter:counter+window-1])
      counter += 1
    end
    
    return unique(motif_vector)
end

motif_vector = []

for vec in fasta_sequence_records
    push!(motif_vector, seqSlider(sequence(vec), 6))
end

using StatsBase

motif_vector_counted = countmap(reduce(vcat, motif_vector))

filtered_motif = filter(x -> x.second >= 20, collect(motif_vector_counted))
sort!(filtered_motif, by=last)

plot([x.first for x in filtered_motif], [x.second for x in filtered_motif], xtick=false, yscale=:log)

# Do same thing on a random sequence - see if there is an enrichment of certain motifs above expectation

using Combinatorics

saturation of motifs? - at what length do all motifs occur in a string? assuming base frequency is equal - likelihood of seeing 20 ACGA motifs in string of length P … etc

collect(permutations("acgt"))
join.(collect(permutations("acgt")))

# Create vector from permutation where all elements are "DNA" type
[LongDNA{4}(x) for x in join.(collect(permutations("acgt")))] 

# as ExactSearchQuery which can be used to match a biosequence 
[ExactSearchQuery(LongDNA{4}(x)) for x in join.(collect(permutations("acgt")))] 

```


* Unbiased PWM/enrichment using all possible nmers of size k 
* Normalize for nucleotide content e.g. A/T rich?

### Montecarlo of motif on randomly generated sequence of size k

Original idea taken from [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html)


