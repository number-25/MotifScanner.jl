module MotifScanner

using Pkg, CodecZlib, BioSequences, FASTX, ArgParse

# Write your package code here.

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




end


Typos to push?

# SamplerWeighted 
Weighted sampler of type T. Instantiate with a collection of eltype T containing the
elements to sample, and an **ordered** collection of probabilities to sample each element
except the last. The last probability is the remaining probability up to 1.

# src/biosequence/predicates.hl 
!!! note
Using the [`reverse_complement`](@ref) of a DNA sequence will **give** this
reverse complement. 

