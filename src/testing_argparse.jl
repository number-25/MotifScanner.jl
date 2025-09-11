using Pkg, ArgParse

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
        "--output_directory"
            help = "path to output directory, will be created if it doesn't exist"
            required = true
            arg_type = String
        "--plot"
            help = "produce a plot of motif frequency across transcript"
            action = :store_true
        "--plot-end-range"
            help = "end range of plot (3000bp by default)"
            arg_type = Int64
            default = 3000
            required = false

## TODO - montecarlo option to run motif detection on simulated sequences and
# plot results alongside motif frequency
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

args = parse_commandline()

@show args["fasta_file"]

