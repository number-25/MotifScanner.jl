# ScanSequence.jl 
# Only modularize it if it'll be reused 
#module ScanSequence

# Add more info so that these packages are installed if they are not found on
# the users system
#
using Pkg, CodecZlib, BioSequences, FASTX #, ArgParse

export ScanSequence 

# Accept and manage argument using the ArgParse package https://argparsejl.readthedocs.io/en/latest/argparse.html 

sequence = $1
motif = $2

function ScanSequence(fasta, motif::AbstractString)
    try 
        catch(e)
    end 
    validate_fasta








end 
