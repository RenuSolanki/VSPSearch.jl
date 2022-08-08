module VSPSearch

using FASTX
using BioSequences
using Scratch
using Base.Threads
using DocStringExtensions


export vsp_az, download_and_save, search_marked_orfs,
        run_blastp, run_cdhit


include("searchvsp.jl")



end
