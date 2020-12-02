#!/bin/env julia

using Glob
using IDPSequences

function main(fasta_file, motif_file)
    @info fasta_file
    records = load_data(fasta_file, motif_file)
    aligned = msa(records)
    save_msf(replace(motif_file, "yaml" => "msf"))
end

if abspath(PROGRAM_FILE) === @__FILE__
    if isdir(ARGS[1])
        for each in glob("$(ARGS[1])/*.in_tfa")
            main(each, replace(each, "in_tfa" => "yaml"))
        end
    else
        main(ARGS[1], ARGS[2])
    end
end
