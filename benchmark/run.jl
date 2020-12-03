#!/bin/env julia

using Printf
using Glob
using FASTX
using IDPSequences

function main(fasta_file, motif_file)
    @info "Aligning $fasta_file"
    records = load_data(fasta_file, motif_file)
    aligned = msa(records)
    save_msf(replace(motif_file, "yaml" => "msf"), aligned, records)
    open(replace(motif_file, ".yaml" => ".aligned"), "w") do fo
        for (each, record) in zip(aligned, records)
            @printf fo "%-27s" identifier(record)
            println(fo, each)
        end
    end
end

if abspath(PROGRAM_FILE) === @__FILE__
    if isdir(ARGS[1])
        for each in glob("$(ARGS[1])/*.in_tfa")
            main(each, replace(each, "in_tfa" => "yaml"))
        end
        for each in glob("$(ARGS[1])/*.fasta")
            main(each, replace(each, "fasta" => "yaml"))
        end
    else
        main(ARGS[1], ARGS[2])
    end
end
