using IDPSequences

@info "Test MSA..."
fasta_file = joinpath(@__DIR__, "test.fasta")
motif_file = joinpath(@__DIR__, "test.yaml")
if !isfile(motif_file)
    ELM.prepare_motifs(fasta_file, motif_file)
end
records = load_data(fasta_file, motif_file)

aligned = msa(records)

println(aligned)

msf_file = joinpath(@__DIR__, "test.msf")
save_msf(msf_file, aligned, records)
