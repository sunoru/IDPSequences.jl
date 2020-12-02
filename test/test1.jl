using IDPSequences

@info "Test 1..."
fasta_file = joinpath(@__DIR__, "test1.fasta")
motif_file = joinpath(@__DIR__, "test1.yaml")
if !isfile(motif_file)
    ELM.prepare_motifs(fasta_file, motif_file)
end
records = load_data(fasta_file, motif_file)

aligned = msa(records)

for each in aligned
    println(each)
end

msf_file = joinpath(@__DIR__, "test1.msf")
save_msf(msf_file, aligned, records)
