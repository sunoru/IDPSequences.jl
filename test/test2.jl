using IDPSequences

@info "Test 2..."
fasta_file = joinpath(@__DIR__, "test1.fasta")
motif_file = joinpath(@__DIR__, "test1.yaml")
records = load_data(fasta_file, motif_file)

aligned = msa(records; first_gapped = true)

for each in aligned
    println(each)
end

msf_file = joinpath(@__DIR__, "test2.msf")
save_msf(msf_file, aligned, records)
