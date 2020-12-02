using IDPSequences

fasta_file = joinpath(@__DIR__, "test1.fasta")
motif_file = joinpath(@__DIR__, "test1.yaml")
if !isfile(motif_file)
    fastas = load_fasta(fasta_file)
    records = ELM.get_data.([identifier(x) for x in fastas])
    ELM.save_yaml(motif_file, records)
end
records = load_data(fasta_file, motif_file)

aligned1, aligned2 = msa(records)

println(aligned1)
println(aligned2)
