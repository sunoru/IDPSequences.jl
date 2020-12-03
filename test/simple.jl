using IDPSequences
using BioAlignments

@info "Test Simple..."
fasta_file = joinpath(@__DIR__, "simple.fasta")
motif_file = joinpath(@__DIR__, "simple.yaml")
records = load_data(fasta_file, motif_file)

@info "Alignment without motif information:"
aligned1 = msa(
    records;
    motif_weight = 0,
    # one_round = true,
    gap_opening_penalty = -3.0,
    gap_extension_penalty = -1.0,
    ending_gaps_penalty = -0.2,
)
println(aligned1)

@info "Alignment with motif information:"
aligned2 = msa(
    records;
    motif_weight = 30,
    # one_round = true,
    gap_opening_penalty = -3.0,
    gap_extension_penalty = -1.0,
    ending_gaps_penalty = -0.2,
)
println(aligned2)

msf_file1 = joinpath(@__DIR__, "simple-1.msf")
msf_file2 = joinpath(@__DIR__, "simple-2.msf")
save_msf(msf_file1, aligned1, records)
save_msf(msf_file2, aligned2, records)

using Test

for (x, y) in zip(aligned1, [
    "SLPFCSVPFSIPF",
    "SLPF---------",
    "SLPF-----SIPF",
    "---------SIPF",
    "SLPF-----SIPF"
])
    @test String(x) ≡ y
end

for (x, y) in zip(aligned2, [
    "SLPFCSVPFSIPF",
    "-----SLPF----",
    "-----SLPFSIPF",
    "-----SIPF----",
    "SLPF-SIPF----"
])
    @test String(x) ≡ y
end

