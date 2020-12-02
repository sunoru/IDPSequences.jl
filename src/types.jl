using FASTX
using BioSequences

struct Feature
    id::String
    range::UnitRange{Int64}
end

struct Record
    identifier::String
    description::String
    sequence::LongAminoAcidSeq
    motifs::Vector{Feature}
end
Record(fasta::FASTA.Record, motifs) = Record(
    identifier(fasta), description(fasta), sequence(fasta),
    motifs
)

identifier(r::Record) = r.identifier
description(r::Record) = r.description
sequence(r::Record) = r.sequence
Base.length(r::Record) = length(sequence(r))

const Records = AbstractVector{Record}
