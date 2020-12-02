using FASTX

load_fasta(fasta_file) = collect(open(FASTA.Reader, fasta_file))

function load_data(fasta_file, motif_file)
    motifs = ELM.load_yaml(motif_file)
    records = Record[]
    open(fasta_file) do fi
        for fasta in load_fasta(fasta_file)
            push!(records, Record(
                fasta, get(motifs, identifier(fasta), [])
            ))
        end
    end
    records
end

function Base.print(io::IO, seq::FeaturedSequence)
    if sum(seq.lowered) ≡ 0
        print(io, String(seq.sequence))
    else
        print(io, join(let c = Char(x)
            seq.lowered[i] ? lowercase(c) : c
        end for (i, x) in enumerate(seq.sequence)))
    end
end
