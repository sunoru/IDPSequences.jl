using BioSequences

function calc_identity(s1, s2, seq1)
    n = length(s1)
    @assert length(s2) ≡ n
    gaps = 0
    identicals = 0
    @inbounds for i in 1:n
        if s1[i] ≡ AA_Gap
            gaps += 1
        elseif s2[i] ≡ seq1[i - gaps]
            identicals += 1
        end
    end
    identicals / n
end

function get_identities(sequences::SequenceList, model::DisorderMSAModel)
    seq1 = sequence(sequences[1])
    gapped = true
    identities = [if i ≡ 1
        1.0
    else
        s1, s2 = pairwise_align(seq, model, gapped)
        identity = calc_identity(s1, s2, seq1)
    end for (i, seq) in enumerate(sequences)]
end

function perform_msa_round_ungapped(
    (aligned1, aligned2)::NTuple{2, Vector{FeaturedSequence}},
    aligned_number::Int,
    sequences::SequenceList, model::DisorderMSAModel,
    identities, cutoff
)
    next_alignments = sum(identities .≥ cutoff)
    if next_alignments ≤ aligned_number
        return (aligned1, aligned2)
    end
    seq1 = sequences[1]
    aligned1 = [seq1]
    aligned2 = [seq1]
    nseq = length(sequences)
    first_gapped = false
    for i in 2:nseq
        if identities[i] ≥ cutoff
            @inbounds seq = sequences[i]
            s1, s2 = pairwise_align(seq, model, first_gapped)
            push!(aligned1, s1)
            push!(aligned2, s2)
        end
    end
    aligned1, aligned2
end

function merge_alignments(
    (aligned1, aligned2)::NTuple{2, Vector{FeaturedSequence}}
)
# TODO
    aligned1, aligned2
end

function perform_msa_round_gapped(
    (aligned1, aligned2)::NTuple{2, Vector{FeaturedSequence}},
    aligned_number::Int,
    sequences::SequenceList, model::DisorderMSAModel,
    identities, cutoff
)
    next_alignments = sum(identities .≥ cutoff)
    if next_alignments ≤ aligned_number
        return (aligned1, aligned2)
    end
    aligned1 = FeaturedSequence[]
    aligned2 = FeaturedSequence[]
    nseq = length(sequences)
    first_gapped = true
    for i in 1:nseq
        if identities[i] ≥ cutoff
            @inbounds seq = sequences[i]
            s1, s2 = pairwise_align(seq, model, first_gapped)
            push!(aligned1, s1)
            push!(aligned2, s2)
        end
    end
    aligned1, aligned2 = merge_alignments(aligned1, aligned2)
end

function msa_init(
    records::RecordList,
    gap_opening_penalty,
    gap_extension_penalty,
    ending_gaps_penalty,
    motif_weight,
    submat
)
    @assert length(records) > 0
    n = length(records[1])
    original_sequences = FeaturedSequence.(records)
    init_sequences = [original_sequences[1]]
    score_profile = ScoreProfile(init_sequences, submat)
    feature_profile = FeatureProfile(records, motif_weight)
    update!(feature_profile, init_sequences)
    model = DisorderMSAModel(
        Ref(score_profile), feature_profile, submat,
        gap_opening_penalty, gap_extension_penalty, ending_gaps_penalty
    )
    model, original_sequences
end


function msa(
    records::RecordList;
    first_gapped = false,
    gap_opening_penalty = -5.0,
    gap_extension_penalty = -1.0,
    ending_gaps_penalty = -0.1,
    motif_weight = 3,
    submat = DisorderScoringMatrix,
    one_round = false
)
    model, original_sequences = msa_init(
        records,
        gap_opening_penalty, gap_extension_penalty, ending_gaps_penalty,
        motif_weight, submat
    )
    feature_profile = model.feature_profile

    # Align all to the first
    identities = get_identities(original_sequences, model)
    aligned = (FeaturedSequence[], FeaturedSequence[])
    alignments = 0
    perform_msa_round = first_gapped ? perform_msa_round_gapped : perform_msa_round_ungapped
    if !one_round
        for i in 1:9
            cutoff = (9 - i) / 10
            aligned = perform_msa_round(aligned, alignments, original_sequences, model, identities, cutoff)
            new_alignments = length(aligned[1])
            if alignments < new_alignments
                alignments = new_alignments
                update!(feature_profile, aligned[1])
                model.score_profile[] = ScoreProfile(aligned[1], submat)
            end
        end
    end
    alignments = 0
    cutoff = 0
    aligned = perform_msa_round(aligned, alignments, original_sequences, model, identities, cutoff)
    # update!(feature_profile, aligned[1])
    # model.score_profile[] = ScoreProfile(aligned[1])
    # if optimize ...
    aligned[2]
end
