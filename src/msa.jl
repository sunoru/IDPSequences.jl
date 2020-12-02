using FASTX

struct DisorderMSAModel{
    S <: AbstractSubstitutionMatrix{Float64}
} <: AbstractScoreModel{Float64}
    score_profile::ScoreProfile
    feature_profile::FeatureProfile
    submat::S
    gap_opening_penalty::Float64
    gap_extension_penalty::Float64
    ending_gaps_penalty::Float64
end

function pairwise_align(record::Record, model::DisorderMSAModel)
    n1 = size(model.score_profile.data, 2)
    seq = sequence(record)
    n2 = length(seq)
    # Three matrices: V, H, G
    V = zeros(n1 + 1, n2 + 1)
    H = zeros(n1 + 1, n2 + 1)
    G = zeros(n1 + 1, n2 + 1)

    for i = 2:n1 + 1
        V[i, 1] = G[i, 1] = -∞
        H[i, 1] = (i - 1) * model.ending_gaps_penalty
    end
    for i = 2:n2 + 1
        V[1, i] = H[1, i] = -∞
        G[1, i] = (i - 1) * model.ending_gaps_penalty
    end
end

function get_identities(records::Records, model::DisorderMSAModel)
    nrecords = length(records)

end

function msa(
    records::Records;
    gap_opening_penalty = -5.0,
    gap_extension_penalty = -1.0,
    ending_gaps_penalty = -0.1,
    motif_weight = 3,
    submat = DisorderScoringMatrix
)
    n = length(records[1])
    init_records = [records[1]]
    score_profile = ScoreProfile(init_records, submat)
    feature_profile = FeatureProfile(n, motif_weight)
    update!(feature_profile, init_records)
    model = DisorderMSAModel(
        score_profile, feature_profile, submat,
        gap_opening_penalty, gap_extension_penalty, ending_gaps_penalty
    )
    # Align all to the first
    identities = [1.0]
end

msa(filename_or_records; kwargs...) = msa(parse_fasta(filename_or_records); kwargs...)
