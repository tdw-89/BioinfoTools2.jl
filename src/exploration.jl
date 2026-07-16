module Exploration

using Distributions
using KernelDensity
using StatsBase

using ..Reference
using ..Studies

function coverage(
    genome::Genome, 
    feature::Union{String, Symbol}, 
    data::BedData
    ; 
    filter_zeros::Bool=false
    )::Dict{String, Vector{Float64}}
    
    intersection = Studies.intersect(genome, data, feature)
    scaffolds = Dict{String, Vector{Float64}}()
    for (scaffold_name, scaffold) in genome.scaffolds
        if !haskey(intersection, scaffold_name)
            continue
        end

        left_tree = scaffold.features
        right_tree = intersection[scaffold_name]
        iter = leftjoin(left_tree, right_tree)
        frac_coverage = Float64[]
        sizehint!(frac_coverage, length(left_tree))
        for pair in iter
            if isnothing(pair[2])
                push!(frac_coverage, 0)
            else
                length_subject = pair[1].last - pair[1].first + 1
                length_object = pair[2].last - pair[2].first + 1
                @assert length_subject >= length_object
                push!(frac_coverage, length_object / length_subject)
            end
        end

        if filter_zeros 
            frac_coverage = frac_coverage |> filter(x -> x != 0)
        end

        scaffolds[scaffold_name] = frac_coverage
    end
    return scaffolds
end

function kde(genome::Genome, feature::Union{String, Symbol}, data::BedData; filter_zeros::Bool=false)
    frac_coverage = coverage(genome, feature, data, filter_zeros=filter_zeros)
    coverage_kde = Dict{String, Union{Nothing, UnivariateKDE}}()
    for k in keys(frac_coverage)
        coverage_kde[k] = isempty(frac_coverage[k]) ? nothing : KernelDensity.kde(frac_coverage[k])
    end
    return coverage_kde
end

end