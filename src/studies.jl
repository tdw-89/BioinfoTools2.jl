module Studies

using ..Reference
using Dates
using IntervalTrees

mutable struct Data
end

mutable struct BioSample
    sample_id::String
    tissue_type::String
    species::Species
end

mutable struct AssayMethod
    name::String
    description::String
end

mutable struct Measurement
    file_path::String
    format::String
    data::Data
end

mutable struct Assay
    id::String
    type::String
    description::String
    measurement::Measurement
    biosample::BioSample
    method::AssayMethod
end

mutable struct Study
    id::String
    title::String
    date::Date
    assays::Vector{Assay}
end

end






