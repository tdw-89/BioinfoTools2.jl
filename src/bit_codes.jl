module BitCodes

using GFF3

const Strand = GFF3.GenomicFeatures.Strand

"""
Canonical 2-bit strand encoding, shared by every packed metadata code in this
package (gene features in `Reference`, BED intervals in `Data`, and methylation
calls in `Methylation`).

| Code   | Strand           |
|--------|------------------|
| `0x00` | forward (`+`)    |
| `0x01` | reverse (`-`)    |
| `0x02` | both (`.`)       |
| `0x03` | unknown / NA     |

Two bits is the narrowest field any of the code layouts allocates (the
methylation payload has only 2 bits to spare), so it is the width every layout
encodes to. Wider strand fields — `Reference`/`Data` keep an 8-bit slot — simply
leave the upper bits of the field zero.

**NOTE:** `0x00` means *forward*, not *unknown*: a zeroed code decodes to `+`.
Always write a strand explicitly (use [`STRAND_NA`](@ref) when there isn't one)
rather than relying on a default-initialized code.
"""
const STRAND_FWD = UInt8(0)
const STRAND_REV = UInt8(1)
const STRAND_BOTH = UInt8(2)
const STRAND_NA = UInt8(3)

"""Width, in bits, of the canonical strand encoding."""
const STRAND_WIDTH = 2

#= Generic bit-field helpers =#

"""Mask covering the low `width` bits of `T`."""
@inline field_mask(::Type{T}, width::Integer) where {T<:Unsigned} =
    (one(T) << width) - one(T)

"""
Read the `width`-bit field based at bit `shift` (0-based) out of `code`.
"""
@inline get_field(code::T, shift::Integer, width::Integer) where {T<:Unsigned} =
    (code >> shift) & field_mask(T, width)

"""
Return `code` with the `width`-bit field based at bit `shift` (0-based) replaced
by `value`. Bits of `value` above `width` are discarded, so clamp first (see
[`clamp_field`](@ref)) when overflow should saturate rather than wrap.
"""
@inline function set_field(
    code::T,
    value::Integer,
    shift::Integer,
    width::Integer,
) where {T<:Unsigned}
    mask = field_mask(T, width)
    return (code & ~(mask << shift)) | ((T(value) & mask) << shift)
end

"""
Clamp `value` to the range representable in a `width`-bit unsigned field, so
that oversized counts saturate at the field maximum instead of wrapping into a
neighbouring field.
"""
@inline function clamp_field(value::Integer, ::Type{T}, width::Integer) where {T<:Unsigned}
    value <= 0 && return zero(T)
    mask = field_mask(T, width)
    return value >= mask ? mask : T(value)
end

#= Strand conversions =#

"""Encode a `GenomicFeatures.Strand` as its canonical 2-bit code."""
function strand_code(strand::Strand)
    if strand == GFF3.GenomicFeatures.STRAND_POS
        STRAND_FWD
    elseif strand == GFF3.GenomicFeatures.STRAND_NEG
        STRAND_REV
    elseif strand == GFF3.GenomicFeatures.STRAND_BOTH
        STRAND_BOTH
    else
        STRAND_NA
    end
end

"""Encode a strand character (`+`, `-`, `.`) as its canonical 2-bit code."""
function strand_code(strand::Char)
    if strand == '+'
        STRAND_FWD
    elseif strand == '-'
        STRAND_REV
    elseif strand == '.'
        STRAND_BOTH
    else
        STRAND_NA
    end
end

"""Decode a canonical 2-bit strand code back to a `GenomicFeatures.Strand`."""
function decode_strand(code::Integer)
    bits = UInt8(code & 0x03)
    if bits == STRAND_FWD
        GFF3.GenomicFeatures.STRAND_POS
    elseif bits == STRAND_REV
        GFF3.GenomicFeatures.STRAND_NEG
    elseif bits == STRAND_BOTH
        GFF3.GenomicFeatures.STRAND_BOTH
    else
        GFF3.GenomicFeatures.STRAND_NA
    end
end

"""Render a canonical 2-bit strand code as a character (`+`, `-`, `.`, `?`)."""
function strand_char(code::Integer)
    bits = UInt8(code & 0x03)
    return bits == STRAND_FWD ? '+' :
           bits == STRAND_REV ? '-' : bits == STRAND_BOTH ? '.' : '?'
end

"""Convert a strand character (`+`, `-`, `.`) to a `GenomicFeatures.Strand`."""
get_strand(strand::Char) = decode_strand(strand_code(strand))

export STRAND_FWD,
    STRAND_REV,
    STRAND_BOTH,
    STRAND_NA,
    STRAND_WIDTH,
    clamp_field,
    decode_strand,
    field_mask,
    get_field,
    get_strand,
    set_field,
    strand_char,
    strand_code

end
