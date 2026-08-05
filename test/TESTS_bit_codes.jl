using BioinfoTools2.BitCodes
using GFF3
using Test

const GF = GFF3.GenomicFeatures

@testset "BitCodes" begin

    # -------------------------------------------------------------------------
    @testset "strand code constants" begin
        @test STRAND_FWD == 0
        @test STRAND_REV == 1
        @test STRAND_BOTH == 2
        @test STRAND_NA == 3
        @test STRAND_WIDTH == 2
        # The canonical strand codes are all 2-bit UInt8 values.
        @test (STRAND_FWD, STRAND_REV, STRAND_BOTH, STRAND_NA) isa NTuple{4,UInt8}
    end

    # -------------------------------------------------------------------------
    @testset "field_mask" begin
        @test field_mask(UInt64, 0) == 0x0000000000000000
        @test field_mask(UInt64, 8) == 0x00000000000000ff
        @test field_mask(UInt32, 2) == 0x00000003
        @test field_mask(UInt16, 16) == typemax(UInt16)
        @test field_mask(UInt64, 32) == 0x00000000ffffffff
        @test field_mask(UInt64, 64) == typemax(UInt64)
        # The mask type follows the requested unsigned type.
        @test field_mask(UInt32, 8) isa UInt32
    end

    # -------------------------------------------------------------------------
    @testset "get_field / set_field" begin
        # A value written into a field reads back unchanged.
        code = set_field(UInt64(0), UInt32(0xdeadbeef), 0, 32)
        @test get_field(code, 0, 32) == 0xdeadbeef

        # Bits above the field width are discarded on write.
        @test get_field(set_field(UInt64(0), 0xfff, 0, 8), 0, 8) == 0xff

        # Writing one field leaves its neighbours untouched.
        packed = UInt64(0)
        packed = set_field(packed, 0x11, 0, 8)
        packed = set_field(packed, 0x22, 8, 8)
        packed = set_field(packed, 0x33, 16, 8)
        @test get_field(packed, 0, 8) == 0x11
        @test get_field(packed, 8, 8) == 0x22
        @test get_field(packed, 16, 8) == 0x33

        # Overwriting a field clears its old bits rather than OR-ing into them.
        packed = set_field(packed, 0x00, 8, 8)
        @test get_field(packed, 8, 8) == 0x00
        @test get_field(packed, 0, 8) == 0x11
        @test get_field(packed, 16, 8) == 0x33
    end

    # -------------------------------------------------------------------------
    @testset "clamp_field" begin
        @test clamp_field(-5, UInt8, 8) == 0
        @test clamp_field(0, UInt8, 8) == 0
        @test clamp_field(100, UInt8, 8) == 100
        @test clamp_field(255, UInt8, 8) == 255
        # Oversized values saturate at the field maximum instead of wrapping.
        @test clamp_field(256, UInt16, 8) == 255
        @test clamp_field(70_000, UInt32, 16) == 65535
        @test clamp_field(5, UInt32, 16) == 5
        @test clamp_field(-1, UInt32, 16) == 0
        # The result type follows the requested unsigned type.
        @test clamp_field(10, UInt32, 16) isa UInt32
    end

    # -------------------------------------------------------------------------
    @testset "strand_code" begin
        @test strand_code(GF.STRAND_POS) == STRAND_FWD
        @test strand_code(GF.STRAND_NEG) == STRAND_REV
        @test strand_code(GF.STRAND_BOTH) == STRAND_BOTH
        @test strand_code(GF.STRAND_NA) == STRAND_NA

        @test strand_code('+') == STRAND_FWD
        @test strand_code('-') == STRAND_REV
        @test strand_code('.') == STRAND_BOTH
        @test strand_code('?') == STRAND_NA
        @test strand_code('x') == STRAND_NA
    end

    # -------------------------------------------------------------------------
    @testset "decode_strand" begin
        @test decode_strand(STRAND_FWD) == GF.STRAND_POS
        @test decode_strand(STRAND_REV) == GF.STRAND_NEG
        @test decode_strand(STRAND_BOTH) == GF.STRAND_BOTH
        @test decode_strand(STRAND_NA) == GF.STRAND_NA

        # Only the low two bits matter; higher bits are ignored.
        @test decode_strand(0xff) == GF.STRAND_NA
        @test decode_strand(0xfc) == GF.STRAND_POS

        # strand_code and decode_strand are inverse over the four strands.
        strands = (GF.STRAND_POS, GF.STRAND_NEG, GF.STRAND_BOTH, GF.STRAND_NA)
        @test all(decode_strand(strand_code(s)) == s for s in strands)
    end

    # -------------------------------------------------------------------------
    @testset "strand_char" begin
        @test strand_char(STRAND_FWD) == '+'
        @test strand_char(STRAND_REV) == '-'
        @test strand_char(STRAND_BOTH) == '.'
        @test strand_char(STRAND_NA) == '?'
        # Low two bits only.
        @test strand_char(0xff) == '?'
        # Round-trips against strand_code for the three real strand characters.
        @test all(strand_char(strand_code(c)) == c for c in ('+', '-', '.'))
    end

    # -------------------------------------------------------------------------
    @testset "get_strand" begin
        @test get_strand('+') == GF.STRAND_POS
        @test get_strand('-') == GF.STRAND_NEG
        @test get_strand('.') == GF.STRAND_BOTH
        @test get_strand('?') == GF.STRAND_NA
    end

    # -------------------------------------------------------------------------
    @testset "pack_metadata / parse_* round-trip" begin
        code = pack_metadata(UInt32(7), STRAND_REV, UInt16(0x0141))
        @test parse_index(code) == UInt32(7)
        @test parse_strand(code) == GF.STRAND_NEG
        @test parse_so_term(code) == UInt16(0x0141)
        @test parse_strand(code) isa GF.Strand

        # A zeroed code decodes to index 0, forward strand (0x00 == +) and SO
        # term 0 — a reminder that a default-initialised code is *not* "unknown".
        @test parse_index(UInt64(0)) == 0
        @test parse_strand(UInt64(0)) == GF.STRAND_POS
        @test parse_so_term(UInt64(0)) == 0

        # Field maxima are representable and independent of one another.
        maxed = pack_metadata(typemax(UInt32), STRAND_BOTH, typemax(UInt16))
        @test parse_index(maxed) == typemax(UInt32)
        @test parse_strand(maxed) == GF.STRAND_BOTH
        @test parse_so_term(maxed) == typemax(UInt16)
        # The top 8 bits are reserved and must stay clear.
        @test (maxed >> 56) == 0
    end

    # -------------------------------------------------------------------------
    @testset "pack_metadata round-trips over many values" begin
        indices = UInt32[0, 1, 42, 1000, 0x7fffffff, typemax(UInt32)]
        strands = UInt8[STRAND_FWD, STRAND_REV, STRAND_BOTH, STRAND_NA]
        so_terms = UInt16[0, 1, 0x0141, 0x00ff, typemax(UInt16)]

        # Every (index, strand, so_term) combination, packed then unpacked.
        codes = [
            (i, s, t, pack_metadata(i, s, t)) for i in indices, s in strands, t in so_terms
        ]

        @test all(parse_index(c) == i for (i, s, t, c) in codes)
        @test all(parse_strand(c) == decode_strand(s) for (i, s, t, c) in codes)
        @test all(parse_so_term(c) == t for (i, s, t, c) in codes)
        # No combination ever disturbs the reserved top byte.
        @test all((c >> 56) == 0 for (i, s, t, c) in codes)
    end
end
