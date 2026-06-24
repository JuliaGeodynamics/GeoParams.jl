using Test
using GeoParams, Unitful
import GeoParams: fastpow, pow_check, find_ind, find_max_tuple, max_length_tuple,
    ntuple_idx, nreduce, value_and_partial, make_tuple, ptr2string

@testset "Utils.jl" begin

    # fastpow(::Number, ::Integer)
    @test fastpow(2.0, 3) == 8.0

    # fastpow(::Number, ::AbstractFloat)
    @test fastpow(2.0, 3.0) == 8.0                     # integer-valued float -> x^Int(n)
    @test fastpow(2.0, 2.5) ≈ exp(log(2.0) * 2.5)      # positive base -> exp(log(x)*n)
    @test fastpow(0.0, 2.5) == 0.0                     # non-positive base -> x^n fallback

    # fastpow(::Quantity, ::AbstractFloat)
    @test fastpow((2.0)u"m", 2.0) == (4.0)u"m^2"       # integer-valued float -> x^Int(n)
    @test fastpow((4.0)u"m^2", 0.5) == (2.0)u"m"       # non-integer -> x^n fallback

    # pow_check (used by the @pow macro)
    @test pow_check(2.0, 3) == 8.0
    @test pow_check(1.0, 5) == 1.0                     # isone(x)
    @test pow_check(2.0, 1) == 2.0                     # isone(n)
    @test pow_check(2.0, 0) == 1.0                     # iszero(n)

    # find_ind: allocation-free index search in a tuple
    @test find_ind((3, 1, 4, 1, 5), 4) == 3
    @test find_ind((3, 1, 4, 1, 5), 5) == 5
    @test find_ind((3, 1, 4, 1, 5), 9) == 0            # not found -> 0

    # find_max_tuple
    @test find_max_tuple((1, 5, 2, 4)) == 5
    @test find_max_tuple((7,)) == 7

    # max_length_tuple
    @test max_length_tuple(((1, 2), (3,), (4, 5, 6))) == 3

    # ntuple_idx: broadcast getindex over a NamedTuple of arrays
    nt = (a = [10, 20, 30], b = [40, 50, 60])
    @test ntuple_idx(nt, 2) == (a = 20, b = 50)

    # nreduce: sum f(v[i]) over a tuple without allocation
    v = (1.0, 2.0, 3.0)
    @test nreduce(x -> x^2, v) ≈ 14.0
    @test nreduce(identity, (4.0, 6.0)) ≈ 10.0

    # value_and_partial: returns f(x) and f'(x) via ForwardDiff
    val, dval = value_and_partial(x -> x^3, 2.0)
    @test val ≈ 8.0
    @test dval ≈ 12.0

    # make_tuple: wraps scalars, passes tuples through
    @test make_tuple(5) == (5,)
    @test make_tuple((1, 2)) == (1, 2)

    # ptr2string: String input -> interned string
    s = ptr2string("hello")
    @test s == "hello"

    # ptr2string: unsupported type throws ArgumentError
    @test_throws ArgumentError ptr2string(42)

    # @print macro: println returns nothing; false branch returns nothing
    @test GeoParams.@print(true, "x") === nothing
    @test GeoParams.@print(false, "x") === nothing
end
