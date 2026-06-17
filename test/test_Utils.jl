using Test
using GeoParams, Unitful
import GeoParams: fastpow, pow_check

@testset "Utils.jl" begin

    # fastpow(::Number, ::Integer)
    @test fastpow(2.0, 3) == 8.0

    # fastpow(::Number, ::AbstractFloat)
    @test fastpow(2.0, 3.0) == 8.0                     # integer-valued float → x^Int(n)
    @test fastpow(2.0, 2.5) ≈ exp(log(2.0) * 2.5)      # positive base → exp(log(x)*n)
    @test fastpow(0.0, 2.5) == 0.0                     # non-positive base → x^n fallback

    # fastpow(::Quantity, ::AbstractFloat)
    @test fastpow((2.0)u"m", 2.0) == (4.0)u"m^2"       # integer-valued float → x^Int(n)
    @test fastpow((4.0)u"m^2", 0.5) == (2.0)u"m"       # non-integer → x^n fallback

    # pow_check (used by the @pow macro)
    @test pow_check(2.0, 3) == 8.0
    @test pow_check(1.0, 5) == 1.0                     # isone(x)
    @test pow_check(2.0, 1) == 2.0                     # isone(n)
    @test pow_check(2.0, 0) == 1.0                     # iszero(n)
end
