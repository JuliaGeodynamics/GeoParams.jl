struct Dual{T}
    val::T
    partial::T
end
Dual(x::T) where T = Dual(x, one(T))

import Base: +, -, *, /, ^, exp, log

for op in (:+, :-)
    @eval begin
        ($op)(a::Dual, b::Dual) = Dual(($op)(a.val, b.val), ($op)(a.partial, b.partial))
        ($op)(a::Dual, b::Number) = Dual( ($op)(a.val, b), a.partial)
        ($op)(a::Number, b::Dual) = Dual( ($op)(b.val, a), b.partial)
    end
end

(*)(a::Dual, b::Dual)   = Dual(a.val*b.val, a.val*b.partial + b.val*a.partial)
(*)(a::Dual, b::Number) = Dual(a.val*b, b*a.partial)
(*)(a::Number, b::Dual) = Dual(a*b.val, a*b.partial)

(^)(a::Dual, n::Integer) = Dual(a.val^n, n*a.partial)
(^)(a::Dual, n::Float64) = Dual(a.val^n, n*a.partial)

(/)(a::Dual, b::Dual)   = Dual(a.val/b.val, (b.val*a.partial - a.val*b.partial) / (b.val^2))
(/)(a::Dual, b::Number) = Dual(a.val/b, b*a.partial / (b^2))
(/)(a::Number, b::Dual) = Dual(a/b.val, - a.val*b.partial / (b.val^2))

function exp(a::Dual) 
    expa = exp(a.val)
    Dual(expa, expa*a.partial)
end

log(a::Dual) = Dual( log(a.val), a.partial / a.val)
log10(a::Dual) = Dual( log10(a.val), a.partial * 0.43429448190325176 / a.val)

sin(a::Dual) = Dual( sin(a.val), a.partial * cos(a.val))
sind(a::Dual) = Dual( sind(a.val), a.partial * cosd(a.val))
cos(a::Dual)  = Dual( cos(a.val), -a.partial * sin(a.val))
cosd(a::Dual) = Dual( cosd(a.val), -a.partial * sind(a.val))

n = 1.1
ForwardDiff.derivative(x->GeoParams.fastpow(x, n), 1.0)


ff(x,n) = exp(log( Dual(x))*n)
gg(x,n) = x^n, exp(log( x)*n)

@btime ff($2.0,$n)
@btime gg($2.0,$n)

##############################
a = Dual(2)

f1(x) = x + x
f2(x) = 2*x 

df(fn::F, x::Number) where F = fn(Dual(x))

df(f1, 2)
df(f2, 2)

a^1.1
Dual(2)/3

compute_εII(
    a::DiffusionCreep, TauII::_T; T, P=0.0, f=1.0, d=1.0, kwargs...
) 


# Define a linear viscous creep law ---------------------------------
x1 = DiffusionCreep()
TauII = 1e6
args = (T=1e3, P=1e9)

function foo1(x1, TauII, args)
    f = compute_εII(x1, TauII, args)
    dfdTau = dεII_dτII(x1, TauII, args)
    return f, dfdTau
end

function foo2(x1, TauII, args)
    f = compute_εII(x1, TauII, args)
    dfdTau = ForwardDiff.derivative(x-> compute_εII(x1, x, args), TauII)
    return f, dfdTau
end

function foo3(x1, TauII, args)
    DualTau = compute_εII(x1, Dual(TauII), args)

    f, dfdTau = DualTau.val, DualTau.partial
    return f, dfdTau
end


foo1(x1, TauII, args) .≈ foo2(x1, TauII, args) .≈ foo3(x1, TauII, args)

@btime foo($x1, $TauII, $args)
@btime foo2($x1, $TauII, $args)
@btime foo3($x1, $TauII, $args)

a = Dual(TauII)
compute_εII(x1, TauII, args)
compute_εII(x1, a, args)

dεII_dτII(x1, TauII, args)

@btime foo($x1, $TauII, $args)
@btime foo2($x1, $TauII, $args)
@btime compute_εII($x1, $a, $args)
@btime ForwardDiff.derivative(x-> compute_εII($x1, x, $args), $TauII)


#########################

f=1.0
d=1.0
@unpack_val n, r, p, A, E, V, R = x1
FT, FE = x1.FT, x1.FE

A * fastpow(TauII * FT, n) *
    fastpow(f, r) *
    fastpow(d, p) *
    exp(-(E + P * V) / (R * T)) / FE

A * GeoParams.fastpow(a * FT, n) *
    GeoParams.fastpow(f, r) *
    GeoParams.fastpow(d, p) *
    exp(-(E + P * V) / (R * T)) / FE

    A * n * FT * GeoParams.fastpow(FT * TauII, -1 + n) *
           GeoParams.fastpow(f, r) *
           GeoParams.fastpow(d, p) *
           exp((-E - P * V) / (R * T)) *
           (1 / FE)