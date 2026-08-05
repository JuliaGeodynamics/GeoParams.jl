using UnPack
export @unpack_val, @unpack_units

"""
    This unpacks the numerical value of `GeoUnit` in a structure, without the units.
    All requested variables must be GeoUnits.

    This is a modification of the `@unpack` macro as implemented in the `UnPack.jl` package, which can be used to retrieve the full variables.

# Example
```julia
julia> struct Density{T} 
        ρ::GeoUnit{T} 
        α::GeoUnit{T} 
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_val ρ,α = r
julia> α
4.0e-5
julia> typeof(α)
Float64
```    
"""
macro unpack_val(args)
    args.head != :(=) && error("Expression needs to be of form `a, b = c`")
    items, suitecase = args.args
    items = isa(items, Symbol) ? [items] : items.args
    suitecase_instance = gensym()

    # This extracts the value, but not the units
    kd = [
        :($key = $UnPack.unpack($suitecase_instance, Val{$(Expr(:quote, key))}()).val) for
            key in items
    ]

    kdblock = Expr(:block, kd...)

    expr = quote
        local $suitecase_instance = $suitecase # handles if suitecase is not a variable but an expression
        $kdblock
        $suitecase_instance # return RHS of `=` as standard in Julia
    end
    return esc(expr)
end

"""
    This unpacks the numerical value with units of `GeoUnit` parameters in a structure
    All requested variables must be GeoUnits.

    This is a modification of the `@unpack` macro as implemented in the `UnPack.jl` package, which can be used to retrieve the full variables.

# Example
```julia
julia> struct Density{T} 
        ρ::GeoUnit{T} 
        α::GeoUnit{T} 
       end
julia> r = Density(GeoUnit(100kg/m^3),GeoUnit(4e-5/K));
julia> @unpack_units ρ,α = r
julia> α
4.0e-5 K⁻¹·⁰
julia> typeof(α)
Quantity{Float64, 𝚯⁻¹·⁰, Unitful.FreeUnits{(K⁻¹·⁰,), 𝚯⁻¹·⁰, nothing}}
```    
"""
macro unpack_units(args)
    args.head != :(=) && error("Expression needs to be of form `a, b = c`")
    items, suitecase = args.args
    items = isa(items, Symbol) ? [items] : items.args
    suitecase_instance = gensym()

    # This extracts the value with units
    kd = [
        :(
                $key =
                $UnPack.unpack($suitecase_instance, Val{$(Expr(:quote, key))}()).val .*
                $UnPack.unpack($suitecase_instance, Val{$(Expr(:quote, key))}()).unit
            ) for key in items
    ]

    kdblock = Expr(:block, kd...)

    expr = quote
        local $suitecase_instance = $suitecase # handles if suitecase is not a variable but an expression
        $kdblock
        $suitecase_instance # return RHS of `=` as standard in Julia
    end
    return esc(expr)
end
