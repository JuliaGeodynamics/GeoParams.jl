function Base.show(io::IO, units::GeoUnits{T}) where {T}
    temperature = if T === GEO
        "$(round(ustrip(C, units.temperature), digits = 4)) °C"
    else
        "$(round(ustrip(units.temperature), digits = 4)) $(unit(units.temperature))"
    end
    return print(
        io,
        "Employing $T units \n",
        "Characteristic values: \n",
        "         length:      $(units.length)\n",
        "         time:        $(units.time)\n",
        "         stress:      $(units.stress)\n",
        "         temperature: $temperature\n",
    )
end

function superscript(value::Rational{Int64}; io = nothing)
    fancy = if io === nothing
        nothing
    else
        get(io, :fancy_exponent, nothing)
    end
    if !(fancy isa Bool)
        fancy = something(tryparse(Bool, lowercase(get(ENV, "UNITFUL_FANCY_EXPONENTS", "false"))), false)
    end
    return fancy ? superscript(float(value)) : "^" * string(float(value))
end

function superscript(value::Float64)
    return map(repr(value)) do character
        character == '-' && return '\u207b'
        character == '1' && return '¹'
        character == '2' && return '²'
        character == '3' && return '³'
        character == '4' && return '⁴'
        character == '5' && return '⁵'
        character == '6' && return '⁶'
        character == '7' && return '⁷'
        character == '8' && return '⁸'
        character == '9' && return '⁹'
        character == '0' && return '⁰'
        character == '.' && return '·'
        error("unexpected character")
    end
end
