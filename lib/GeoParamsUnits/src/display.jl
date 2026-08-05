# Define a view for the GEO_Units structure
function show(io::IO, g::GeoUnits{T}) where {T}
    return print(
        io,
        "Employing $T units \n",
        "Characteristic values: \n",
        "         length:      $(g.length)\n",
        "         time:        $(round(ustrip(g.time), digits = 4)) $(unit(g.time))\n",
        "         stress:      $(g.stress)\n",
        "         temperature: $(Float64(g.temperature))\n",
    )
end

# Format floating-point exponents without extending Unitful's internal display methods.
function superscript(i::Rational{Int64}; io = nothing)
    if io === nothing
        iocontext_value = nothing
    else
        iocontext_value = get(io, :fancy_exponent, nothing)
    end
    if iocontext_value isa Bool
        fancy_exponent = iocontext_value
    else
        v = get(ENV, "UNITFUL_FANCY_EXPONENTS", Sys.isapple() ? "true" : "false")
        t = tryparse(Bool, lowercase(v))
        fancy_exponent = (t === nothing) ? false : t
    end
    if fancy_exponent
        return superscript(float(i))
    else
        return "^" * string(float(i))
    end
end

function superscript(i::Float64)
    return map(repr(i)) do c
        if c == '-'
            '\u207b'
        elseif c == '1'
            '\u00b9'
        elseif c == '2'
            '\u00b2'
        elseif c == '3'
            '\u00b3'
        elseif c == '4'
            '\u2074'
        elseif c == '5'
            '\u2075'
        elseif c == '6'
            '\u2076'
        elseif c == '7'
            '\u2077'
        elseif c == '8'
            '\u2078'
        elseif c == '9'
            '\u2079'
        elseif c == '0'
            '\u2070'
        elseif c == '.'
            '\u0387'
        else
            error("unexpected character")
        end
    end
end
