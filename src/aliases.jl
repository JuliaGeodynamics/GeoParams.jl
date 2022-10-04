export @use

const VALID_KWARGS = (
    :density,
    :gravity,
    :heatcapacity,
    :conductivity,
    :latent_heat,
    :radioactive_heat,
    :shearheating,
    :pwave_velocity,
    :swave_velocity,
    :meltfraction,
)

macro use(args...)
    checkargs(args...)
    esc(create_module(__module__, args...))
    return esc(use_module(__module__, args[1]))
end

function checkargs(args...)
    if (length(args) < 2)
        error("arguments missing.")
    end # NOTE: this needs to be adapted to: 1 + nb math symbols needed
    if (args[1] != :GeoParamsAliases)
        error("the first argument must be `GeoParamsAliases`")
    end
    if !all(is_kwarg.(args[2:end]))
        error(
            "all arguments starting from the second must be keyword arguments (see macro documentation).",
        )
    end
end

is_kwarg(arg) = isa(arg, Expr) && (arg.head == :(=))
split_kwargs(kwargs) = Dict(x.args[1] => x.args[2] for x in kwargs)

function validate_kwargkeys(kwargs::Dict, macroname::String)
    for k in keys(kwargs)
        if !(k in VALID_KWARGS)
            error(
                "Invalid keyword argument in $macroname call: `$k`. Valid keyword arguments are: `$(join(valid_kwargs, "`, `"))`.",
            )
        end
    end
end

function create_module(caller::Module, modulename::Symbol, kwargs_expr::Expr...)
    if isdefined(caller, modulename)
        error("Module $modulename already exists in $caller.")
    end
    kwargs = split_kwargs(kwargs_expr)
    validate_kwargkeys(kwargs, "@use GeoParamsAliases")
    exports = Symbol[]
    aliases = Expr[]
    for key in keys(kwargs)
        param = kwargs[key]
        param! = Symbol(param, "!")
        fname = Symbol("compute_" * String(key))
        fname! = Symbol("compute_" * String(key) * "!")

        # Exports
        push!(exports, param)
        # Aliases
        push!(
            aliases,
            quote
                $param(args...) = GeoParams.$fname(args...)
            end,
        )
        # Only in case compute_...! method exists set alias for it as well
        if isdefined(GeoParams, fname!)
            push!(exports, param!)
            push!(
                aliases,
                quote
                    $param!(args...) = GeoParams.$fname!(args...)
                end,
            )
        end
    end

    module_expr = :(module $modulename # NOTE: there cannot be any newline before 'module $modulename' or it will create a begin end block and the module creation will fail.
    using GeoParams: GeoParams
    export $(exports...)
    $(aliases...)
    end)
    @eval(caller, $module_expr)
end

function use_module(caller::Module, modulename::Symbol)
    @eval(caller, using $(Symbol(caller)).$modulename)
end
