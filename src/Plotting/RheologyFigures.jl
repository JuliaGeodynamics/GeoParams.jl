# Rheology figures of chapter 5 of Becker & Faccenna, Tectonic Geodynamics (draft, 2026).
# All keyword defaults are the book's choices; any creep law of GeoParams can be passed instead.
# Units are SI throughout (Pa, K, m, s); stresses on the axes are differential stresses σd = √3 τII.

import GeoParams: PlotFlowStressTemperature, PlotWaterFugacity, PlotDeformationMechanismMaps,
    PlotOceanicStrength, PlotContinentalStrength
using SpecialFunctions: erf

const _yr = 365.25 * 24 * 3600
const _g = 9.81

# differential stress σd = √3 τII [Pa] at strain rate ε̇II [1/s], and the inverse
_σd(law, ε̇; args...) = sqrt(3) * compute_τII(law, ε̇; args...)
_ε̇II(law, σd; args...) = compute_εII(law, σd / sqrt(3); args...)

# water variable: a number (C_OH or a fixed fugacity) or a function (P, T) -> f
_water(f::Function, P, T) = f(P, T)
_water(f, P, T) = f

# a law, or a function T -> law (e.g. the two temperature regimes of GBS)
_law(l::Function, T) = l(T)
_law(l, T) = l

_halfspace_T(z, age; Ts = 273.15, Tm = 1673.15, κ = 1.0e-6) = Ts + (Tm - Ts) * erf(z / (2 * sqrt(κ * age * 1.0e6 * _yr)))
_linear_T(z, q; Ts = 273.15, k = 3.0, Tmax = 1623.15) = min(Ts + q * z / k, Tmax)
_lithostatic_p(z; z1 = 0.0, ρ1 = 2750.0, ρ2 = 3300.0) = _g * (ρ1 * min(z, z1) + ρ2 * max(z - z1, 0))

# Byerlee (1978) frictional strength |σd| [Pa] for pore pressure factor λ (eqs. 5.30-5.31)
function _byerlee(p; λ = 0.0, mode = :compression)
    p′ = (1 - λ) * p
    s = mode === :compression ? -1 :
        mode === :extension ? 1 :
        throw(ArgumentError("mode must be :compression or :extension, got $mode"))
    σ(μ, C) = 2 * (μ * p′ + C) / (sqrt(1 + μ^2) + s * μ)
    return min(σ(0.85, 0.0), σ(0.6, 50.0e6))
end

_book_fig5_6_laws() = [
    (label = "qu-J84", law = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Jaoul_1984), color = :blue, light_above = 650),
    (label = "qu-RB04", law = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Rutter_2004), color = :orange, light_above = 650),
    (label = "qu-T19", law = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Tokle_2019), color = :red, light_above = 650),
    (label = "qu-L21", law = SetDislocationCreep(GeoParams.Dislocation.lowP_wet_quartzite_Lusk_2021), color = :lime, light_above = 650),
    (label = "fe-RD06", law = SetDislocationCreep(GeoParams.Dislocation.wet_anorthite_Rybacki_2006), color = :cyan),
    (label = "ol-HK04", law = SetDislocationCreep(GeoParams.Dislocation.wet_olivine_Hirth_2003), color = :magenta),
]

function PlotFlowStressTemperature(;
        laws = _book_fig5_6_laws(),
        P = 0.5e9,
        εII = (1.0e-15, 1.0e-13),
        T_C = range(200, 1200; length = 201),
        f = compute_water_fugacity,
        size = (600, 600),
    )
    fig = Figure(; size)
    ax = Axis(
        fig[1, 1]; xscale = log10, yreversed = true, xaxisposition = :top,
        xticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000], yticks = first(T_C):200:last(T_C),
        limits = (1, 1000, first(T_C) - 50, last(T_C)),
        xlabel = "differential stress σ₁-σ₃ [MPa]", ylabel = "temperature [°C]"
    )
    T = T_C .+ 273.15
    for l in laws
        fl = get(l, :f, f)
        σ(ε̇) = [_σd(l.law, ε̇; T = T[i], P, f = _water(fl, P, T[i])) for i in eachindex(T)] ./ 1.0e6
        σ_lo, σ_hi = σ(first(εII)), σ(last(εII))
        Tl = get(l, :light_above, Inf)
        solid, light = T_C .<= Tl, T_C .>= Tl
        band!(ax, Point2f.(σ_lo[solid], T_C[solid]), Point2f.(σ_hi[solid], T_C[solid]); color = (l.color, 0.6), label = l.label)
        any(light) && band!(ax, Point2f.(σ_lo[light], T_C[light]), Point2f.(σ_hi[light], T_C[light]); color = (l.color, 0.15))
    end
    text!(ax, 1.3, first(T_C) + 30; text = "p = $(P / 1.0e9) GPa", fontsize = 20)
    axislegend(ax; position = :rb)
    return fig
end

function PlotWaterFugacity(;
        P = range(0.1e9, 2.0e9; length = 150),
        T = range(375, 1775; length = 150),
        heatflow = (0.05, 0.1),
        crustal_depth = (20.0e3, 40.0e3),
        k = 3.0,
        ρ_crust = 2750.0,
        ages = (10, 25, 50, 150),
        Tm = 1623.15,
        ρ_mantle = 3000.0,
        size = (1000, 520),
    )
    fw = [compute_water_fugacity(p, t) for p in P, t in T]
    fw_S = [compute_water_fugacity_Shinevar(p, t) for p in P, t in T]
    C_OH = compute_water_content_olivine.(fw, P, T')
    pMPa = P ./ 1.0e6

    fig = Figure(; size)
    ax1 = Axis(fig[1, 1]; xaxisposition = :top, xlabel = "p [MPa]", ylabel = "T [K]")
    cf = contourf!(ax1, pMPa, T, fw ./ 1.0e6; levels = 0:250:5000, colormap = :roma, extendhigh = :auto)
    contour!(ax1, pMPa, T, fw_S ./ 1.0e6; levels = 0:250:3000, color = :black, linewidth = 0.5)
    contour!(ax1, pMPa, T, fw_S ./ 1.0e6; levels = 1000:1000:3000, color = :black, labels = true, labelsize = 14)
    z = range(crustal_depth...; length = 50)
    for q in heatflow
        lines!(ax1, _lithostatic_p.(z; ρ1 = ρ_crust, ρ2 = ρ_crust) ./ 1.0e6, _linear_T.(z, q; k); color = :lightgray, linewidth = 3)
    end
    Colorbar(fig[2, 1], cf; vertical = false, label = "f_w [MPa]", width = Relative(0.6))

    ax2 = Axis(fig[1, 2]; xaxisposition = :top, yaxisposition = :right, xlabel = "p [MPa]", ylabel = "T [K]")
    cf2 = contourf!(ax2, pMPa, T, C_OH; levels = 0:100:2500, colormap = Reverse(:vik), extendhigh = :auto)
    z = range(0, 75.0e3; length = 200)
    for age in ages
        lines!(ax2, _lithostatic_p.(z; ρ1 = ρ_mantle, ρ2 = ρ_mantle) ./ 1.0e6, _halfspace_T.(z, age; Tm); color = :lightgray, linewidth = 2)
    end
    Colorbar(fig[2, 2], cf2; vertical = false, label = "C_OH [H/10⁶ Si]", width = Relative(0.6))
    for ax in (ax1, ax2)
        limits!(ax, extrema(pMPa)..., extrema(T)...)
    end
    return fig
end

function _book_fig5_10_mechanisms()
    gbs_cold = SetGrainBoundarySliding(GeoParams.GBS.cold_dry_olivine_Hirth_2003)
    gbs_hot = SetGrainBoundarySliding(GeoParams.GBS.hot_dry_olivine_Hirth_2003)
    return [
        (name = "diffusion", law = SetDiffusionCreep(GeoParams.Diffusion.dry_olivine_Hirth_2003; V = 6.0e-6m^3 / mol), color = RGBf(0.85, 0.85, 0.83)),
        (name = "dislocation", law = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003; V = 18.0e-6m^3 / mol), color = RGBf(0.7, 0.63, 0.48)),
        (name = "GBS", law = T -> T < 1523.15 ? gbs_cold : gbs_hot, color = RGBf(0.43, 0.14, 0.14)),
        (name = "Peierls", law = SetNonLinearPeierlsCreep(GeoParams.NonLinearPeierls.dry_olivine_Mei_2010), color = RGBf(0.63, 0.42, 0.35), Tmax = 1523.15),
    ]
end

function PlotDeformationMechanismMaps(;
        mechanisms = _book_fig5_10_mechanisms(),
        T = 1523.15,
        P = 1.0e9,
        d = 3.0e-3,
        age = 100,
        Tm = 1673.15,
        dTdz = 0.15e-3,
        ρ = 3300.0,
        yield_law = SetPeierlsCreep(GeoParams.Peierls.dry_olivine_Goetze_1979),
        yield_εII = 1.0e-15,
        yield_σmin = 300.0e6,
        log_σ = range(-1, 4; length = 251),
        log_d = range(-3, 2; length = 201),
        T_C = range(0, 1650; length = 166),
        depth = range(5.0e3, 600.0e3; length = 239),
        size = (650, 1500),
    )
    function rates(σ, T, P, d)
        return map(mechanisms) do m
            T <= get(m, :Tmax, Inf) ? _ε̇II(_law(m.law, T), σ; T, P, d) : 0.0
        end
    end
    function panel!(ax, x, conditions)
        dom = zeros(Int, length(x), length(log_σ))
        log_ε̇ = zeros(length(x), length(log_σ))
        for (i, xᵢ) in enumerate(x), (j, s) in enumerate(log_σ)
            ε̇s = rates(10^s * 1.0e6, conditions(xᵢ)...)
            dom[i, j] = argmax(ε̇s)          # dashpots in series: the fastest mechanism dominates (eq. 5.29)
            log_ε̇[i, j] = log10(sum(ε̇s))
        end
        heatmap!(ax, x, log_σ, dom; colormap = [m.color for m in mechanisms], colorrange = (1, length(mechanisms)))
        contour!(ax, x, log_σ, log_ε̇; levels = -19:2:-3, color = :black, linewidth = 1, labels = true)
        contour!(ax, x, log_σ, log_ε̇; levels = -20:2:-2, color = :black, linewidth = 0.4)
        return ax
    end

    fig = Figure(; size)
    ylabel = "log₁₀(σ [MPa])"
    ax_a = Axis(fig[1, 1]; xlabel = "log₁₀(d [mm])", ylabel, title = "a) T = $(round(Int, T - 273.15)) °C, p = $(P / 1.0e9) GPa")
    panel!(ax_a, log_d, x -> (T, P, 10^x * 1.0e-3))
    ax_b = Axis(fig[2, 1]; xlabel = "T [°C]", ylabel, title = "b) d = $(d * 1.0e3) mm, p = $(P / 1.0e9) GPa")
    panel!(ax_b, T_C, x -> (x + 273.15, P, d))
    ax_c = Axis(fig[3, 1]; xlabel = "depth [km]", ylabel, title = "c) d = $(d * 1.0e3) mm, $age Myr half-space cooling")
    T_z(z) = _halfspace_T(z, age; Tm) + dTdz * z
    p_z(z) = _lithostatic_p(z; ρ1 = ρ, ρ2 = ρ)
    panel!(ax_c, depth ./ 1.0e3, x -> (T_z(x * 1.0e3), p_z(x * 1.0e3), d))
    if yield_law !== nothing
        σy = [_σd(yield_law, yield_εII; T = T_z(z)) for z in depth]
        σy = [isnan(σ) || σ < yield_σmin ? yield_σmin : σ for σ in σy]
        lines!(ax_c, depth ./ 1.0e3, log10.(σy ./ 1.0e6); color = :black, linestyle = :dash, linewidth = 4)
    end
    ylims!(ax_c, extrema(log_σ)...)
    Legend(fig[1:3, 2], [PolyElement(color = m.color) for m in mechanisms], [m.name for m in mechanisms])
    return fig
end

_book_fig5_12_laws() = [
    # Table 5.2 uses C_OH^1 with C_OH = 1000 H/10⁶ Si; this law has r = 1.2, so f = 1000^(1/1.2) gives the same A f^r
    (label = "wet", law = SetDislocationCreep(GeoParams.Dislocation.wet_olivine1_Hirth_2003), f = 1000^(1 / 1.2), linestyle = :solid),
    (label = "dry", law = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003), f = 1.0, linestyle = :dash),
]

function PlotOceanicStrength(;
        laws = _book_fig5_12_laws(),
        ages = (25, 50, 150),
        colors = (:orange, :red, :darkred),
        εII = 5.0e-15,
        Tm = 1673.15,
        h_crust = 7.0e3,
        ρ_crust = 2750.0,
        ρ_mantle = 3300.0,
        λ = 0.0,
        depth = range(0, 100.0e3; length = 1001),
        size = (600, 700),
    )
    p = _lithostatic_p.(depth; z1 = h_crust, ρ1 = ρ_crust, ρ2 = ρ_mantle)
    σ_ext = _byerlee.(p; λ, mode = :extension)
    σ_comp = _byerlee.(p; λ, mode = :compression)
    zkm = depth ./ 1.0e3

    fig = Figure(; size)
    ax = Axis(
        fig[1, 1]; xaxisposition = :top, yreversed = true, limits = (-20, 1600, 0, last(zkm)),
        xlabel = "differential stress σ₁-σ₃ [MPa]", ylabel = "depth [km]"
    )
    axT = Axis(fig[1, 1]; yreversed = true, limits = (273, Tm + 40, 0, last(zkm)), xlabel = "temperature [K]", xgridvisible = false)
    hideydecorations!(axT)
    hidespines!(axT)
    lines!(ax, σ_comp ./ 1.0e6, zkm; color = :navy, linewidth = 1)
    lines!(ax, σ_ext ./ 1.0e6, zkm; color = :navy, linewidth = 3)
    strength = Dict()        # depth-integrated strength [TN/m]
    for (age, color) in zip(ages, colors)
        T = _halfspace_T.(depth, age; Tm)
        lines!(axT, T, zkm; color, linestyle = :dot)
        for (k, l) in enumerate(laws)
            σ_creep = [_σd(l.law, εII; T = T[i], P = p[i], f = _water(get(l, :f, 1.0), p[i], T[i])) for i in eachindex(depth)]
            σ_min = min.(σ_creep, σ_ext)       # trapezoidal rule, works for any depth vector
            strength[(age, l.label)] = sum((σ_min[1:(end - 1)] .+ σ_min[2:end]) ./ 2 .* diff(depth)) / 1.0e12
            lines!(
                ax, ifelse.(σ_creep .< σ_ext, σ_creep ./ 1.0e6, NaN), zkm; color, linestyle = get(l, :linestyle, :solid),
                linewidth = 2, label = k == 1 ? "$age Myr" : nothing
            )
        end
    end
    axislegend(ax; position = :rt)
    return fig, strength
end

function _book_fig5_13_panels()
    anorthite = SetDislocationCreep(GeoParams.Dislocation.wet_anorthite_Rybacki_2006)
    quartz = SetDislocationCreep(GeoParams.Dislocation.wet_quartzite_Tokle_2019)
    return [
        (title = "a)", layers = [(bottom = 35.0e3, law = anorthite)]),
        (title = "b)", layers = [(bottom = 35.0e3 / 3, law = quartz), (bottom = 35.0e3, law = anorthite)]),
    ]
end

function PlotContinentalStrength(;
        panels = _book_fig5_13_panels(),
        mantle = SetDislocationCreep(GeoParams.Dislocation.dry_olivine_Hirth_2003; V = 18.0e-6m^3 / mol),
        f = compute_water_fugacity,
        heatflow = (0.045, 0.065, 0.085),
        colors = (:brown, :red, :orange),
        k = 3.0,
        Tmax = 1623.15,
        εII = 1.0e-15,
        λ = 0.4,
        ρ_crust = 2750.0,
        ρ_mantle = 3300.0,
        depth = range(0, 60.0e3; length = 1201),
        size = (1000, 650),
    )
    fig = Figure(; size)
    zkm = depth ./ 1.0e3
    for (col, panel) in enumerate(panels)
        h_crust = last(panel.layers).bottom
        p = _lithostatic_p.(depth; z1 = h_crust, ρ1 = ρ_crust, ρ2 = ρ_mantle)
        σ_byerlee = _byerlee.(p; λ, mode = :compression)
        ax_σ = Axis(
            fig[1, col]; xaxisposition = :top, yreversed = true, limits = (0, 2000, 0, last(zkm)),
            xlabel = "differential stress σ₁-σ₃ [MPa]", ylabel = "depth [km]", title = panel.title
        )
        ax_T = Axis(fig[1, col]; yreversed = true, limits = (273, Tmax + 90, 0, last(zkm)), xlabel = "temperature [K]", xgridvisible = false)
        hideydecorations!(ax_T)
        hidespines!(ax_T)
        lines!(ax_σ, σ_byerlee ./ 1.0e6, zkm; color = :navy, linewidth = 2)
        hlines!(ax_σ, [l.bottom for l in panel.layers] ./ 1.0e3; color = :gray, linestyle = :dash)
        for (q, color) in zip(heatflow, colors)
            T = _linear_T.(depth, q; k, Tmax)
            lines!(ax_T, T, zkm; color, linestyle = :dot)
            σ_creep = map(eachindex(depth)) do i
                i_layer = findfirst(l -> depth[i] <= l.bottom, panel.layers)
                layer = isnothing(i_layer) ? (law = mantle,) : panel.layers[i_layer]
                fl = _water(get(layer, :f, f), p[i], T[i])
                _σd(layer.law, εII; T = T[i], P = p[i], f = fl)
            end
            lines!(ax_σ, min.(σ_creep, σ_byerlee) ./ 1.0e6, zkm; color, linewidth = 2.5, label = "$(round(Int, 1000q)) mW/m²")
        end
        axislegend(ax_σ; position = :rt)
    end
    return fig
end
