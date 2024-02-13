export MeltViscosity,
    compute_εII!,
    compute_εII,
    compute_τII!,
    compute_τII

#=----Melt Viscosity correction ---

From

@article{schmeling2012effective,
  title={Effective shear and bulk viscosity of partially molten rock based on elastic moduli theory of a fluid filled poroelastic medium},
  author={Schmeling, Harro and Kruse, Jan Philipp and Richard, Guillaume},
  journal={Geophysical Journal International},
  volume={190},
  number={3},
  pages={1571--1578},
  year={2012},
  publisher={Blackwell Publishing Ltd Oxford, UK}
}

Defines effective shear viscosity and effective bulk viscosity as a function of melt fraction and geometry of the pores following parameterization of *Schmeling, A., Kruse, J. P., Richard, G., Effective shear and bulk viscosity of partially molten rock based on elastic moduli theory of a fluid filled poroelastic medium, Geophysical Journal International, Volume 190, Issue 3, September 2012, Pages 1571–1578, https://doi.org/10.1111/j.1365-246X.2012.05596.x*

α is the aspect ratio (between 1 and 0)
a1 is a fitting constant
a2 is a fitting constant
b1 is a fitting constant
b2 is a fitting constant
k2 is a fitting constant
k3 is a fitting constant
c3 is a fitting constant

=#

struct MeltMorphology{T,N,U1} <: AbstractCreepLaw{T}
    Name::NTuple{N,Char}
    α::GeoUnit{T,U1} # aspect ratio
    a1::GeoUnit{T,U1} # fitting constant
    a2::GeoUnit{T,U1} # fitting constant
    b1::GeoUnit{T,U1} # fitting constant
    b2::GeoUnit{T,U1} # fitting constant
    k2::GeoUnit{T,U1} # fitting constant
    k3::GeoUnit{T,U1} # fitting constant
    c3::GeoUnit{T,U1} # fitting constant

    function MeltMorphology(;
        Name = "",
        α    = 1.0NoUnits,
        a1   = 0.97NoUnits,
        a2   = 0.8NoUnits,
        b1   = 2.2455NoUnits,
        b2   = 3.45NoUnits,
        k2   = 1.25NoUnits,
        k3   = 1.29NoUnits,
        c3   = 2.4NoUnits,
    )

        # Corrections from lab experiments
        αU  = convert(GeoUnit, α)
        a1U = convert(GeoUnit, a1)
        a2U = convert(GeoUnit, a2)
        b1U = convert(GeoUnit, b1)
        b2U = convert(GeoUnit, b2)
        k2U = convert(GeoUnit, k2)
        k3U = convert(GeoUnit, k3)
        c3U = convert(GeoUnit, c3)
        # Extract struct types
        T = typeof(αU).types[1]
        U1 = typeof(αU).types[2]
        # Create struct
        N = length(Name)
        name = ntuple(i -> Name[i], Val(N))
        return new{T,N,U1}(
            name, αU, a1U, a2U, b1U, b2U, k2U, k3U, c3U
        )
    end

    function MeltMorphology(Name, α, a1, a2, b1, b2, k2, k3, c3)
        return MeltMorphology(;
            Name=Name, α=α, a1=a1, a2=a2, b1=b1, b2=b2, k2=k2, k3=k3, c3=c3
        )
    end
end

@inline function corrected_melt_viscosity(melt::MeltViscosity, η, ϕ) 
    @unpack_val α, a1, a2, b1, b2, k2, k3, c3 = melt

    k1 = a1 * (a2 + α * (1 - a2))
    c1 = (b1 * α) / (1 + b2 * fastpow(α, k3))

    η_melt = (ϕ < c1) * η * (1 - ϕ / c1)^k1

    return η_melt
end
