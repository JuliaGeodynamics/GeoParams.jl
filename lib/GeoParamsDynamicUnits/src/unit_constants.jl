struct DynamicUnitSymbol{S} end

const _SymbolicQuantity = typeof(us"Pa")
const _MPa_UNIT = Ref{_SymbolicQuantity}(us"Pa")
const _GPa_UNIT = Ref{_SymbolicQuantity}(us"Pa")
const _KBAR_UNIT = Ref{_SymbolicQuantity}(us"bar")
const _MICROWATT_UNIT = Ref{_SymbolicQuantity}(us"W")

@inline _materialize(::DynamicUnitSymbol{:MPa}) = _MPa_UNIT[]
@inline _materialize(::DynamicUnitSymbol{:GPa}) = _GPa_UNIT[]
@inline _materialize(::DynamicUnitSymbol{:kbar}) = _KBAR_UNIT[]
@inline _materialize(::DynamicUnitSymbol{:μW}) = _MICROWATT_UNIT[]

Base.show(io::IO, ::DynamicUnitSymbol{S}) where {S} = print(io, S)
Base.:*(value::Number, unit::DynamicUnitSymbol) = value * _materialize(unit)
Base.:*(unit::DynamicUnitSymbol, value::Number) = _materialize(unit) * value
Base.:/(value::Number, unit::DynamicUnitSymbol) = value / _materialize(unit)
Base.:/(unit::DynamicUnitSymbol, value::Number) = _materialize(unit) / value
Base.:^(unit::DynamicUnitSymbol, power::Number) = _materialize(unit)^power
Base.:(==)(unit::DynamicUnitSymbol, value::DynamicQuantities.AbstractQuantity) =
    _materialize(unit) == value
Base.:(==)(value::DynamicQuantities.AbstractQuantity, unit::DynamicUnitSymbol) = unit == value

# Symbolic quantities retain the unit scale supplied by the user.
const km = us"km"
const m = us"m"
const cm = us"cm"
const mm = us"mm"
const μm = us"μm"
const Myr = us"Myr"
const yr = us"yr"
const s = us"s"
const kg = us"kg"
const g = us"g"
const Pa = us"Pa"
const MPa = DynamicUnitSymbol{:MPa}()
const GPa = DynamicUnitSymbol{:GPa}()
const bar = us"bar"
const kbar = DynamicUnitSymbol{:kbar}()
const Pas = us"Pa*s"
const K = us"K"
const C = ua"degC"
const mol = us"mol"
const kJ = us"kJ"
const J = us"J"
const Watt = us"W"
const μW = DynamicUnitSymbol{:μW}()

# DynamicQuantities represents dimensionless values as ordinary numbers.
const NoUnits = 1.0
