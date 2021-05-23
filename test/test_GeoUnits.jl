# Tests the GeoUnits
using Test
using GeoParams


CharUnits   =   GEO_units(viscosity=1e19, length=1000km);
@test CharUnits.length==1000km
@test CharUnits.Pa==10000000Pa

CharUnits   =   SI_units();
@test CharUnits.length==1000m
@test CharUnits.Pa==10Pa

CharUnits   =   NO_units();
@test CharUnits.length==1000
@test CharUnits.Pa==10
