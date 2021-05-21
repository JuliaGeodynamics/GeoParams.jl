

@with_kw struct GeoUnits{R<:Float64} @deftype R
    # Seletable input parameters
    type::String    =  "geo"       # "geo", "SI" or "none"
    temperature     =  900         # temperature
    length          =  1000        # length [km if geo, or m if SI]
    stress          =  10          # stress [MPa if geo, or Pa if SI]
    viscosity       =  1e20        # viscosity [Pas]

    # All derived parameters from the above 



    # primary characteristic units
	time;
	time_si;                 # time in SI units for material parameter scaling
	length
	#PetscScalar length_si;   # length in SI units for material parameter scaling
	#PetscScalar area_si;     # area   in SI units for material parameter scaling
	#PetscScalar volume_si;   # volume in SI units for material parameter scaling
	#PetscScalar temperature; # Kelvin (if dimensional)
	#PetscScalar force;       # additional variable for quasi-static case
	#PetscScalar angle;       # radian expressed in degrees (if dimensional)

end


