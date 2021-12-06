# This contains routines to read phase diagrams from disk 
using DelimitedFiles, Interpolations
import Base.show

export PhaseDiagram_LookupTable, Read_LaMEM_Perple_X_Diagram

"""
    Contains data of a Phase Diagram that is regularly spaced in P & T

"""
struct PhaseDiagram_LookupTable
    Type        ::  String
    HeaderText  ::  Any
    Name        ::  String
    meltRho     ::  Interpolations.Extrapolation
    meltFrac    ::  Interpolations.Extrapolation
    rockRho     ::  Interpolations.Extrapolation
    Rho         ::  Interpolations.Extrapolation
    rockVp      ::  Any
    rockVs      ::  Any
    rockVpVs    ::  Any
end


"""
    PD_Data = Read_LaMEM_Perple_X_Diagram(fname::String)


Reads a precomputed phase diagram in the LaMEM/Perple\\_X format (which is a phase diagram computed using Perple\\_X, but formatted in a manner that is readable using LaMEM).
The data is stored in the `PhaseDiagram_LookupTable` structure.

"""
function Read_LaMEM_Perple_X_Diagram(fname::String)
    
    # Read header: 
    #  the first 50 lines are comments (freely useable), followed by data 

    n           = 55
    header      = open(readlines, `head -n $(n) $(fname)`)
    header_text = header[1:50]

    T0      =   parse(Float64,header[50])       # in K
    dT      =   parse(Float64,header[51])
    numT    =   parse(Int64,  header[52])

    P0      =   parse(Float64,header[53])*1e5   # in Pa (file gives it in bar)
    dP      =   parse(Float64,header[54])*1e5
    numP    =   parse(Int64,  header[55])
    
    Tvec    =   T0:dT:(T0+dT*(numT-1))          # 1D vector
    Pvec    =   P0:dP:(P0+dP*(numP-1))
    
    # In the LaMEM/Perple_X file format, the first 50 lines are comments
    data = readdlm(fname,skipstart=55,header=false);        # read numerical data

    # Shape to 2D:
    siz      =   (numP, numT)
    
    # Data that ishould always be present
    meltRho         =       reshape(data[:,1], siz);    # in kg/m3
    meltFrac        =       reshape(data[:,2], siz);    # in wt-%
    rockRho         =       reshape(data[:,3], siz);    # in kg/m3
    totalRho        =       rockRho.*(1.0 .- meltFrac) + meltRho.*meltFrac;

    T_K             =       reshape(data[:,4], siz);    # in K  
    P_Pa            =   1e5*reshape(data[:,5], siz);    # in Pa

    # Create interpolation objects
    intp_meltRho    =   LinearInterpolation((Tvec, Pvec), meltRho,      extrapolation_bc = Line()); 
    intp_meltFrac   =   LinearInterpolation((Tvec, Pvec), meltFrac,     extrapolation_bc = Line()); 
    intp_rockRho    =   LinearInterpolation((Tvec, Pvec), rockRho,      extrapolation_bc = Line()); 
    intp_totalRho   =   LinearInterpolation((Tvec, Pvec), totalRho,     extrapolation_bc = Line()); 
    
    # optional data
    if size(data,2)>=6
        rockVp          =   reshape(data[:,6], siz);    # in km/s
        intp_rockVp     =   LinearInterpolation((Tvec, Pvec), rockVp,  extrapolation_bc = Line()); 
    else
        intp_rockVp     =   nothing
    end
    if size(data,2)>=7
        rockVs          =   reshape(data[:,7], siz);    # in km/s
        intp_rockVs     =   LinearInterpolation((Tvec, Pvec), rockVs,  extrapolation_bc = Line()); 
    else
        intp_rockVs     =   nothing
    end
    if size(data,2)>=8
        rockVpVs        =   reshape(data[:,8], siz);    # ratio []
        intp_rockVpVs   =   LinearInterpolation((Tvec, Pvec), rockVpVs,  extrapolation_bc = Line()); 
    else
        intp_rockVpVs   =   nothing
    end

    # Store in phase diagram structure
    PD_data = PhaseDiagram_LookupTable("PerpleX_LaMEM", header_text, fname, 
                                intp_meltRho,   intp_meltFrac,  intp_rockRho,   intp_totalRho,
                                intp_rockVp,    intp_rockVp,    intp_rockVpVs)

    return PD_data
end


# Print info 
function show(io::IO, d::PhaseDiagram_LookupTable)  
    T = d.rockRho.itp.ranges[1]
    P = d.rockRho.itp.ranges[2]

    println(io, "$(d.Type) Phase Diagram Lookup Table: ")  
    println(io, "   File    :   $(d.Name)")  
    println(io, "   T       :   $(minimum(T)) - $(maximum(T))")  
    println(io, "   P       :   $(minimum(P)) - $(maximum(P))")  
    
   
    lst = fieldnames(typeof(d))
    for i=4:length(lst)
        if !isnothing(getfield(d,lst[i]))
            if i==4
                println(io, "   fields  :   $(lst[i])")  
            else
                println(io, "               $(lst[i])")  
            end
        end
    end
    
end


# Calculation routine
function ComputeDensity(P,T, s::PhaseDiagram_LookupTable)
   
    ρ = s.Rho.(T,P)

    return ρ
end