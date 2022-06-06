using Test
using GeoParams, LinearAlgebra, DelimitedFiles
@testset "ZirconAge.jl" begin

# declare constant variables
s2y 			= 365.0*24.0*3600.0 							# second to year
ZirconData  	= ZirconAgeData();		# default data

# Parameters for Zircon statistical analysis
n_analyses		= 300											# number of synthetic zircon analyses

# read input file (this step should be skipped, when using tracers)		
filename 		= 	"Data/Tt15_st_Bl_rad5_0.0126.txt"				# path to input file 
Tt_paths		=	readdlm(filename,',')
n_paths			= 	size(Tt_paths,2)-1							# -1 ro remove first column that corresponds to the time
time_years 		= 	Tt_paths[:,1]./s2y							# time fron seconds to years

# Here a matrix is constructed that only contains the temperature information:
Tt_paths_Temp 	= 	Tt_paths[:,2:end]	
	
# first: test case in which we provide the Tt-path as matrix 
time_years, prob, ages_eruptible, number_zircons, T_av_time, T_sd_time,zircon_cumulativePDF  = compute_zircons_Ttpath(time_years, Tt_paths_Temp, ZirconData=ZirconData)
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average  = zircon_age_PDF(ages_eruptible, number_zircons, bandwidth=1e5, n_analyses=n_analyses)

@test  prob[100] ≈ 5.536734758661114e-6

# A second way to generate the input is having Vector{Vector} for both time and Tt-path. 
time_years_vecs = Vector{Vector{Float64}}(undef,size(Tt_paths_Temp,2));
Tt_paths_Temp_vecs		= Vector{Vector{Float64}}(undef,size(Tt_paths_Temp,2));
for i=1:size(Tt_paths_Temp,2)
	if i<400
		time_years_vecs[i] 		= time_years
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[:,i]
	else
		# slightly adjust the size of the vectors
		time_years_vecs[i] 		= time_years[3:end]
		Tt_paths_Temp_vecs[i] 	= Tt_paths_Temp[3:end,i]
	end

end

# the calling routine is the same, but we get one additional output vector:
time_years1, prob1, ages_eruptible1, number_zircons1, T_av_time1, T_sd_time1, zircon_cumulativePDF1  = compute_zircons_Ttpath(time_years_vecs, Tt_paths_Temp_vecs)

@test  prob1[100] ≈ 5.536734758661114e-6


# Do the same but with a single routine that also returns the PDF's 
# Note that given the randomness, you'll always get different results
time_Ma, PDF_zircons, time_Ma_average, PDF_zircon_average, time_years, 
	prob2, ages_eruptible, number_zircons2, T_av_time, T_sd_time, zircon_cumulativePDF = compute_zircon_age_PDF(time_years_vecs, Tt_paths_Temp_vecs)

@test  prob2[100] ≈ 5.536734758661114e-6


#=	
    # Plot Zircon age probability distribution
	# these are the ploting routines using Makie, which is currently not a dependency of GeoParams (but may become one @ some stage)

	f = Figure()
	Axis(f[1, 1], xlabel = "Age [Myr]", ylabel = "Kernel density [ ]", title = "Zircon age probability distribution")
	for i in 1:length(PDF_zircons)
		CairoMakie.lines!(time_Ma[i]/1e6, PDF_zircons[i], color="gray66",linewidth=0.25)
	end
	CairoMakie.lines!(time_Ma_average/1e6, PDF_zircon_average, color="grey0",linewidth=2.)
	CairoMakie.xlims!(-1e5/1e6,1.5e6/1e6)
	save("Zircon_probability_plot_800kyrs.png",f)

	# plot evolution of the average temperature since the onset of magma injection
	f = Figure()
	Axis(f[1, 1], ylabel = "Temperature [°C]", xlabel = "Time [Myr]", title = "Evolution of the average temperature of the magmatic system")
	CairoMakie.lines!(time_years./1e6, T_av_time, color="grey0",linewidth=0.5)
	f
	save("Average_T.png",f)
=#




end

