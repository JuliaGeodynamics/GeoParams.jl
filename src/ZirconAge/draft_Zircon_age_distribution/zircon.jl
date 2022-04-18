# Julia translation of the of the R-script provided by Gregor Weber & Luca Caricchi 
# Used in the publication "Estimating the current size and state of subvolcanic magma reservoirs"

# 15/04/2022

using CSV
using DataFrames
using Loess
using Statistics
using StatsBase													
using KernelDensity
using CairoMakie
using Test
using Parameters

"""

	Declare function to calculate zircon fraction as function of a temperature profile
"""
function zircon_fraction(T, max_x_zr) 
	A = (1.62.-1.8*(10^4)*exp.((-10^4)./(T .+ 273.15))).*max_x_zr 
	A[A .<= 0.0] .= 0.0

	return A
end

@with_kw_noshow struct ZirconAgeData
	Tsat::Float64 			= 825.0		# Maximum zircon saturation temperature [C]
	Tmin::Float64 			= 690.0 	# [C] Minimum zircon saturation Temperature [C]
	Tsol::Float64 			= 690.0		# [C] Solidus temperature [C]
	Tcal_max::Float64		= 800.0		# max temperature to calculate zircon fraction
	Tcal_step::Float64 		= 1.0		# temperature step to caclulate zircon fraction (resolution of Zircon saturation curve discretization)
	max_x_zr::Float64 		= 0.001		# max fraction zircons at solidus
	zircon_number::Int64    = 100.0		# number of required zircons 
end

"""
	Loess fit of zircon number over saturation range
"""
function loess_fit_zircon_sat(ZirconData::ZirconAgeData)
	@unpack Tsol, Tcal_step, Tcal_max, Tsat, max_x_zr, zircon_number, Tmin = ZirconData								

	# Get cumulative Zircon fraction (after Tierney et al., 2016; Geology)
	T 				= collect(Float64, Tsol:Tcal_step:Tcal_max)
	m_steps 		= length(T)

	Tfit 			= range(Tsol, stop = Tsat, length = m_steps)

	n_zircon 		= Vector{Float64}(undef,m_steps)				# number of zircon

	x_zircon 		= zircon_fraction(T,max_x_zr)

	n_zircon[2:length(x_zircon)] = -diff(x_zircon)					# differentiate zircon fraction
	n_zircon[1] 	= n_zircon[2]									# set first value

	n_zircon		= ceil.(((n_zircon*zircon_number)/maximum(n_zircon)).-minimum(floor.((n_zircon*zircon_number)/maximum(n_zircon))))

	# fit n_zircon_N with loess regression
	n_zircon_N_fit 	= loess(Tfit, n_zircon, span=1.0)

	return n_zircon_N_fit
end



function compute_zirconsaturation(time_years, Tt_paths, ZirconData::ZirconAgeData)

	@unpack Tmin, Tsat = ZirconData

	#time_step		= time_years[2]-time_years[1]				# timestep in years (might be good to generate an array for this in case timestep is not a constant)
	Δt 				=	diff(time_years)[1]

	# get zircon number / zircon saturation temperature loess fit
	n_zircon_N_fit  = loess_fit_zircon_sat(ZirconData)
	
	time_er_min 	= maximum(time_years)							# backward count
	
	# In the end here they only remove the last entry...
	@show size(time_years), Δt, size(Tt_paths), time_er_min
	ID_row_er 		= findall( time_years .== maximum(time_years[time_years .< time_er_min]) )
	@show ID_row_er

	# Here the matrix is reduced to only contain the temperature informations
	Tt_paths_Temp 	= Tt_paths[1:ID_row_er[1],2:size(Tt_paths,2)]	
	
	# find all the Tt paths that go through the zircon saturation range
	ID_col_er 		= findall( (Tt_paths_Temp[ID_row_er,:] .> Tmin) .& (Tt_paths_Temp[ID_row_er,:] .< Tsat))
	@show size(ID_col_er)
	@test sum(Tt_paths_Temp) == 2.4995483211e8
	@test sum(ID_col_er)[2] == 111995
	
	# find the number of timesteps during which the temperature is > Tmin and < Tsat
	length_trace 	= Vector{Float64}(undef,length(ID_col_er))
	for i in 1:length(ID_col_er)
		length_trace[i]    = length( findall( (Tt_paths_Temp[:,ID_col_er[i][2]] .> Tmin) .& (Tt_paths_Temp[:,ID_col_er[i][2]] .< Tsat)) )
	end
	
	# the next several lines can likely be achieved in a more elegant way...
	id				= findall( length_trace .== maximum(length_trace)) 
	ID_col_lgst_tr 	= ID_col_er[ id[1] ][2]
	
	id 				= findall( Tt_paths_Temp[:,ID_col_lgst_tr] .< Tsat) 
	VALID_min_time 	= findmin(Tt_paths_Temp[id,ID_col_lgst_tr]) 
	ID_min_time		= VALID_min_time[2]
	
	max_age_spread	= maximum(length_trace)*Δt					# This is defined among all selected paths
	n_zr			= similar(Tt_paths_Temp, Float64) .= 0.0
	
	# Had to use 'for' loops because Loess prediction model do not work outside regression bounds
	# This part calculates the number of zircon generated at each timestep
	for i in 1:size(Tt_paths_Temp,2)
		for j in 1:size(Tt_paths_Temp,1)
			if (Tt_paths_Temp[j,i] > Tsol) .& (Tt_paths_Temp[j,i] < Tsat)
				n_zr[j,i]    = floor.(Loess.predict(n_zircon_N_fit, Tt_paths_Temp[j,i]))
			else
				n_zr[j,i]    = 0.0
			end
		end
	end
	
	T_av_time_1 	= replace!(Tt_paths_Temp, 0.0 => NaN)
	T_av_time 		= Vector{Float64}(undef,size(T_av_time_1,1)) .= 0.0
	sd_time 		= Vector{Float64}(undef,size(T_av_time_1,1)) .= 0.0
	
	# get the average temperature of the Tt paths and the standard deviation
	for i in 1:size(Tt_paths_Temp,1)
		T_av_time[i]= mean(filter(!isnan, T_av_time_1[i,:]))
		sd_time[i] 	= std(filter( !isnan, T_av_time_1[i,:]))
	end
	
	
	# I clarified the R function because the minimum step length to grow a zircon is simply a ratio of the maximum trace between Tmin and Tsol
	# This makes sense as we only deal with fractions here. Because no mass is provided the real zircon size cannot possibly be determined
	min_step_n		= floor( (time_zr_growth/max_age_spread)*(max_age_spread/Δt) )
	
	# find the Tt paths that have a number of timestep in the saturation range greater than the defined min_step_n
	id				= findall( length_trace .> min_step_n) 
	ID_col_er_1		= getindex.(ID_col_er[id], [2])
	
	int_zr_sat		= collect(Float64,  ID_min_time:1.0:(time_er_min/Δt)-min_step_n)
	int_zr_sat		= floor.(Int64,int_zr_sat)
	
	T_av_time_slct 	= Vector{Float64}(undef,length(int_zr_sat)-1) .= 0.0
	
	for i in 1:(size(int_zr_sat,1)-1)
		id2				= ID_col_er_1[ findall( (Tt_paths_Temp[int_zr_sat[i],ID_col_er_1[:]] .> Tmin) .& (Tt_paths_Temp[int_zr_sat[i],ID_col_er_1[:]] .< Tsat)) ]
		if isempty(id2) == true
			T_av_time_slct[i] = NaN
		else
			T_av_time_slct[i] = median(filter(!isnan, Tt_paths_Temp[int_zr_sat[i],id2]))
		end
	end
	
	replace!(Tt_paths_Temp, NaN => 0.0)
	ID_col_er			= getindex.(ID_col_er, [2])
	for i in 1:length(ID_col_er)
		k 				= maximum(findall( (Tt_paths_Temp[:,ID_col_er[i]]) .== 0.0 ))
		Tt_paths_Temp[1:k,ID_col_er[i]] .= 0.0
	end

	zr_select			= similar(Tt_paths_Temp, Float64) .= 0.0
	zr_select[Tt_paths_Temp .> 0.0] .= 1.0	
	n_zrc2_0			= zr_select.*n_zr
	n_measurable_ages 	= sum(n_zrc2_0[:,ID_col_er_1], dims=2)
	sz 					= size(n_zrc2_0[:,ID_col_er_1],1)
	ages_eruptible		= collect(Float64,  1.0:Δt:sz*Δt)
	
	return Tt_paths_Temp, n_measurable_ages, ages_eruptible, n_zr, ID_col_er_1, n_zrc2_0, T_av_time
end

# output
saveplot 		= 1

# declare constant variables
s2y 			= 365.0*24.0*3600.0 							# second to year

#=
# Zircons informations
Tsat 			= 825.0 										# Maximum zircon saturation temperature [C]
Tmin 			= 690.0 										# [C] Minimum zircon saturation Temperature [C]
Tsol 			= 690.0 										# [C] Solidus temperature [C]
Tcal_max		= 800.0											# max temperature to calculate zircon fraction
Tcal_step 		= 1.0											# temperature step to caclulate zircon fraction (resolution of Zircon saturation curve discretization)
max_x_zr 		= 0.001 										# max fraction zircons at solidus
zircon_number	= 100.0											# number of wanted zircons 
=#

ZirconData  = ZirconAgeData();		# default data

# Clarified way to define the zircon growth as a minimum time within T saturation range (This is what the method used in the R script, boils down too)
# -> remain in the Zr saturation zone more than 1/3 of the time the Tt path with the longest time in the saturation zone
time_zr_growth 	= 0.7e6											# ref should be 0.7

# Parameters for Zircon statistical analysis
n_analyses		= 300											# number of synthetic zircon analyses
n				= 100											# number of repetitions
n_samples		= 80											# number of sampled Tt paths


# read input file (this step should be skipped, when using tracers)		
filename 		= "Tt15_st_Bl_rad5_0.0126.txt"					# path to input file 
Tt_paths		= Matrix(CSV.read(filename, DataFrame, header=0, copycols=true))
n_paths			= size(Tt_paths,2)-1							# -1 ro remove first column that corresponds to the time
time_years 		= Tt_paths[:,1]./s2y							# time fron seconds to years

#=
#@unpack Tmin, Tsat = ZirconData

time_step		= (Tt_paths[2,1]-Tt_paths[1,1])/s2y				# timestep in years (might be good to generate an array for this in case timestep is not a constant)


# get zircon number / zircon saturation temperature loess fit
n_zircon_N_fit = loess_fit_zircon_sat(ZirconData)


time_er_min 	= maximum(time_years)							# backward count

# In the end here they only remove the last entry...
ID_row_er 		= findall( time_years .== maximum(time_years[time_years .< time_er_min]) )

# Here the matrix is reduced to only contain the temperature informations
Tt_paths_Temp 	= Tt_paths[1:ID_row_er[1],2:size(Tt_paths,2)]	

# find all the Tt paths that go through the zircon saturation range
ID_col_er 		= findall( (Tt_paths_Temp[ID_row_er,:] .> Tmin) .& (Tt_paths_Temp[ID_row_er,:] .< Tsat))

@test sum(Tt_paths_Temp) == 2.4995483211e8
@test sum(ID_col_er) == 111995

# find the number of timesteps during which the temperature is > Tmin and < Tsat
length_trace 	= Vector{Float64}(undef,length(ID_col_er))
for i in 1:length(ID_col_er)
	length_trace[i]    = length( findall( (Tt_paths_Temp[:,ID_col_er[i][2]] .> Tmin) .& (Tt_paths_Temp[:,ID_col_er[i][2]] .< Tsat)) )
end

# the next several lines can likely be achieved in a more elegant way...
id				= findall( length_trace .== maximum(length_trace)) 
ID_col_lgst_tr 	= ID_col_er[ id[1] ][2]

id 				= findall( Tt_paths_Temp[:,ID_col_lgst_tr] .< Tsat) 
VALID_min_time 	= findmin(Tt_paths_Temp[id,ID_col_lgst_tr]) 
ID_min_time		= VALID_min_time[2]

max_age_spread	= maximum(length_trace)*time_step				# This is defined among all selected paths
n_zr			= similar(Tt_paths_Temp, Float64) .= 0.0

# Had to use 'for' loops because Loess prediction model do not work outside regression bounds
# This part calculates the number of zircon generated at each timestep
for i in 1:size(Tt_paths_Temp,2)
	for j in 1:size(Tt_paths_Temp,1)
		if (Tt_paths_Temp[j,i] > Tsol) .& (Tt_paths_Temp[j,i] < Tsat)
			n_zr[j,i]    = floor.(Loess.predict(n_zircon_N_fit, Tt_paths_Temp[j,i]))
		else
			n_zr[j,i]    = 0.0
		end
	end
end



T_av_time_1 	= replace!(Tt_paths_Temp, 0.0 => NaN)
T_av_time 		= Vector{Float64}(undef,size(T_av_time_1,1)) .= 0.0
sd_time 		= Vector{Float64}(undef,size(T_av_time_1,1)) .= 0.0

# get the average temperature of the Tt paths and the standard deviation
for i in 1:size(Tt_paths_Temp,1)
	T_av_time[i]= mean(filter(!isnan, T_av_time_1[i,:]))
	sd_time[i] 	= std(filter( !isnan, T_av_time_1[i,:]))
end


# I clarified the R function because the minimum step length to grow a zircon is simply a ratio of the maximum trace between Tmin and Tsol
# This makes sense as we only deal with fractions here. Because no mass is provided the real zircon size cannot possibly be determined
min_step_n		= floor( (time_zr_growth/max_age_spread)*(max_age_spread/time_step) )

# find the Tt paths that have a number of timestep in the saturation range greater than the defined min_step_n
id				= findall( length_trace .> min_step_n) 
ID_col_er_1		= getindex.(ID_col_er[id], [2])

int_zr_sat		= collect(Float64,  ID_min_time:1.0:(time_er_min/time_step)-min_step_n)
int_zr_sat		= floor.(Int64,int_zr_sat)

T_av_time_slct 	= Vector{Float64}(undef,length(int_zr_sat)-1) .= 0.0

for i in 1:(size(int_zr_sat,1)-1)
	id2				= ID_col_er_1[ findall( (Tt_paths_Temp[int_zr_sat[i],ID_col_er_1[:]] .> Tmin) .& (Tt_paths_Temp[int_zr_sat[i],ID_col_er_1[:]] .< Tsat)) ]
	if isempty(id2) == true
		T_av_time_slct[i] = NaN
	else
		T_av_time_slct[i] = median(filter(!isnan, Tt_paths_Temp[int_zr_sat[i],id2]))
	end
end

replace!(Tt_paths_Temp, NaN => 0.0)
ID_col_er			= getindex.(ID_col_er, [2])
for i in 1:length(ID_col_er)
	k 				= maximum(findall( (Tt_paths_Temp[:,ID_col_er[i]]) .== 0.0 ))
	Tt_paths_Temp[1:k,ID_col_er[i]] .= 0.0
end


zr_select			= similar(Tt_paths_Temp, Float64) .= 0.0
zr_select[Tt_paths_Temp .> 0.0] .= 1.0	
n_zrc2_0			= zr_select.*n_zr
n_measurable_ages 	= sum(n_zrc2_0[:,ID_col_er_1], dims=2)
sz 					= size(n_zrc2_0[:,ID_col_er_1],1)
ages_eruptible		= collect(Float64,  1.0:time_step:sz*time_step)
age_resampled       = Matrix{Float64}(undef,n_analyses,n) .= 0.0
=#


Tt_paths_Temp, n_measurable_ages, ages_eruptible, n_zr, ID_col_er_1,n_zrc2_0,T_av_time 		= compute_zirconsaturation(time_years, Tt_paths, ZirconData)

prob = n_measurable_ages/sum(n_measurable_ages)
prob = prob[:,1]


age_resampled       = Matrix{Float64}(undef,n_analyses,n) .= 0.0
for i in 1:n
	age_resampled[:,i] = sample(ages_eruptible, Weights(prob), n_analyses, replace=true)
end

# calculate standard deviation of the age span
sd_age_resampled	= Vector{Float64}(undef,n) .= 0.0

for i in 1:n
	sd_age_resampled[i] = std(age_resampled[:,i])
end

lower_bound 		= (mean(sd_age_resampled*2)-2.0*std(sd_age_resampled*2))/1000.0
upper_bound 		= (mean(sd_age_resampled*2)+2.0*std(sd_age_resampled*2))/1000.0
mean_sd_resampled 	= (mean(sd_age_resampled*2.0))/1000.
Sigma_age_span 		= (mean_sd_resampled,lower_bound,upper_bound)


# add tests to check that results remain consistent
@test sum(n_zr[:,200]) == 28117.0
@test  mean(Tt_paths_Temp) == 487.0858188593903
#@test  mean(age_resampled) == 996532.3057680113

#=
@test Sigma_age_span[1] ≈ 635.4927816736695
@test Sigma_age_span[2] ≈ 604.8156882669368
@test Sigma_age_span[3] ≈ 672.7800193854439
@test sum(age_resampled) ≈ 2.990029563628869e10
=#



# Plot Zircon age propability distribution
if saveplot == 1
	f = Figure()
	Axis(f[1, 1], xlabel = "Age [Myr]", ylabel = "Kernel density [?]", title = "Zircon age propability distribution")
	for i in 1:length(ID_col_er_1)
		n_meas 			= n_zrc2_0[:,ID_col_er_1[i]]
		px	  			= n_meas/sum(n_meas)
		smp				= sample( (maximum(ages_eruptible) .- ages_eruptible)./1e6, Weights(px), n_analyses, replace=true)
		y 				= kde(smp, bandwidth=1e5/1e6)
		CairoMakie.plot!(y, color="gray66",linewidth=0.25)
	end

	pxAv	  			= n_measurable_ages[:,1]./sum(n_measurable_ages[:,1])
	smpAv				= sample( (maximum(ages_eruptible) .- ages_eruptible)./1e6, Weights(pxAv), n_analyses, replace=true)
	yAv 				= kde(smpAv, bandwidth=1e5/1e6)
	CairoMakie.plot!(yAv, color="grey0",linewidth=2.)
	CairoMakie.xlims!(-1e5/1e6,1.5e6/1e6)
	f
	save("Zircon_propability_plot_800kyrs.png",f)

	# plot evolution of the average temperature since the onset of magma injection
	f = Figure()
	Axis(f[1, 1], xlabel = "Temperature [°C]", ylabel = "Time [Myr]", title = "Evolution of the average temperature of the magmatic system")
	CairoMakie.lines!(time_years[1:length(time_years)-1]./1e6, T_av_time, color="grey0",linewidth=0.5)
	f
	save("Average_T.png",f)
end
