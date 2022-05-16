using Test
using GeoParams, LinearAlgebra
@testset "TASclassification.jl" begin

	ClassTASdata    = TASclassificationData();	 # default data

	# list of test points S-TA
	point_list =[
		60 5;
		50 14;
		70 2;
		41 2;
		42 6.5;
		56 15;
		90 7;
		57 5.9;
		50 8.1;
		58 11.7;
		49 7.3;
		54 4.5;
		50 6;
		49 3.4;
		54 12;
		60 7;
		40 10;
		67 9.4;
		53.1 6.4;
		70.5 12.1;
		38 4.4		];

	field = [13 1 14 2 3 6 15 9 4 11 3 12 8 7 5 10 1 11 9 15 1];
	name  = ["andesite"	"foidite"	"dacite"	"picrobasalt"	"basanite"	"phonolite"	"rhyolite"	"basaltic trachyandesite"	"phonotephrite"	"trachyte"	"basanite"	"basaltic andesite"	"trachybasalt"	"basalt"	"tephriphonolite"	"trachyandesite"	"foidite"	"trachyte"	"basaltic trachyandesite"	"rhyolite"	"foidite"];
	# add tests to check that results are consistent
	for i = 1:size(point_list,1)
		index = computeTASclassification(point_list[i,:], ClassTASdata=ClassTASdata);
		@test  index == field[i];
		@test retrieveTASrockType(index, ClassTASdata=ClassTASdata) == name[i];
	end

end

