using LaTeXStrings
using DataFrames
using Unitful
using Parameters
using ..Units
using ..MaterialParameters.MaterialParams
using GeoParams: AbstractMaterialParam, AbstractMaterialParamsStruct

function ExpDataAsTable(p::AbstractMaterialParamsStruct)
    ### table1.jl script
    T[:,:Name] = "\$\\" * T[:,:Description] * "\$";  # replace symbol with latex math

    # Generate table header
    Table  = "\\begin{tabular}{| " * "c |" ^ (ncol(T)+1) * "}\n";
    Table *= "    \\hline\n";
    Table *= "    % Table header\n";
    Table *= "    \\rowcolor[gray]{0.9}\n";
    Table *= "  Row"; for i in 1:ncol(T); Table *= " & " * string(names(T)[i]); end; Table *= " \\\\\n";

    # Generate table body (with nice alternating row colours)
    toggleRowColour(x) = x == "0.8" ? "0.7" : "0.8";
    rowcolour = toggleRowColour(0.7);

    Table *= "    % Table body\n";
    for row in 1 : nrow(T)
        Table *= "  \\rowcolor[gray]{" * (rowcolour = toggleRowColour(rowcolour); rowcolour) * "}\n";
        Table *= "  " * string(row); for col in 1 : ncol(T) Table *= " & " * string(T[row,col]); end; Table *= " \\\\\n";
        Table *= "  \\hline\n"; 
    end
    Table *= "\\end{tabular}\n";

    # Export result to .tex file
    write("table1.tex", Table);
end