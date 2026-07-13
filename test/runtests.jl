using GeoParams
using Test, Statistics, ParallelTestRunner
using LaTeXStrings


function runtests(args)
    printstyled("Testing package GeoParams.jl\n"; bold = true, color = :white)

    testsuite = find_tests(@__DIR__)
    pt_args = ParallelTestRunner.parse_args(args)
    nfail = 0

    try
        ParallelTestRunner.runtests(GeoParams, pt_args; testsuite)
    catch ex
        nfail += 1
    end

    return nfail
end

exit(runtests(ARGS))
