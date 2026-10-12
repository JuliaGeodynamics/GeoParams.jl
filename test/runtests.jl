using GeoParams
using Test, Statistics, ParallelTestRunner
using LaTeXStrings


function runtests(args)
    printstyled("Testing package GeoParams.jl\n"; bold = true, color = :white)

    testsuite = find_tests(@__DIR__)
    pt_args = ParallelTestRunner.parse_args(args; custom = ["backend"])
    backend = something(pt_args.custom["backend"], "CPU")
    backend in ("CPU", "CUDA", "AMDGPU", "Metal") ||
        error("unknown --backend=$backend; use CPU, CUDA, AMDGPU or Metal")
    # GPU agents run only the GPU suite; test_GPU.jl reads the backend in its worker
    ENV["GEOPARAMS_TEST_BACKEND"] = backend
    backend == "CPU" ? delete!(testsuite, "test_GPU") : filter!(p -> p.first == "test_GPU", testsuite)
    nfail = 0

    try
        ParallelTestRunner.runtests(GeoParams, pt_args; testsuite)
    catch ex
        nfail += 1
    end

    return nfail
end

exit(runtests(ARGS))
