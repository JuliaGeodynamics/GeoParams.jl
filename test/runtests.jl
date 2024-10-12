using Test, Statistics


function runtests()
    files = readdir(@__DIR__)
    test_files = filter(startswith("test_"), files)

    for f in test_files[1]
        if !isdir(f)
            include(f)
        end
    end
end

runtests()

