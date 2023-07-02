
using UVITTools
using Documenter

#push!(LOAD_PATH,"../src")
# println(LOAD_PATH)

makedocs(
         sitename = "UVITTools.jl",
         modules  = [UVITTools],
         pages=[
                "Home" => "index.md"
               ])
deploydocs(;
    repo="github.com/gulabd/UVITTools.jl.git",
    devbranch = "main"
)