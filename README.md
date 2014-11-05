# eiscor - Eigensolvers based on unitary CORe transformations
This package is a Julia wrapper of code accompanying the article written by Jared L. Aurentz, Thomas Mach, Raf Vandebril and David S. Watkins. 

## Installation
```
Pkg.clone("https://github.com/andreasnoackjensen/AMVW.jl")
Pkg.build("AMVW")
```
## Example: Roots of a polynomial of degree 10,000
```julia
julia> using AMVW

julia> p = Poly(randn(10000));

julia> @time AMVW.rootsAMVW(p);
elapsed time: 48.162882002 seconds (1987560 bytes allocated)
```
Don't try `roots(p)`
