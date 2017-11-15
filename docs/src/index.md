# NonlinearOptics.jl Documentation

This is a package for simulating nonlinear optics problems in Julia. The
purpose of this package is to supply a flexible interface for designing nonlinear
optic systems with full accounting for all physical effects, if desired.

The functionality provided includes support for:

* Supercontinuum generation
* Frequency combs (both Lugatio-Lefever and Ikeda Map formulation)

with plans to include more (Raman lasing, four-wave mixing, etc). The design
philosophy is to provide the user with a generic, powerful interface which
can be loaded with the user's concrete needs.

## Installation and First steps

To install the package, run the following command in the Julia REPL:
```julia
Pkg.clone("https://github.com/platawiec/NonlinearOptics.jl")
```

A variety of tutorials are provided in this documentation as well.
