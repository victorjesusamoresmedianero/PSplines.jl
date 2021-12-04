# PSplines.jl

The source code of the package is available in: [PSplines](https://github.com/victorjesusamoresmedianero/PSplines.jl).

## Objective of the PSplines.jl
The objective of this package is to provide basic functionality on calculation and manipulation of cubic periodic PSplines (Penalized BSplines). For the theoretical fundamentals, see the section on PSplines [Article](https://www.sciencedirect.com/science/article/pii/S0965997818310779).

## Using PSplines.jl
Let's imagine we have some exprimental whose x and y coordinates are defined by the vectors $\bm{x}$ and $\bm{y}$. We would like to construct a BSpline which passes through the experimental points and at the same time has some smoothing conditions. The first step is to define the knots of the BSpline that will fix the set of basis functions ($N_i\left(x\right)$) for the considered interval.