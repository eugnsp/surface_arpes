# ARPES spectra calculator

## Synopsis

This is a code to calculate ARPES spectra of slabs with an accumulation
or a depletion surface layer. The calculation proceeds in two steps:

1. The Poisson and the Schr&ouml;dinger equations are solved self-consistently
to obtain eigenpairs {<code>E<sub>n</sub></code>, <code>&psi;<sub>n</sub>(z)</code>}.
2. ARPES spectra are then calculated using a Fourier transform of
<code>&psi;<sub>n</sub>(z)</code>, followed by a convolution with normal
distributions to account for instrumental broadening.

1D finite elements ([`es_fe` library](https://github.com/eugnsp/es_fe))
are used to discretize the Poisson and the Schr&ouml;dinger equations.
The Intel MKL library is used to solve linear systems, generalized eigenvalue
system, compute Fourier transforms and convolution. The mathematical details
can be found in [this PDF file](doc/model.pdf).

The results are exported into Matlab/Octave MAT-files and Gnuplot binary
matrix files.

GCC 8.3 was used to compile the code. Clang doesn't work due to an `auto`
non-type template parameter deduction [bug](https://stackoverflow.com/questions/56125811/auto-non-type-template-parameter-ambiguous-partial-specializations-in-clang).
It can easily be overcome, if need be.

## Results

All images below are presented for exposition only, no attempt has been
made here to fit any experimental data.

Accumulation layer with one bound state:

![Accumulation layer with one bound state](example/accum1_sm.png)
[Large image](example/accum1.png)

Accumulation layer with two bound states:

![Accumulation layer with two bound states](example/accum2_sm.png)
[Large image](example/accum2.png)

Depletion layer:

![Depletion layer](example/depl_sm.png)
[Large image](example/depl.png)

## How to run

**NB**: this code is under development, and is not a ready-to-use tool for an
end-user.

Set `MKLROOT` environment variable to point to the MKL installation directory.
Then:

```sh
git clone --recursive https://github.com/eugnsp/surface_arpes.git
cd surface_arpes
mkdir build
cd build
cmake .. && make

./arpes < ../example/accum1.txt
../plot/plot_all.sh
```

## References

1. V.N.Strocov. *Photoemission response of 2D states.*\
[arXiv preprint (2018)](https://arxiv.org/abs/1801.07505).
2. S.Moser et al. *How to extract the surface potential profile
from the ARPES signature of a 2DEG.*\
[J. Electron. Spectrosc. **225**, 16 (2018)](https://doi.org/10.1016/j.elspec.2018.01.008).
3. A.Trellakis et al. *Iteration scheme for the solution of the
two-dimensional Schr&ouml;dinger-Poisson equations in quantum
structures*.\
[J. Appl. Phys. 81, 7880 (1997)](https://doi.org/10.1063/1.365396).
