#!/usr/bin/gnuplot

set terminal pngcairo size 1800,600 
set output 'arpes.png'
set multiplot layout 1,3 title "ARPES spectra" font "sans,20"
set pm3d map
set cbrange [0:]
set autoscale fix
unset cbtics
unset colorbox
set size square
set palette rgb -30,-31,-32

set xrange [-.3:.3]
set yrange [-.6:.4]
set xlabel "k_x [Å^{-1}]"
set ylabel "E [eV]"
splot 'arpes_kx_e.dat' binary matrix using 2:1:3 notitle

set xrange [-.3:.3]
set yrange [-.6:.4]
set xlabel "k_z [Å^{-1}]"
set ylabel "E [eV]"
splot 'arpes_e_kz.dat' binary matrix using 1:2:3 notitle

set xrange [-.3:.3]
set yrange [-.3:.3]
set xlabel "k_x [Å^{-1}]"
set ylabel "k_z [Å^{-1}]"
set colorbox
splot 'arpes_kx_kz.dat' binary matrix using 2:1:3 notitle
