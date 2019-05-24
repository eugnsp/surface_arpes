#!/usr/bin/gnuplot

set terminal pngcairo size 1000,800 
set output 'arpes_kz_e.png'
set pm3d map
set title "ARPES spectrum"
set xlabel "k_z [Ang^{-1}]"
set ylabel "E [eV]"
set cbrange [0:]
set autoscale fix
unset cbtics
set palette rgb -30,-31,-32
splot 'arpes_e_kz.dat' binary matrix using 1:2:3 notitle
