#!/usr/bin/gnuplot

set terminal pngcairo size 1000,800 
set output 'arpes_kx_kz.png'
set pm3d map
set title "ARPES spectrum"
set xlabel "k_x [Å^{-1}]"
set ylabel "k_z [Å^{-1}]"
set cbrange [0:]
set autoscale fix
unset cbtics
set palette rgb -30,-31,-32
splot 'arpes_kx_kz.dat' binary matrix using 2:1:3 notitle
