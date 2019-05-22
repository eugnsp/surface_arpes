#!/usr/bin/gnuplot

set terminal pngcairo size 1000,800 
set output 'arp.png'
set pm3d map
set title "ARPES spectrum"
set xlabel "k_x [Ang^{-1}]"
set ylabel "E [eV]"
set tics out nomirror
set cbrange [0:]
set autoscale fix
unset cbtics
set palette rgb -21,-22,-23
splot 'arp.dat' binary matrix
