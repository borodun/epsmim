#!/usr/bin/gnuplot
nx=ARG1
ny=ARG2
filename="${ARG3}"
set terminal png size 1000,1000
set output filename.".png"
set xrange[-1:nx]
set yrange[-1:ny]
set palette gray
set title filename
plot filename binary array=(ny,nx) format="%lf" with image