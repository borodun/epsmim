#!/bin/sh

mkdir build
mkdir results

make

Nx=100
Ny=100
Nt=100
Sx=50
Sy=50
filename="results/out"

./build/wavesim $Nx $Ny $Nt $Sx $Sy $filename

for FILE in results/*.bin; do
    gnuplot <<- EOF
        nx=${Nx}
        ny=${Ny}
        filename="${FILE}"
        set terminal png size 1000,1000
        set output filename.".png"
        set xrange[-1:nx]
        set yrange[-1:ny]
        set palette gray
        set title filename
        plot filename binary array=(ny,nx) format="%lf" with image
EOF
done



