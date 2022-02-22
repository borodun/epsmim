#!/bin/sh

mkdir build
mkdir results

echo Building
make

Nx=600
Ny=600
Nt=1500
Sx=1
Sy=1
output="results/out"

echo Simulating
./build/wavesim $Nx $Ny $Nt $Sx $Sy $output

for FILE in results/*.bin; do
    echo Plotting "$FILE"
    gnuplot <<- EOF
        nx=${Nx}
        ny=${Ny}
        filename="${FILE}"
        set terminal png size 1000,1000
        set output filename.".png"
        set xrange[-1:nx]
        set yrange[-1:ny]
        set size ratio 1
        set autoscale fix
        set palette gray
        set title filename
        plot filename binary array=(ny,nx)  format="%lf" with image
EOF
done

echo ""
echo Finished



