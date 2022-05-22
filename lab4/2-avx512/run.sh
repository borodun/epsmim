#!/bin/sh

rm -r build
mkdir build
mkdir results -p

echo Building
make

Nx=10000
Ny=10000
Nt=120
Sx=5000
Sy=5000
output="results/out"

for threads in 1 2 4 8
do
  ./build/wavesim $Nx $Ny $Nt $Sx $Sy $threads $output
done

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



