set term pngcairo size 1000,700
set output 'absoluteMagnitude.png'
set grid
set xlabel "prime"
set ylabel "magnitude over all chi and lambda in {0,1,2,3}"
set key left top
plot 'matho/absoluteMagnitude.dat' using 1:2 with points

set output 'lambda1.png'
plot 'matho/lambda-1.dat' using 1:2 with points
set output 'lambda2.png'
plot 'matho/lambda-2.dat' using 1:2 with points
set output 'lambda3.png'
plot 'matho/lambda-3.dat' using 1:2 with points
set output 'lambda4.png'
plot 'matho/lambda-4.dat' using 1:2 with points

quit