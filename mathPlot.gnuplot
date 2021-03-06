set term png size 1000,700
set output 'matho/absoluteMagnitude.png'
set grid
set xlabel "prime"
set ylabel "magnitude over all chi and lambda in {0,1,2,3}"
set key left top
set pointsize .5
plot 'matho/lambda-0.dat' using 1:2 with points pt 7, 'matho/lambda-1.dat' using 1:2 with points pt 7, 'matho/lambda-2.dat' using 1:2 with points pt 7, 'matho/lambda-3.dat' using 1:2 with points pt 7, 'matho/absoluteMagnitude.dat' using 1:2 with points pt 7

set output 'matho/lambda0.png'
plot 'matho/lambda-0.dat' using 1:2 with points
set output 'matho/lambda1.png'
plot 'matho/lambda-1.dat' using 1:2 with points
set output 'matho/lambda2.png'
plot 'matho/lambda-2.dat' using 1:2 with points
set output 'matho/lambda3.png'
plot 'matho/lambda-3.dat' using 1:2 with points

# set output 'matho/lambda3D-1.png'
# unset hidden3d
# set view 60,300
# set autoscale
# set style data points
# set key box
# splot "matho/lambda3D-1.dat"

quit