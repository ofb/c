set term pngcairo size 1000,700
set output 'absoluteMagnitude.png'
set grid
set xlabel "prime"
set ylabel "magnitude over all chi and lambda in {0,1,2,3}"
set key left top
plot 'matho/absoluteMagnitude.dat' using 1:2 with points
quit