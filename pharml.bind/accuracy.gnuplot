#/usr/bin/gnuplot

set title "Correct Binds Found"
set xlabel "Percent of Dataset"
set ylabel "Number of Actual Binds Identified"

set datafile missing "-"
plot 'accuracy.txt' using 1:3 with lines ti 'Actual Binds Identified', x*379.12 with lines ti 'Uniform'

pause mouse
