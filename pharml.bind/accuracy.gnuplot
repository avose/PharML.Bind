#/usr/bin/gnuplot

set title "Correct Binds Found"
set xlabel "Percent of Dataset"
set ylabel "Number of Actual Binds Identified"

set datafile missing "-"
plot 'accuracy.txt' using 1:3 with lines ti 'Actual Binds Identified', x*379.12 with lines ti 'Uniform'

pause mouse

reset

set xrange[0:1]
set yrange[0:1]

set title "ROC Curve"
set xlabel "FPR"
set ylabel "TPR"

set datafile missing "-"
plot 'roc.txt' using 2:1 with lines ti 'ROC', x with lines ti 'y=x'

pause mouse
