#!/usr/bin/gnuplot
#
#set title "Correct Binds Found"
#set xlabel "Percent of Dataset"
#set ylabel "Number of Actual Binds Identified"
#
#set datafile missing "-"
#plot 'accuracy.txt' using 1:3 with lines ti 'Actual Binds Identified', x*379.12 with lines ti 'Uniform'
#
#pause mouse
#
#reset

set xrange[0:1]
set yrange[0:1]

set title "PharML.Bind Ensemble ROC Curves on DUD-E Dataset\n(Receiver Operating Characteristics)"
set xlabel "False Positive Rate (FP/(FP+TN))"
set ylabel "True Positive Rate (TP/(TP+FN))"

set xtics .1
set ytics .1

set key outside right

set datafile missing "-"
plot x with lines lw 1.5 lt 3 lc rgb 'red' ti 'Random', \
     'roc1.txt' using 2:1 with lines lw 1.0 lc rgb '#FF5500' ti 'e_1 ROC (AUC: .669)', \
     'roc2.txt' using 2:1 with lines lw 1.0 lc rgb '#DD7700' ti 'e_2 ROC (AUC: .709)', \
     'roc3.txt' using 2:1 with lines lw 1.0 lc rgb '#BB9900' ti 'e_3 ROC (AUC: .727)', \
     'roc4.txt' using 2:1 with lines lw 1.0 lc rgb '#99AA00' ti 'e_4 ROC (AUC: .739)', \
     'roc5.txt' using 2:1 with lines lw 2.0 lc rgb '#00FF00' ti 'e_5 ROC (AUC: .741)'

pause mouse

reset

set yrange[0:1]

set title "PharML.Bind Ensemble AUC-ROC Scores on DUD-E Dataset"
set ylabel "AUC-ROC (Higher Is Better)"
set xlabel "DUD-E Protein ID"
set ytics .1

set key outside right

set datafile missing "-"
set style data histogram
set boxwidth 1.5 relative
set style fill solid border -1
set xtics rotate by 90 right

plot 'srt_per_protein.txt' using 4:xtic(1) noti

pause mouse
