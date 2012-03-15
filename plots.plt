unset xtics
unset ytics
#
set style line 1 lt rgb "#ff0000" lw 2 ps 1 pt 0
set style line 2 lt rgb "#ee1122" lw 2 ps 1 pt 0
set style line 3 lt rgb "#dd2244" lw 2 ps 1 pt 0
set style line 4 lt rgb "#cc3366" lw 2 ps 1 pt 0
set style line 5 lt rgb "#bb4488" lw 2 ps 1 pt 0
set style line 6 lt rgb "#aa55aa" lw 2 ps 1 pt 0
set style line 7 lt rgb "#9966cc" lw 2 ps 1 pt 0
set style line 8 lt rgb "#8877ee" lw 2 ps 1 pt 0
set style line 9 lt rgb "#7788ff" lw 2 ps 1 pt 0
#
plot "pp0_1.dat" using 1:2 ls 1 title "pp = 0.1" w linespoints, \
"pp0_2.dat" using 1:2 ls 2 title "pp = 0.2" w linespoints, \
"pp0_3.dat" using 1:2 ls 3 title "pp = 0.3" w linespoints, \
"pp0_4.dat" using 1:2 ls 4 title "pp = 0.4" w linespoints, \
"pp0_5.dat" using 1:2 ls 5 title "pp = 0.5" w linespoints, \
"pp0_6.dat" using 1:2 ls 6 title "pp = 0.6" w linespoints, \
"pp0_7.dat" using 1:2 ls 7 title "pp = 0.7" w linespoints, \
"pp0_8.dat" using 1:2 ls 8 title "pp = 0.8" w linespoints, \
"pp0_9.dat" using 1:2 ls 9 title "pp = 0.9" w linespoints
#
pause -1