set output "filename.ps"
set size xscale[-10,10],yscale[-10,10]
set xlabel "x axis"
set tics out
set ticslevel 10
splot "endVal.dat" every ::1::361 using 1:2:3 with lines


