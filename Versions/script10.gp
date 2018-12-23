set terminal wxt

dataFile = 'val3.dat'

set view 67,18
set hidden3d
set palette defined (0 "purple",0.2 "blue", 0.3 "green",0.55 "yellow",0.7 "orange",0.85"red",1 "brown")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0

splot dataFile using 1:2:3:4 with pm3d lc rgb 1


