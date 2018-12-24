do for [t=0:3] {
set terminal wxt t+1

dataFile = 'val0.dat'
dataFile = sprintf('val%d.dat',t)

set view 67,18
set hidden3d
set palette defined (0 "purple",0.2 "blue", 0.3 "green",0.55 "yellow",0.7 "orange",0.85"red",1 "brown")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0

set xlabel "x-plane"
set ylabel "y-plane"
set xlabel "vertical-plane"

splot dataFile using 1:2:3:4 with pm3d lc rgb 1
}

a=4

{
set terminal wxt a+1

dataFile = 'initAll.dat'

set view 83,87
set hidden3d
set palette defined (0 "purple",0.2 "blue", 0.3 "green",0.55 "yellow",0.7 "orange",0.85"red",1 "brown")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0

set xlabel "x-plane"
set ylabel "y-plane"
set xlabel "vertical-plane"

splot dataFile using 1:2:3:4 with pm3d lc rgb 1
}

{
set terminal wxt a+2

dataFile = 'endAll.dat'

set view 83,87
set hidden3d
set palette defined (0 "purple",0.2 "blue", 0.3 "green",0.55 "yellow",0.7 "orange",0.85"red",1 "brown")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0

set xlabel "x-plane"
set ylabel "y-plane"
set xlabel "vertical-plane"

splot dataFile using 1:2:3:4 with pm3d lc rgb 1
}
