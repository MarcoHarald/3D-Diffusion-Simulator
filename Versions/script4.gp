{
set terminal wxt
set dgrid3d 30,30
dataFile1='slice1.dat'
dataFile2='slice2.dat'
dataFile3='slice3.dat'

set view 60,45
set hidden3d
set palette defined (0 "blue", 0.5 "yellow", 1 "red")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0

splot'< paste slice1.dat slice2.dat slice3.dat' u 1:2:4, 
}


