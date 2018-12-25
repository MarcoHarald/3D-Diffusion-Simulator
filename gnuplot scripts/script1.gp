{
set pm3d depth order

set terminal wxt
set dgrid3d 30,30
dataFile='slice1.dat'

set table dataFile1.'.grid'
splot dataFile u 1:2:3
unset table

set table dataFile.'.color'
splot dataFile u 1:2:4
unset table

set view 60,45
set hidden3d
set palette defined (0 "blue", 0.5 "yellow", 1 "red")
set autoscale cbfix
set pm3d
unset dgrid3d
set ticslevel 0
splot sprintf('< paste %s.grid %s.color', dataFile, dataFile) u 1:2:3:7 with pm3d notitle

