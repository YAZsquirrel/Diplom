unset key
set palette functions gray, gray, gray 
set cntrparam linear
set dgrid 50,50
set pm3d at b
set ticslevel 0.8
set samples 1000,1000
set pm3d map
splot "ResultZ.txt" u 1:2:3
set size ratio 1

