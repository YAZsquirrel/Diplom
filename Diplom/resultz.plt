unset key
set palette functions gray, gray, gray 
set cntrparam linear
set cntrparam points 100
# smooth unique
set dgrid3d 200,200, 1
set dgrid3d splines 
# set dgrid3d 11,11
set pm3d at b
set ticslevel 0.8
set contour both
set samples 50,50
set pm3d map
#show pm3d
set size ratio 1
splot "ResultZ.txt" u 1:2:3 #, "ResultZ.txt" u 1:2:3:("") w labels point pt 1
unset pm3d
unset contour
unset dgrid3d