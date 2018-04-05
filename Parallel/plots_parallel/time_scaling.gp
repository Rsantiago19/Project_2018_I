set ylabel 'time per task'
set xlabel '# processors'
set logscale x
set logscale y
set key outside
plot 'time_scaling.dat' u 1:($2/$1) w lp t 'total time'
replot 'time_scaling.dat' u 1:($3/$1) w lp t 'Vel_Verlet time'
replot 'time_scaling.dat' u 1:($4/$1) w lp t 'therm. time'
replot 'time_scaling.dat' u 1:($5/$1) w lp t 'forces time'
replot 'time_scaling.dat' u 1:($6/$1) w lp t 'moment time'
replot 'time_scaling.dat' u 1:($7/$1) w lp t 'kinetic time'
set term png
set output 'scaling.png'
replot
set output
set term x11












