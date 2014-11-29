set terminal x11

set title "StandardModel renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='StandardModel_rgflow.dat'

plot for [i=2:33+1] filename using 1:(column(i)) title columnhead(i)
