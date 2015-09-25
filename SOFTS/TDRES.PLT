set title "MIAMI -- Single Freq., TD Residuals"
set xlabel "Time Past 0h  (s)"
set ylabel "TD Residuals  (m)"
set nokey
# *** this is a comment line
set yrange [-0.04:0.04]
plot "tdxl.res" using 1:2 with points 8 2
