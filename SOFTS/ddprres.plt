set title "MIAMI -- Single Freq., DD PR Residuals"
set xlabel "Time Past 0h  (s)"
set ylabel "DD PR Residuals  (m)"
set nokey
#
# *** next line commented out
# *** set yrange [-10.0:10.0]
#
plot "ddpr1.res" using 1:2 with points 8 2
