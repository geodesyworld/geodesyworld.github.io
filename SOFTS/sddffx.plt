set title "OHAT to STKR -- L1 Freq., DD Phase Residuals"
set xlabel "Time Past 0h  (s)"
set ylabel "DD Phase Residuals  (m)"
set nokey
### *** set xrange [43800.:50400.]
### *** set yrange [-0.03:0.03]
plot "sddffx.res" using 1:2 with points 8 2
