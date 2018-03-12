#!/bin/sh
rm cord.pdf
for src in 1 2
do
  for prb in 1 2 3
  do
    rm plot.p
    echo 'set autoscale' >> plot.p
    echo 'set logscale xy' >> plot.p
    echo 'unset label' >> plot.p
    echo 'set xtic auto' >> plot.p
    echo 'set ytic auto' >> plot.p
    echo 'set grid xtic' >> plot.p
    echo 'set grid ytic' >> plot.p
    echo 'set ylabel "Error (L2 norm)" enhanced' >> plot.p
    echo 'set xlabel "Spatial mesh (cm)" enhanced' >> plot.p
    echo 'set format y "10^{%L}"' >> plot.p
    echo 'set format x "10^{%L}"' >> plot.p
    echo 'set title "Problem '${prb}'"' >> plot.p
    echo 'set size square' >> plot.p
    echo 'set pointsize 0.5' >> plot.p
    echo 'set xr [1.0e-4:10]' >> plot.p
    echo 'set key outside' >> plot.p
    echo 'set style line 1 lt 1 lc rgb "red" lw 1' >> plot.p
    echo 'set style line 2 lt 1 lc rgb "green" lw 1' >> plot.p
    echo 'set style line 3 lt 1 lc rgb "blue" lw 1' >> plot.p
    echo 'set style line 4 lt 1 lc rgb "magenta" lw 1' >> plot.p
    echo 'plot    "p'${prb}'-'${src}'.p" using 1:2 title "Diamond Difference" with points pt 1, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:3 title "Step Characteristics" with points pt 2, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:4 title "Linear Discontinuous FEM" with points pt 3, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:5 title "Linear Characteristics" with points pt 4, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:6 title "2nd-order" with lines ls 1, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:7 title "2nd-order" with lines ls 2, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:8 every ::1::8 title "3rd-order" with lines ls 3, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:9 every ::1::6 title "4th-order" with lines ls 4' >> plot.p
    echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
    echo 'set output "plot-'${src}'-'${prb}'.pdf"' >> plot.p
    echo 'replot' >> plot.p
    echo 'set terminal x11' >> plot.p
    gnuplot plot.p
  done  
done
for src in 3 4
do
  rm plot.p
  echo 'set autoscale' >> plot.p
  echo 'set logscale xy' >> plot.p
  echo 'unset label' >> plot.p
  echo 'set xtic auto' >> plot.p
  echo 'set ytic auto' >> plot.p
  echo 'set grid xtic' >> plot.p
  echo 'set grid ytic' >> plot.p
  echo 'set ylabel "Error (L2 norm)" enhanced' >> plot.p
  echo 'set xlabel "Spatial mesh (cm)" enhanced' >> plot.p
  echo 'set format y "10^{%L}"' >> plot.p
  echo 'set format x "10^{%L}"' >> plot.p
  echo 'set title "Problem 1"' >> plot.p
  echo 'set size square' >> plot.p
  echo 'set pointsize 0.5' >> plot.p
  echo 'set xr [1.0e-4:10]' >> plot.p
  echo 'set key outside' >> plot.p
  echo 'set style line 1 lt 1 lc rgb "red" lw 1' >> plot.p
  echo 'set style line 2 lt 1 lc rgb "green" lw 1' >> plot.p
  echo 'set style line 3 lt 1 lc rgb "blue" lw 1' >> plot.p
  echo 'set style line 4 lt 1 lc rgb "magenta" lw 1' >> plot.p
  echo 'plot    "p1-'${src}'.p" using 1:2 title "Diamond Difference" with points pt 1, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:3 title "Step Characteristics" with points pt 2, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:4 title "Linear Discontinuous FEM" with points pt 3, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:5 title "Linear Characteristics" with points pt 4, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:6 title "2nd-order" with lines ls 1, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:7 title "2nd-order" with lines ls 2, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:8 every ::1::8 title "3rd-order" with lines ls 3, \' >> plot.p
  echo '        "p1-'${src}'.p" using 1:9 every ::1::6 title "4th-order" with lines ls 4' >> plot.p
  echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
  echo 'set output "plot-'${src}'-1.pdf"' >> plot.p
  echo 'replot' >> plot.p
  echo 'set terminal x11' >> plot.p
  gnuplot plot.p
done
pdfunite plot*.pdf cord.pdf
rm plot*.pdf
