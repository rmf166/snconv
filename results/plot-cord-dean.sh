#!/bin/sh
rm cord.pdf
for src in 1 1
do
  for prb in 4
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
    if [ ${src} == "1" ] && [ ${prb} == "4" ]; then
    echo 'set xr [1.0e-6:1.0e+1]' >> plot.p
    else
    echo 'set xr [1.0e-4:10]' >> plot.p
    fi
    echo 'set key top left' >> plot.p
    echo 'set key width -1' >> plot.p
    echo 'set key spacing 0.75' >> plot.p
    echo 'set key font ",8"' >> plot.p
    echo 'set style line 1 lt 1 lc rgb "red" lw 1' >> plot.p
    echo 'set style line 2 lt 1 lc rgb "green" lw 1' >> plot.p
    echo 'set style line 3 lt 1 lc rgb "blue" lw 1' >> plot.p
    echo 'set style line 4 lt 1 lc rgb "magenta" lw 1' >> plot.p
    echo 'plot    "p'${prb}'-'${src}'.p" using 1:2 title "DD" with points pt 1, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:3 title "SC" with points pt 2, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:4 title "LD" with points pt 3, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:5 title "LC" with points pt 4, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:6 title "O(h^2)" with lines ls 1, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:7 title "O(h^2)" with lines ls 2, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:8 title "O(h^3)" with lines ls 3, \' >> plot.p
    echo '        "p'${prb}'-'${src}'.p" using 1:9 title "O(h^4)" with lines ls 4' >> plot.p
    echo 'set terminal pdfcairo enhanced color dashed' >> plot.p
    echo 'set output "plot-'${src}'-'${prb}'.pdf"' >> plot.p
    echo 'replot' >> plot.p
    echo 'set terminal x11' >> plot.p
    gnuplot plot.p
  done  
done
pdfunite plot*.pdf cord.pdf
rm plot*.pdf
