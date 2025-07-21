reset
#set terminal windows
#set terminal aqua
set terminal x11
#set terminal qt


set style line 1 lt 6 linewidth   3 lc 1
set style line 2 lt 4 linewidth   3 lc 5
set style line 3 lt 8 linewidth   3 lc 2
set style line 4 lt 1 linewidth   3 lc 3
set style line 5 lt 2 linewidth   3 lc 5
set style line 6 lt 3 linewidth   3 lc 6
set style line 7 lt 3 linewidth   3 lc rgbcolor '#800000'

set border linewidth 2

set style arrow 1 nohead ls 5

set xzeroaxis
set yzeroaxis


#set pm3d
#set palette rgbformulae 33,13,10
#set palette rgbformulae 21,22,23

#set title ' 16 X 16'
#set xrange [-1.0:1.0]
#set yrange [-1.0:1.0]

#set xlabel '{/Helvetica-Oblique P} (MPa)' enhanced font ',28'

set xlabel 'x' enhanced font 'Helvetica-Oblique,18'
set ylabel 'A' enhanced font 'Helvetica-Oblique,18'

#set ylabel 'Mass (tons)' enhanced font ',28'
#set ylabel 'Subsidence (ft)' enhanced font ',28'
#set ylabel 'Mass Rate (tons/day)' enhanced font ',28' rotate by 0 center


#set key center top
set key right top
#set key left top  Left
#set key right top  Right
set key box

set size .7,1 # grafico quadrado


#set data style lines
# malha 11X11 nos ou 10X10 volumes
#set dgrid3d 11, 11, 2
#set dgrid3d
set hidden3d
#set cntrparam levels 10
#set contour


set title 'Metodos de Volumes Finitos' font ',26'

hsuperior(x)= x>=xcentro-raio? x<=xcentro+raio? ycentro + sqrt(raio**2 - (x-xcentro)**2) : 0 : 0
hinferior(x) = x>=xcentro-raio? x<=xcentro+raio? ycentro - sqrt(raio**2 - (x-xcentro)**2) : 0 : 0


plot "data.100" using 1:2 with lines title 'Dados reais' linestyle 4
replot f(x) with lines title 'Ajuste' linestyle 7
#replot "a.out" using 1:4 with lines title 'zb' linestyle 7
#replot "fort.100" using 1:2 with lines title 'DG1' linestyle 4

pause -1

set terminal postscript eps enhanced
set termoption enhanced
set output "x-a-pb.eps"

replot

set terminal postscript eps enhanced color font ',20'
set termoption enhanced
set output "x-a.eps"

replot

set terminal png giant size 900,600 enhanced
set termoption enhanced
set output "x-a.png"

replot
