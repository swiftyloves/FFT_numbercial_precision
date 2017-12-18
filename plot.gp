set output outputfile
set terminal png size 1024,980
set xlabel "x"
set ylabel "RMS"

if(!small){
plot for [i=0:0] 'error_data/error'.i using 1:2 title 'n = 4, x^7' with linespoints ls i
} else {
plot for [i=0:0] 'error_data/error_small'.i using 1:2 title 'x^('.i.')' with linespoints ls i
}

