set terminal svg size 900, 500
set output "bench_spmv.svg"
set xlabel "Matrix dim A: M=N"
set ylabel "Runtime"
set title "SPMV"
set key outside
set pointsize 0.5
plot "hpc_bench_80.out" using 1:3 with linespoints lw 3 title "COO", \
     "hpc_bench_80.out" using 1:4 with linespoints lw 3 title "JDS", \
     "hpc_bench_80.out" using 1:5 with linespoints lw 3 title "SED", \
     "hpc_bench_80.out" using 1:6 with linespoints lw 3 title "SKY"
