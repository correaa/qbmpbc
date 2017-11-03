#!/bin/bash
gnuplot -persist <<EOF
fit a*x+b "<grep -h econst $*" u 0:2 via a,b
fit c "<grep -h econst $*" u 0:2 via c
p "<grep -h econst $*" u 2 w l, a*x+b, c
print "Econst_avg=",c,"   dE/dt=",a
EOF
