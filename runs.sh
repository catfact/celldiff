#! /bin/bash
#i=0
#while [ "$i" -le 2 ]
#do
	# rm pp0_$i.dat
	./celldiff -x1 -n60 -c3000 -b1 -rrp30fp9.dat -ssp30fp9.dat -p0.3 -f.9 -t1000
  ./celldiff -x1 -n60 -c3000 -b1 -rrp30fp1.dat -ssp30fp1.dat -p0.3 -f.1 -t1000
#	i=`expr $i + 1`
#done
gnuplot plots.plt
