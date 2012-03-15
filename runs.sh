#! /bin/bash
i=1
while [ "$i" -le 9 ]
do
	rm pp0_$i.dat
	./celldiff -x1 -n64 -c20000 -b1.5 -rpp0_$i.dat -p0.$i
	i=`expr $i + 1`
done
gnuplot plots.plt