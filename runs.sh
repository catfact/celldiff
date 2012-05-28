#! /bin/bash
i=0
while [ "$i" -le 9 ]
do
	rm pp0_$i.dat
	./celldiff -x1 -n64 -c2000 -b1.5 -rpp0_$i.dat -p0.$i
	i=`expr $i + 1`
done
gnuplot plots.plt
