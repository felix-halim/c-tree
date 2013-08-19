mkdir bin 2> /dev/null
(cd src; make)

N=100000000
Q=1000000000
out=results.js

printf "var NOUP = [\n" > $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: $N, q: $q, " >> $out
	bin/comb_noup $N $q >> $out
	bin/ctree_noup $N $q >> $out
	bin/comb2_noup $N $q >> $out
	bin/sort_noup $N $q >> $out
	bin/btree_stx_noup $N $q >> $out
	bin/btree_google_noup $N $q >> $out
	printf " },\n" >> $out
done

printf "];\n\nvar LFHV = [\n" >> $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: $N, q: $q, " >> $out
	bin/comb_lfhv $N $q >> $out
	bin/ctree_lfhv $N $q >> $out
	bin/comb2_lfhv $N $q >> $out
	bin/btree_stx_lfhv $N $q >> $out
	bin/btree_google_lfhv $N $q >> $out
	printf " },\n" >> $out
done

printf "];\n" >> $out
