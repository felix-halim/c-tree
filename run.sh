mkdir bin 2> /dev/null
N=100000000
Q=1
out=results.js

case $1 in

test) make -C src "../bin/$2_tests" && bin/$2_tests;;

noup)

	if [[ -n $3 ]]; then
		Q=$3
	fi
	make -s -C src "../bin/$2_noup"
	bin/$2_noup `hostname` $N $Q
	;;

batch)

	S=3
	for (( Q = 1; Q <= 1000000000; Q++ ))
	do
		./run.sh noup comb $Q | tee -a $out; sleep $S
		./run.sh noup ctree $Q | tee -a $out; sleep $S
		./run.sh noup comb2 $Q | tee -a $out; sleep $S
		./run.sh noup btree_google $Q | tee -a $out; sleep $S
		./run.sh noup sort $Q | tee -a $out; sleep $S
		./run.sh noup btree_stx $Q | tee -a $out; sleep $S
	done
	;;

bench)

	(cd src; make)

	N=100000000
	Q=1000000000

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
	;;

esac
