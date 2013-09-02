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
	make -s -C src "../bin/$2_noup" || exit;
	bin/$2_noup `hostname` $N $Q
	;;


lfhv)

	if [[ -n $3 ]]; then
		Q=$3
	fi
	make -s -C src "../bin/$2_lfhv" || exit;
	bin/$2_lfhv `hostname` $N $Q
	;;


batch_noup)

	./run.sh noup comb | tee -a $out
	./run.sh noup ctree | tee -a $out
	./run.sh noup ctree_uniform | tee -a $out
	./run.sh noup ctree_exp_leafsize | tee -a $out
	./run.sh noup ctree_eager | tee -a $out
	./run.sh noup comb2 | tee -a $out
	./run.sh noup btree_google | tee -a $out
	./run.sh noup sort | tee -a $out
	./run.sh noup btree_stx | tee -a $out
	;;


batch_lfhv)

	./run.sh lfhv comb | tee -a $out
	./run.sh lfhv ctree | tee -a $out
	./run.sh lfhv ctree_uniform | tee -a $out
	./run.sh lfhv comb2 | tee -a $out
	./run.sh lfhv btree_google | tee -a $out
	./run.sh lfhv btree_stx | tee -a $out
	;;


fhnet)

	scp index.html felixhalim@felix-halim.net:~/public_html/research/ctree
	scp results.csv felixhalim@felix-halim.net:~/public_html/research/ctree/results.csv
	;;

clean)

	rm -fr bin/*;;

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
