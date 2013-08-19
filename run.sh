# g++ -std=c++11 -Isrc -O3 src/comb_tests.cc && ./a.out
# g++ -std=c++11 -Isrc -O3 src/comb_tests.cc && ./a.out

STD=c++0x
N=100000000
Q=100000000
out=results.js

mkdir bin 2> /dev/null

# Compile with "NOUP" will use test_noup.h
g++ -O3 -std=$STD -DNOUP -o bin/comb -Isrc src/comb.cc
g++ -O3 -std=$STD -DNOUP -o bin/ctree -Isrc src/ctree.cc
g++ -O3 -std=$STD -DNOUP -o bin/comb2 -Isrc src/comb2.cc
g++ -O3 -std=$STD -DNOUP -o bin/sort -Isrc src/sort.cc
g++ -O3 -std=$STD -DNOUP -o bin/btree_stx -Isrc src/btree_stx.cc
g++ -O3 -std=$STD -DNOUP -o bin/btree_google -Isrc src/btree_google.cc

# Compile without "NOUP" will use test_lfhv.h
# g++ -O3 -std=$STD -Isrc -o bin/update_comb src/comb.cc
# g++ -O3 -std=$STD -Isrc -o bin/update_comb2 src/comb2.cc
# g++ -O3 -std=$STD -Isrc -o bin/update_btree_stx src/btree_stx.cc
# g++ -O3 -std=$STD -Isrc -o bin/update_btree_google src/btree_google.cc

printf "var NOUP = [\n" > $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: $N, q: $q, " >> $out
	bin/comb $N $q >> $out
	bin/ctree $N $q >> $out
	bin/comb2 $N $q >> $out
	bin/sort $N $q >> $out
	bin/btree_stx $N $q >> $out
	bin/btree_google $N $q >> $out
	printf " },\n" >> $out
done

printf "];\n\nvar LFHV = [\n" >> $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: $N, q: $q, " >> $out
	# bin/update_comb $N $q >> $out
	# bin/update_comb2 $N $q >> $out
	# bin/update_btree_stx $N $q >> $out
	# bin/update_btree_google $N $q >> $out
	printf " },\n" >> $out
done

printf "];\n" >> $out
