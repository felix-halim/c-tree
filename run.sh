STD=c++0x
Q=100000000
out=results.js

mkdir bin 2> /dev/null

# Compile with "NOUP" will use test_noup.h
g++ -O3 -std=$STD -DNOUP -o bin/comb -Isrc src/comb.cc
g++ -O3 -std=$STD -DNOUP -o bin/sort -Isrc src/sort.cc
g++ -O3 -std=$STD -DNOUP -o bin/btree_stx -Isrc src/btree_stx.cc
g++ -O3 -std=$STD -DNOUP -o bin/btree_google -Isrc src/btree_google.cc

# Compile without "NOUP" will use test_lfhv.h
g++ -O3 -std=$STD -Isrc -o bin/update_comb src/comb.cc
g++ -O3 -std=$STD -Isrc -o bin/update_btree_stx src/btree_stx.cc
g++ -O3 -std=$STD -Isrc -o bin/update_btree_google src/btree_google.cc

printf "var NOUP = [\n" > $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: 100000000, q: $q, " >> $out
	bin/comb 100000000 $q >> $out
	bin/sort 100000000 $q >> $out
	bin/btree_stx 100000000 $q >> $out
	bin/btree_google 100000000 $q >> $out
	printf " },\n" >> $out
done

printf "];\n\nvar LFHV = [\n" >> $out

for (( q=1; q<=$Q; q*=10 ))
do
	printf "\t{ n: 100000000, q: $q, " >> $out
	bin/update_comb 100000000 $q >> $out
	bin/update_btree_stx 100000000 $q >> $out
	bin/update_btree_google 100000000 $q >> $out
	printf " },\n" >> $out
done

printf "];\n" >> $out
