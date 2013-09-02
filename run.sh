mkdir bin 2> /dev/null
N=100000000
out=results.js

case $1 in

test)					make -C src "../bin/$2_tests" && bin/$2_tests;;

noup)					make -s -C src "../bin/$2_noup" || exit;
							bin/$2_noup `hostname` $N;;

lfhv)					make -s -C src "../bin/$2_lfhv" || exit;
							bin/$2_lfhv `hostname` $N;;

batch_noup)		./run.sh noup ctree_uniform | tee -a $out
							./run.sh noup comb | tee -a $out
							./run.sh noup ctree | tee -a $out
							./run.sh noup ctree_exp_leafsize | tee -a $out
							./run.sh noup ctree_eager | tee -a $out
							./run.sh noup btree_google | tee -a $out
							./run.sh noup sort | tee -a $out
							./run.sh noup btree_stx | tee -a $out
							;;

batch_lfhv)		#./run.sh lfhv comb | tee -a $out
							./run.sh lfhv ctree_64_128 | tee -a $out
							./run.sh lfhv btree_google | tee -a $out
							./run.sh lfhv btree_stx | tee -a $out
							;;

batch_ctree)	./run.sh noup ctree_32_32 | tee -a $out
							./run.sh noup ctree_32_64 | tee -a $out
							./run.sh noup ctree_32_128 | tee -a $out
							./run.sh noup ctree_32_256 | tee -a $out
							./run.sh noup ctree_32_1024 | tee -a $out
							./run.sh noup ctree_32_4096 | tee -a $out
							./run.sh noup ctree_64_32 | tee -a $out
							./run.sh noup ctree_64_64 | tee -a $out
							./run.sh noup ctree_64_128 | tee -a $out
							./run.sh noup ctree_64_256 | tee -a $out
							./run.sh noup ctree_64_1024 | tee -a $out
							./run.sh noup ctree_64_4096 | tee -a $out
							./run.sh noup ctree_128_32 | tee -a $out
							./run.sh noup ctree_128_64 | tee -a $out
							./run.sh noup ctree_128_128 | tee -a $out
							./run.sh noup ctree_128_256 | tee -a $out
							./run.sh noup ctree_128_1024 | tee -a $out
							./run.sh noup ctree_128_4096 | tee -a $out
							./run.sh noup ctree_256_32 | tee -a $out
							./run.sh noup ctree_256_64 | tee -a $out
							./run.sh noup ctree_256_128 | tee -a $out
							./run.sh noup ctree_256_256 | tee -a $out
							./run.sh noup ctree_256_1024 | tee -a $out
							./run.sh noup ctree_256_4096 | tee -a $out
							./run.sh lfhv ctree_32_32 | tee -a $out
							./run.sh lfhv ctree_32_64 | tee -a $out
							./run.sh lfhv ctree_32_128 | tee -a $out
							./run.sh lfhv ctree_32_256 | tee -a $out
							./run.sh lfhv ctree_32_1024 | tee -a $out
							./run.sh lfhv ctree_32_4096 | tee -a $out
							./run.sh lfhv ctree_64_32 | tee -a $out
							./run.sh lfhv ctree_64_64 | tee -a $out
							./run.sh lfhv ctree_64_128 | tee -a $out
							./run.sh lfhv ctree_64_256 | tee -a $out
							./run.sh lfhv ctree_64_1024 | tee -a $out
							./run.sh lfhv ctree_64_4096 | tee -a $out
							./run.sh lfhv ctree_128_32 | tee -a $out
							./run.sh lfhv ctree_128_64 | tee -a $out
							./run.sh lfhv ctree_128_128 | tee -a $out
							./run.sh lfhv ctree_128_256 | tee -a $out
							./run.sh lfhv ctree_128_1024 | tee -a $out
							./run.sh lfhv ctree_128_4096 | tee -a $out
							./run.sh lfhv ctree_256_32 | tee -a $out
							./run.sh lfhv ctree_256_64 | tee -a $out
							./run.sh lfhv ctree_256_128 | tee -a $out
							./run.sh lfhv ctree_256_256 | tee -a $out
							./run.sh lfhv ctree_256_1024 | tee -a $out
							./run.sh lfhv ctree_256_4096 | tee -a $out
							;;

fhnet) 				scp index.html felixhalim@felix-halim.net:~/public_html/research/ctree
							scp results.csv felixhalim@felix-halim.net:~/public_html/research/ctree/results.csv
							;;

clean)				rm -fr bin/*;;

esac
