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
							./run.sh noup comb2 | tee -a $out
							./run.sh noup btree_google | tee -a $out
							./run.sh noup sort | tee -a $out
							./run.sh noup btree_stx | tee -a $out
							;;

batch_lfhv)		./run.sh lfhv comb | tee -a $out
							./run.sh lfhv ctree | tee -a $out
							./run.sh lfhv ctree_uniform | tee -a $out
							./run.sh lfhv comb2 | tee -a $out
							./run.sh lfhv btree_google | tee -a $out
							./run.sh lfhv btree_stx | tee -a $out
							;;

fhnet) 				scp index.html felixhalim@felix-halim.net:~/public_html/research/ctree
							scp results.csv felixhalim@felix-halim.net:~/public_html/research/ctree/results.csv
							;;

clean)				rm -fr bin/*;;

esac
