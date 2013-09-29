mkdir bin 2> /dev/null
N=data/100000000.data
Q=1000000000
out=results.js

case $1 in

test)         make -C src "../bin/$2_tests" && bin/$2_tests;;

algo)         make -s -C src "../bin/$2" || exit;
              make -s -C src "../$N" || exit;
              bin/$2 $N $Q 0.1 1 $3;;

skew)         make -s -C src "../bin/$2" || exit;
              bin/$2 `hostname` $N 2 $Q $3;;

data)         echo "var data = 'timestamp,hostname,algorithm,N,Q,selectivity\
insert_time,query_time,checksum,version,nLeaves,nCap,nInternals,\
max_depth,slack,in_size,ln_size,ia_free,ia_size,la_free,la_size\n\\" > data.js;
              sed -e 's/$/\\n\\/' results.js >> data.js;
              echo "';" >> data.js
              ;;

batch_noup)   ./run.sh noup comb | tee -a $out
              ./run.sh noup ctree_32_64 | tee -a $out
              ./run.sh noup art | tee -a $out
              ./run.sh noup art_crack | tee -a $out
              # ./run.sh noup ctree_exp_leafsize | tee -a $out
              ./run.sh noup ctree_eager | tee -a $out
              ./run.sh noup sort | tee -a $out
              ./run.sh noup btree_google | tee -a $out
              ./run.sh noup btree_stx | tee -a $out
              ;;

batch_lfhv)   ./run.sh lfhv comb | tee -a $out
              ./run.sh lfhv ctree_32_64 | tee -a $out
              ./run.sh lfhv art | tee -a $out
              ./run.sh lfhv art_crack | tee -a $out
              ./run.sh lfhv btree_google | tee -a $out
              ./run.sh lfhv btree_stx | tee -a $out
              ;;

batch_skew)
              for (( selectivity = 1; selectivity <= $N; selectivity *= 10))
              do
                # ./run.sh skew comb $selectivity | tee -a $out
                ./run.sh skew ctree_32_64 $selectivity | tee -a $out
                # ./run.sh skew btree_google $selectivity | tee -a $out
                # ./run.sh skew btree_stx $selectivity | tee -a $out
              done
              ;;

batch_ctree)  ./run.sh noup ctree_32_32 | tee -a $out
              ./run.sh noup ctree_32_64 | tee -a $out
              ./run.sh noup ctree_32_128 | tee -a $out
              ./run.sh noup ctree_32_256 | tee -a $out
              ./run.sh noup ctree_64_32 | tee -a $out
              ./run.sh noup ctree_64_64 | tee -a $out
              ./run.sh noup ctree_64_128 | tee -a $out
              ./run.sh noup ctree_64_256 | tee -a $out
              ./run.sh noup ctree_128_32 | tee -a $out
              ./run.sh noup ctree_128_64 | tee -a $out
              ./run.sh noup ctree_128_128 | tee -a $out
              ./run.sh noup ctree_128_256 | tee -a $out
              ./run.sh noup ctree_256_32 | tee -a $out
              ./run.sh noup ctree_256_64 | tee -a $out
              ./run.sh noup ctree_256_128 | tee -a $out
              ./run.sh noup ctree_256_256 | tee -a $out
              ./run.sh lfhv ctree_32_32 | tee -a $out
              ./run.sh lfhv ctree_32_64 | tee -a $out
              ./run.sh lfhv ctree_32_128 | tee -a $out
              ./run.sh lfhv ctree_32_256 | tee -a $out
              ./run.sh lfhv ctree_64_32 | tee -a $out
              ./run.sh lfhv ctree_64_64 | tee -a $out
              ./run.sh lfhv ctree_64_128 | tee -a $out
              ./run.sh lfhv ctree_64_256 | tee -a $out
              ./run.sh lfhv ctree_128_32 | tee -a $out
              ./run.sh lfhv ctree_128_64 | tee -a $out
              ./run.sh lfhv ctree_128_128 | tee -a $out
              ./run.sh lfhv ctree_128_256 | tee -a $out
              ./run.sh lfhv ctree_256_32 | tee -a $out
              ./run.sh lfhv ctree_256_64 | tee -a $out
              ./run.sh lfhv ctree_256_128 | tee -a $out
              ./run.sh lfhv ctree_256_256 | tee -a $out

              ./run.sh noup ctree_8_8 | tee -a $out
              ./run.sh noup ctree_8_16 | tee -a $out
              ./run.sh noup ctree_8_32 | tee -a $out
              ./run.sh noup ctree_8_64 | tee -a $out
              ./run.sh noup ctree_8_128 | tee -a $out
              ./run.sh noup ctree_8_256 | tee -a $out

              ./run.sh noup ctree_16_16 | tee -a $out
              ./run.sh noup ctree_16_32 | tee -a $out
              ./run.sh noup ctree_16_64 | tee -a $out
              ./run.sh noup ctree_16_128 | tee -a $out
              ./run.sh noup ctree_16_256 | tee -a $out

              ./run.sh lfhv ctree_8_8 | tee -a $out
              ./run.sh lfhv ctree_8_16 | tee -a $out
              ./run.sh lfhv ctree_8_32 | tee -a $out
              ./run.sh lfhv ctree_8_64 | tee -a $out
              ./run.sh lfhv ctree_8_128 | tee -a $out
              ./run.sh lfhv ctree_8_256 | tee -a $out

              ./run.sh lfhv ctree_16_16 | tee -a $out
              ./run.sh lfhv ctree_16_32 | tee -a $out
              ./run.sh lfhv ctree_16_64 | tee -a $out
              ./run.sh lfhv ctree_16_128 | tee -a $out
              ./run.sh lfhv ctree_16_256 | tee -a $out

              ;;

fhnet)        scp index.html felixhalim@felix-halim.net:~/public_html/research/ctree
              scp results.csv felixhalim@felix-halim.net:~/public_html/research/ctree/results.csv
              ;;

clean)        rm -fr bin/*;;

esac
