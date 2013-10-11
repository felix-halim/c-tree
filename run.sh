mkdir bin 2> /dev/null
N=data/100000000.data
Q=1000000000
out=results.js

case $1 in

test)         make -C src "../bin/$2_tests" && bin/$2_tests;;

compile)      make -s -C src "../bin/$2" || exit;
              make -s -C src "../$N" || exit;;

algo)         ./run.sh compile $2
              bin/$2 $N $Q 0.1 1 $3;;

sky_append)   ./run.sh compile $2
              bin/$2 ../scrack/data/skyserver.data $Q 0.1 1 6;;

skew)         make -s -C src "../bin/$2" || exit;
              bin/$2 `hostname` $N 2 $Q $3;;

data)         echo "var data = 'timestamp,algorithm,query_workload,update_workload,\
N,Q,selectivity,verified,insert_time,update_time,query_time,checksum,\
n_leaves,n_capacity,n_internals,max_depth,slack,in_size,ln_size,ia_free,ia_size,la_free,la_size\n\\" > data.js;
              sed -e 's/$/\\n\\/' results.js >> data.js;
              echo "';" >> data.js
              ;;

batch_noup)   ./run.sh algo comb 0 | tee -a $out
              ./run.sh algo ctree_32_64 0 | tee -a $out
              ./run.sh algo art 0 | tee -a $out
              ./run.sh algo arto 0 | tee -a $out
              # ./run.sh algo art_crack 0 | tee -a $out
              # ./run.sh algo ctree_exp_leafsize 0 | tee -a $out
              ./run.sh algo ctree_eager 0 | tee -a $out
              ./run.sh algo sort 0 | tee -a $out
              ./run.sh algo btree_google 0 | tee -a $out
              ./run.sh algo btree_stx 0 | tee -a $out
              ;;

batch_lfhv)   ./run.sh algo comb 1 | tee -a $out
              ./run.sh algo ctree_32_64 1 | tee -a $out
              ./run.sh algo art 1 | tee -a $out
              ./run.sh algo art_crack 1 | tee -a $out
              ./run.sh algo btree_google 1 | tee -a $out
              ./run.sh algo btree_stx 1 | tee -a $out
              ;;

batch_sky)    ./run.sh append comb | tee -a $out
              ./run.sh append ctree_32_64 | tee -a $out
              ./run.sh append btree_google | tee -a $out
              ./run.sh append btree_stx | tee -a $out
              ./run.sh append art | tee -a $out
              ./run.sh append arto | tee -a $out
              ./run.sh append art_best | tee -a $out
              ;;

art_batch)    for U in {0..5}
              do
                     ./run.sh algo art $U | tee -a $out
                     ./run.sh algo art_best $U | tee -a $out
                     ./run.sh algo art_best_eager $U | tee -a $out
              done
              ./run.sh sky_append art | tee -a $out
              ./run.sh sky_append art_best | tee -a $out
              ./run.sh sky_append art_best_eager | tee -a $out;;

batch)        ./run.sh batch_noup
              for U in {1..5}
              do
                     ./run.sh algo comb $U | tee -a $out
                     ./run.sh algo ctree_32_64 $U | tee -a $out
                     ./run.sh algo art $U | tee -a $out
                     # ./run.sh algo art_crack $U | tee -a $out
                     ./run.sh algo arto $U | tee -a $out
                     ./run.sh algo btree_google $U | tee -a $out
                     ./run.sh algo btree_stx $U | tee -a $out
              done
              ./run.sh batch_sky
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
