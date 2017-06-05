#todo: measure transition time, skew workload.

mkdir inputs 2> /dev/null
mkdir outputs 2> /dev/null
N=data/10000000.data
Q=10000000
out=results.js

case $1 in

test)         make -C src "../bin/$2_tests" && bin/$2_tests;;

compile)      make -s -C src "../bin/$2" || exit;
              make -s -C src "../$N" || exit;;

algo)         ./run.sh compile $2 || exit;
              S=$3
              QW=$4
              UW=$5
              bin/$2 $N $Q $S $QW $UW | tee res/$2_${S}_${QW}_${UW}.csv;;
              # valgrind --leak-check=yes bin/$2 $N $Q 0.1 1 $3;;

skew)
              ./run.sh algo comb_art_skew_0_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_0_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_0_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_0_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_0_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_0_100 0.1 12 0 | tee -a $out

              ./run.sh algo comb_art_skew_10_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_10_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_10_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_10_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_10_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_10_100 0.1 12 0 | tee -a $out

              ./run.sh algo comb_art_skew_20_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_20_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_20_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_20_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_20_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_20_100 0.1 12 0 | tee -a $out

              ./run.sh algo comb_art_skew_30_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_30_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_30_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_30_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_30_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_30_100 0.1 12 0 | tee -a $out

              ./run.sh algo comb_art_skew_40_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_40_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_40_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_40_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_40_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_40_100 0.1 12 0 | tee -a $out

              ./run.sh algo comb_art_skew_50_0 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_50_20 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_50_40 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_50_60 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_50_80 0.1 12 0 | tee -a $out
              ./run.sh algo comb_art_skew_50_100 0.1 12 0 | tee -a $out

              ;;

domain)       ./run.sh compile $2 || exit;
              bin/$2 $N $Q 0.1 17 $3;;
              # valgrind --leak-check=yes bin/$2 $N $Q 0.1 1 $3;;

sky)          ./run.sh compile $2 || exit;
              bin/$2 data/skyserver.u5data $Q 0.000001 0 $3;;
              # valgrind --leak-check=yes bin/$2 $N $Q 0.1 1 $3;;

skies)
              ./run.sh sky btree_google 6 | tee -a $out
              ./run.sh sky btree_stx 6 | tee -a $out
              ./run.sh sky ctree_eager 6 | tee -a $out
              # ./run.sh sky comb_art_1000000000_1000000000 6 | tee -a $out
              # ./run.sh sky art_best_eager 6 | tee -a $out
              # ./run.sh sky mdd1r 6 | tee -a $out
              ./run.sh sky comb_art_0_0 6 | tee -a $out
              ./run.sh sky comb_art_1_0 6 | tee -a $out
              ./run.sh sky comb_art_10_10 6 | tee -a $out
              ./run.sh sky comb_art_1000_1000 6 | tee -a $out
              # ./run.sh sky art 6 | tee -a $out
              # ./run.sh sky comb_art_0_0 0 | tee -a $out
              # ./run.sh sky comb_art_1_0 0 | tee -a $out
              # ./run.sh sky comb_art_10_0 0 | tee -a $out
              # ./run.sh sky comb_art_10_1 0 | tee -a $out
              # ./run.sh sky comb_art_10_10 0 | tee -a $out
              # ./run.sh sky comb_art_100_0 0 | tee -a $out
              # ./run.sh sky comb_art_100_1 0 | tee -a $out
              # ./run.sh sky comb_art_100_10 0 | tee -a $out
              # ./run.sh sky comb_art_100_100 0 | tee -a $out
              # ./run.sh sky comb_art_1000_0 0 | tee -a $out
              # ./run.sh sky comb_art_1000_1 0 | tee -a $out
              # ./run.sh sky comb_art_1000_10 0 | tee -a $out
              # ./run.sh sky comb_art_1000_100 0 | tee -a $out
              # ./run.sh sky comb_art_1000_1000 0 | tee -a $out
              ;;

custom)       ./run.sh algo ctree_32_64 1
              ./run.sh algo art 1
              ./run.sh algo art_best_eager 1
              ./run.sh algo comb_art 1
              ;;

seq)          ./run.sh compile $2 || exit;
              bin/$2 $N $Q 0.000001 2 $3;;
              # valgrind --leak-check=yes bin/$2 $N $Q 0.001 1 $3;;

split)        make -s -C src "../bin/split" || exit;
              for ith in {0..100}
              do
                bin/split $N $2 $ith
              done;;

append)       ./run.sh compile $2 || exit;
              bin/$2 $N $Q $3 1 7;;
              # bin/$2 data/skyserver.data $Q 0.001 $3 6;;

append8)      ./run.sh compile $2 || exit;
              bin/$2 $N $Q $3 1 8;;
              # bin/$2 data/skyserver.data $Q 0.001 $3 6;;

data)         echo "var data = 'timestamp,algorithm,query_workload,update_workload,\
N,Q,selectivity,verified,insert_time,update_time,query_time,checksum,\
n_leaves,n_capacity,n_internals,max_depth,slack,in_size,ln_size,ia_free,ia_size,la_free,la_size\n\\" > data.js;
              sed -e 's/$/\\n\\/' results.js >> data.js;
              echo "';" >> data.js
              ;;

batch_sela)
              ./run.sh batch_sel mdd1r 0
              ./run.sh batch_sel comb_art_10_10 0
              ./run.sh batch_sel art_best_eager 0
              ./run.sh batch_sel mdd1r 1
              ./run.sh batch_sel comb_art_10_10 1
              ./run.sh batch_sel art_best_eager 1
              ;;

batch_seq)
              ./run.sh seq sort 0
              ./run.sh seq mdd1r 0
              ./run.sh seq btree_google 0
              ./run.sh seq comb_art_1000000000_1000000000 0
              ./run.sh seq art_best_eager 0

              ./run.sh seq mdd1r 1
              ./run.sh seq btree_google 1
              ./run.sh seq comb_art_1000000000_1000000000 1
              ./run.sh seq art_best_eager 1
              ;;

batch_sel) 
              ./run.sh algo $2 0.000001 $3 | tee -a $out
              ./run.sh algo $2 0.00001 $3 | tee -a $out
              ./run.sh algo $2 0.0001 $3 | tee -a $out
              ./run.sh algo $2 0.001 $3 | tee -a $out
              ./run.sh algo $2 0.01 $3 | tee -a $out
              ./run.sh algo $2 0.1 $3 | tee -a $out
              ./run.sh algo $2 0.5 $3 | tee -a $out
              ./run.sh algo $2 0.9 $3 | tee -a $out
              ;;

batch_noup)   ./run.sh algo comb_art_1000000000_1000000000 0.1 0 | tee -a $out
              ./run.sh algo sort 0.1 0 | tee -a $out
              ./run.sh algo mdd1r 0.1 0 | tee -a $out
              ./run.sh algo art_best_eager 0.1 0 | tee -a $out
              ;;

batch_sky)
              ./run.sh algo comb  1 | tee -a $out
              ./run.sh algo mdd1r 0 | tee -a $out
              ./run.sh algo ctree_32_64 0 | tee -a $out
              ./run.sh algo ctree_32_1024 0 | tee -a $out
              ./run.sh algo ctree_32_4096 0 | tee -a $out
              ./run.sh algo ctree_eager 0 | tee -a $out
              ./run.sh algo btree_google 0 | tee -a $out
              ./run.sh algo btree_stx 0 | tee -a $out
              ./run.sh algo ctree_32_64 1 | tee -a $out
              ./run.sh algo ctree_32_1024 1 | tee -a $out
              ./run.sh algo ctree_32_4096 1 | tee -a $out
              ./run.sh algo ctree_eager 1 | tee -a $out
              ./run.sh algo btree_google 1 | tee -a $out
              ./run.sh algo btree_stx 1 | tee -a $out
              ./run.sh algo crack 1 | tee -a $out
              ./run.sh algo mdd1r 1 | tee -a $out
              ./run.sh algo crack 0 | tee -a $out
              # ./run.sh append art | tee -a $out
              # ./run.sh append arto | tee -a $out
              # ./run.sh append art_best | tee -a $out
              ;;

batch)        #./run.sh batch_noup
              for U in {1..5}
              do
                     #./run.sh algo ctree_32_1024 $U | tee -a $out
                     #./run.sh algo ctree_32_4096 $U | tee -a $out
                     #./run.sh algo comb $U | tee -a $out
                     ./run.sh algo crack $U | tee -a $out
                     ./run.sh algo ctree_eager $U | tee -a $out
                     #./run.sh algo ctree_32_64 $U | tee -a $out
                     ./run.sh algo art $U | tee -a $out
                     # ./run.sh algo art_crack $U | tee -a $out
                     ./run.sh algo art_best $U | tee -a $out
                     ./run.sh algo art_best_eager $U | tee -a $out
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

fhnet2)
              scp graphs.html linechart2.js results.csv jquery-2.0.3.min.js trimmer.html trimmer.js felixhalim@felix-halim.net:~/public_html/research/ctree/
              ;;

fhnet)        scp experiments.html felixhalim@felix-halim.net:~/public_html/research/ctree/
              scp illustration.js felixhalim@felix-halim.net:~/public_html/research/ctree/
              scp data.js felixhalim@felix-halim.net:~/public_html/research/ctree/
              scp demo.html felixhalim@felix-halim.net:~/public_html/research/ctree/
              scp linechart.js felixhalim@felix-halim.net:~/public_html/research/ctree/
              ;;

clean)        rm -fr bin/*;;

esac
