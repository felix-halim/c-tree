N=100000000
Q=1000000000
S=10

for (( i = 0; i < 10; i++ ))
do
  bin/comb_noup $N $Q
  echo ""
  sleep $S
  bin/ctree_noup $N $Q
  echo ""
  sleep $S
  bin/comb2_noup $N $Q
  echo ""
  sleep $S
  bin/btree_google_noup $N $Q
  echo ""
  sleep $S
  bin/sort_noup $N $Q
  echo ""
  sleep $S
  bin/btree_stx_noup $N $Q
  echo ""
  sleep $S
  echo ""
done
