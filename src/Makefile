CXX      = g++
# CXX      = clang++-mp-3.4
ODIR     = ../bin
DEPFILE  = .depends
CXXFLAGS = -O3 -g -DNDEBUGx -std=c++0x# -stdlib=libc++
# CXXFLAGS = -O0 -g -pg -std=c++0x# -stdlib=libc++
FLAGS    = -I.

all: \
	$(ODIR)/comb_noup \
	$(ODIR)/comb_lfhv \
	$(ODIR)/comb2_noup \
	$(ODIR)/comb2_lfhv \
	$(ODIR)/ctree_exp_leafsize_noup \
	$(ODIR)/ctree_exp_leafsize_lfhv \
	$(ODIR)/ctree_uniform_noup \
	$(ODIR)/ctree_uniform_lfhv \
	$(ODIR)/ctree_eager_noup \
	$(ODIR)/ctree_noup \
	$(ODIR)/ctree_lfhv \
	$(ODIR)/btree_google_noup \
	$(ODIR)/btree_google_lfhv \
	$(ODIR)/btree_stx_noup \
	$(ODIR)/btree_stx_lfhv \
	$(ODIR)/sort_noup \
	$(ODIR)/comb2_tests \
	$(ODIR)/ctree_tests

$(ODIR)/comb_noup: comb.cc comb.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/comb_lfhv: comb.cc comb.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/comb2_noup: comb2.cc comb2.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/comb2_lfhv: comb2.cc comb2.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/ctree_exp_leafsize_noup: ctree_exp_leafsize.cc ctree_exp_leafsize.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/ctree_exp_leafsize_lfhv: ctree_exp_leafsize.cc ctree_exp_leafsize.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_with_crack_noup: ctree_with_crack.cc ctree_with_crack.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/ctree_with_crack_lfhv: ctree_with_crack.cc ctree_with_crack.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/sort_noup: sort.cc random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/sort_lfhv: sort.cc random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/btree_google_noup: btree_google.cc google/btree_set.h google/btree.h google/btree_container.h test_noup.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/btree_google_lfhv: btree_google.cc google/btree_set.h google/btree.h google/btree_container.h test_lfhv.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/btree_stx_noup: btree_stx.cc stx/btree_multiset stx/btree_multiset.h stx/btree.h test_noup.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/btree_stx_lfhv: btree_stx.cc stx/btree_multiset stx/btree_multiset.h stx/btree.h test_lfhv.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/comb2_tests: comb2_tests.cc comb2.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_tests: ctree_tests.cc ctree.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<


$(ODIR)/ctree_eager_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=128 -DEAGER -Wall -o $@ $<

$(ODIR)/ctree_32_32_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_32_64_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_32_128_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_32_256_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_32_1024_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_32_4096_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_64_32_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_64_64_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_64_128_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_64_256_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_64_1024_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_64_4096_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_128_32_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_128_64_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_128_128_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_128_256_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_128_1024_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_128_4096_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_256_32_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_256_64_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_256_128_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_256_256_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_256_1024_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_256_4096_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=4096 -Wall -o $@ $<




$(ODIR)/ctree_32_32_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_32_64_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_32_128_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_32_256_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_32_1024_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_32_4096_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_64_32_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_64_64_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_64_128_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_64_256_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_64_1024_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_64_4096_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_128_32_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_128_64_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_128_128_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_128_256_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_128_1024_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_128_4096_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_256_32_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_256_64_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_256_128_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_256_256_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_256_1024_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_256_4096_lfhv: ctree.cc ctree.h random.h test_lfhv.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=4096 -Wall -o $@ $<



depend:
	$(CXX) $(CXXFLAGS) -MM *.cc $(FLAGS) | sed -e 's@^\(.*\)\.o:@$(ODIR)/\1.o:@' > $(DEPFILE)

.PHONY: clean

clean:
	-@rm $(ODIR)/* 2> /dev/null || true

sinclude $(DEPFILE)
