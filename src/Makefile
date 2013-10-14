CXX      = g++
# CXX      = clang++-mp-3.4
ODIR     = ../bin
DEPFILE  = .depends
CXXFLAGS = -O3 -g -DNDEBUGx -std=c++0x# -stdlib=libc++
# CXXFLAGS = -O0 -g -pg -std=c++0x# -stdlib=libc++
FLAGS    = -I.

all: \
	$(ODIR)/art \
	$(ODIR)/art_best \
	$(ODIR)/art_crack \
	$(ODIR)/crack \
	$(ODIR)/comb \
	$(ODIR)/ctree \
	$(ODIR)/ctree_eager \
	$(ODIR)/ctree_with_crack \
	$(ODIR)/ctree_exp_leafsize \
	$(ODIR)/btree_google \
	$(ODIR)/btree_stx \
	$(ODIR)/sort \
	$(ODIR)/ctree_tests

$(ODIR)/dataset_gen: dataset_gen.cc
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

../data/10000000.data: $(ODIR)/dataset_gen
	mkdir -p ../data
	$(ODIR)/dataset_gen 10000000
	mv 10000000.data ../data

../data/100000000.data: $(ODIR)/dataset_gen
	mkdir -p ../data
	$(ODIR)/dataset_gen 100000000
	mv 100000000.data ../data

$(ODIR)/split: split.cc update.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DBSIZE=4096 -Wall -o $@ $<

$(ODIR)/comb: comb.cc comb.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/crack: crack.cc crack.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/art: art.cc art.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/arto: arto.cc arto.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/art_best: art_best.cc art_best.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/art_best_eager: art_best.cc art_best.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DEAGER -Wall -o $@ $<

$(ODIR)/art_crack: art_crack.cc art_crack.h art.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/artd_crack: artd_crack.cc artd_crack.h art.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_exp_leafsize: ctree_exp_leafsize.cc ctree_exp_leafsize.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_with_crack: ctree_with_crack.cc ctree_with_crack.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/sort: sort.cc test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/btree_google: btree_google.cc google/btree_set.h google/btree.h google/btree_container.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/btree_stx: btree_stx.cc stx/btree_multiset stx/btree_multiset.h stx/btree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_tests: ctree_tests.cc ctree.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_eager: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -DEAGER -Wall -o $@ $<

$(ODIR)/ctree_8_8: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_8_16: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_8_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_8_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_8_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_8_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_8_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_8_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_16_8: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_16_16: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_16_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_16_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_16_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_16_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_16_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_16_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_32_8: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_32_16: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_32_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_32_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_32_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_32_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_32_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_32_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_64_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_64_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_64_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_64_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_64_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_64_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_128_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_128_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_128_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_128_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_128_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_128_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_256_32: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_256_64: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_256_128: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_256_256: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_256_1024: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_256_4096: ctree.cc ctree.h test.h update.h query.h checksum.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=4096 -Wall -o $@ $<


depend:
	$(CXX) $(CXXFLAGS) -MM *.cc $(FLAGS) | sed -e 's@^\(.*\)\.o:@$(ODIR)/\1.o:@' > $(DEPFILE)

.PHONY: clean

clean:
	-@rm $(ODIR)/* 2> /dev/null || true

sinclude $(DEPFILE)
