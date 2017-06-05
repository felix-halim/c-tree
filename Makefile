CXX = clang++
CXX_FLAGS = -std=c++1y -stdlib=libc++ -O2 -DNDEBUG \
-Wfatal-errors -Wall -Wextra \
-Wpedantic -Wshadow -Wconversion \
-Wno-unused-parameter \
-Wno-sign-conversion

BDIR = ./build
IDIR = ./inputs
ODIR = ./outputs
CPP = $(wildcard src/*.cc)
OBJ = $(CPP:%.cc=$(BDIR)/%.o)
DEP = $(OBJ:%.o=%.d)

$(BDIR)/% : $(BDIR)/src/%.o
	$(CXX) $(CXX_FLAGS) $^ -o $@

-include $(DEP)

$(BDIR)/%.o : %.cc
	@mkdir -p $(@D)
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@


$(ODIR)/art: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/art
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/art > $@

$(ODIR)/multiset: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/multiset
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/multiset > $@

$(ODIR)/comb: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/comb
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/comb > $@

$(ODIR)/readonly: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/readonly
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/readonly > $@


$(IDIR)/sorted_1000000_1000.gz: $(BDIR)/igen_sorted
	$(BDIR)/igen_sorted 1000000 1000 | gzip > $(IDIR)/sorted_1000000_1000.gz

$(IDIR)/sorted_100000000_1000.gz: $(BDIR)/igen_sorted
	$(BDIR)/igen_sorted 100000000 1000 | gzip > $(IDIR)/sorted_100000000_1000.gz

$(ODIR)/std_multiset_sorted_100000000_1000: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/std_multiset
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/std_multiset > $@

$(ODIR)/multiset_sorted_100000000_1000: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/multiset
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/multiset > $@

$(ODIR)/multiset_sorted_1000000_1000: $(IDIR)/sorted_1000000_1000.gz $(BDIR)/multiset
	gunzip -c $(IDIR)/sorted_1000000_1000.gz | $(BDIR)/multiset > $@

$(BDIR)/comb1600: $(SDIR)/comb.cc $(SDIR)/comb.h $(SDIR)/tester.h $(SDIR)/common.h
	$(CXX) -DNDEBUG -DCOMB1600 -o $@ $<

$(ODIR)/comb1600_sorted_100000000_1000: $(IDIR)/sorted_100000000_1000.gz $(BDIR)/comb1600
	gunzip -c $(IDIR)/sorted_100000000_1000.gz | $(BDIR)/comb1600 > $@


all: \
	$(ODIR)/art \
	$(ODIR)/art_best \
	$(ODIR)/art_crack \
	$(ODIR)/crack \
	$(ODIR)/comb \
	$(ODIR)/multiset \
	$(ODIR)/ctree \
	$(ODIR)/ctree_eager \
	$(ODIR)/ctree_with_crack \
	$(ODIR)/ctree_exp_leafsize \
	$(ODIR)/btree_google \
	$(ODIR)/btree_stx \
	$(ODIR)/sort \
	$(ODIR)/ctree_tests

$(ODIR)/dataset_gen: dataset_gen.cc
	$(CXX) $(CXXFLAGS) -Wall -o $@ $<

../data/1000000.data: $(ODIR)/dataset_gen
	mkdir -p ../data
	$(ODIR)/dataset_gen 1000000
	mv 1000000.data ../data

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

$(ODIR)/comb3: comb3.cc comb3.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DCOMB6400 -Wall -o $@ $<

$(ODIR)/comb_art_1000000000_1000000000: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_0_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=0 -DSKEW_Y=100 -Wall -o $@ $<


$(ODIR)/comb_art_skew_10_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_10_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_10_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_10_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_10_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_10_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=10 -DSKEW_Y=100 -Wall -o $@ $<


$(ODIR)/comb_art_skew_20_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_20_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_20_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_20_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_20_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_20_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=20 -DSKEW_Y=100 -Wall -o $@ $<


$(ODIR)/comb_art_skew_30_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_30_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_30_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_30_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_30_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_30_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=30 -DSKEW_Y=100 -Wall -o $@ $<


$(ODIR)/comb_art_skew_40_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_40_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_40_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_40_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_40_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_40_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=40 -DSKEW_Y=100 -Wall -o $@ $<


$(ODIR)/comb_art_skew_50_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=0 -Wall -o $@ $<

$(ODIR)/comb_art_skew_50_20: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=20 -Wall -o $@ $<

$(ODIR)/comb_art_skew_50_40: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=40 -Wall -o $@ $<

$(ODIR)/comb_art_skew_50_60: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=60 -Wall -o $@ $<

$(ODIR)/comb_art_skew_50_80: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=80 -Wall -o $@ $<

$(ODIR)/comb_art_skew_50_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000000000 -DSMALL_TOUCH=1000000000 -DSKEW_X=50 -DSKEW_Y=100 -Wall -o $@ $<



$(ODIR)/comb_art_0_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=0 -DSMALL_TOUCH=0 -Wall -o $@ $<

$(ODIR)/comb_art_1_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1 -DSMALL_TOUCH=0 -Wall -o $@ $<

$(ODIR)/comb_art_1_1: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1 -DSMALL_TOUCH=1 -Wall -o $@ $<

$(ODIR)/comb_art_1_10: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1 -DSMALL_TOUCH=10 -Wall -o $@ $<

$(ODIR)/comb_art_1_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1 -DSMALL_TOUCH=100 -Wall -o $@ $<

$(ODIR)/comb_art_1_1000: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1 -DSMALL_TOUCH=1000 -Wall -o $@ $<

$(ODIR)/comb_art_10_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=10 -DSMALL_TOUCH=0 -Wall -o $@ $<

$(ODIR)/comb_art_10_1: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=10 -DSMALL_TOUCH=1 -Wall -o $@ $<

$(ODIR)/comb_art_10_10: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=10 -DSMALL_TOUCH=10 -Wall -o $@ $<

$(ODIR)/comb_art_10_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=10 -DSMALL_TOUCH=100 -Wall -o $@ $<

$(ODIR)/comb_art_10_1000: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=10 -DSMALL_TOUCH=1000 -Wall -o $@ $<

$(ODIR)/comb_art_100_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=100 -DSMALL_TOUCH=0 -Wall -o $@ $<

$(ODIR)/comb_art_100_1: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=100 -DSMALL_TOUCH=1 -Wall -o $@ $<

$(ODIR)/comb_art_100_10: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=100 -DSMALL_TOUCH=10 -Wall -o $@ $<

$(ODIR)/comb_art_100_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=100 -DSMALL_TOUCH=100 -Wall -o $@ $<

$(ODIR)/comb_art_100_1000: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=100 -DSMALL_TOUCH=1000 -Wall -o $@ $<

$(ODIR)/comb_art_1000_0: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000 -DSMALL_TOUCH=0 -Wall -o $@ $<

$(ODIR)/comb_art_1000_1: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000 -DSMALL_TOUCH=1 -Wall -o $@ $<

$(ODIR)/comb_art_1000_10: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000 -DSMALL_TOUCH=10 -Wall -o $@ $<

$(ODIR)/comb_art_1000_100: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000 -DSMALL_TOUCH=100 -Wall -o $@ $<

$(ODIR)/comb_art_1000_1000: comb_art.cc comb_art.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DLARGE_TOUCH=1000 -DSMALL_TOUCH=1000 -Wall -o $@ $<

$(ODIR)/comb800: comb.cc comb.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DCOMB800 -Wall -o $@ $<

$(ODIR)/comb1600: comb.cc comb.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DCOMB1600 -Wall -o $@ $<

$(ODIR)/comb3200: comb.cc comb.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DCOMB3200 -Wall -o $@ $<

$(ODIR)/comb_count: comb.cc comb.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -DCOUNT_QUERY -o $@ $<

$(ODIR)/crack: crack.cc crack.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/crack_count: crack.cc crack.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -DCOUNT_QUERY -o $@ $<

$(ODIR)/mdd1r: mdd1r.cc crack.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/dd1r: dd1r.cc crack.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/sort_art: sort_art.cc art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/arto: arto.cc arto.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/art_best: art_best.cc art_best.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/art_best_eager: art_best.cc art_best.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DEAGER -Wall -o $@ $<

$(ODIR)/art_crack: art_crack.cc art_crack.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/artd_crack: artd_crack.cc artd_crack.h art.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_exp_leafsize: ctree_exp_leafsize.cc ctree_exp_leafsize.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_with_crack: ctree_with_crack.cc ctree_with_crack.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/sort: sort.cc test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/btree_google: btree_google.cc google/btree_set.h google/btree.h google/btree_container.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/btree_stx: btree_stx.cc stx/btree_multiset stx/btree_multiset.h stx/btree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_tests: ctree_tests.cc ctree.h random.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -Wall -o $@ $<

$(ODIR)/ctree_eager: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -DEAGER -Wall -o $@ $<

$(ODIR)/ctree_8_8: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_8_16: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_8_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_8_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_8_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_8_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_8_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_8_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=8 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_16_8: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_16_16: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_16_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_16_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_16_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_16_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_16_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_16_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=16 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_32_8: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=8 -Wall -o $@ $<
$(ODIR)/ctree_32_16: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=16 -Wall -o $@ $<
$(ODIR)/ctree_32_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_32_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_32_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_32_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_32_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_32_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_64_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_64_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_64_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_64_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_64_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_64_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=64 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_128_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_128_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_128_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_128_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_128_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_128_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=128 -DLEAF_BSIZE=4096 -Wall -o $@ $<

$(ODIR)/ctree_256_32: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=32 -Wall -o $@ $<
$(ODIR)/ctree_256_64: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=64 -Wall -o $@ $<
$(ODIR)/ctree_256_128: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=128 -Wall -o $@ $<
$(ODIR)/ctree_256_256: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=256 -Wall -o $@ $<
$(ODIR)/ctree_256_1024: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=1024 -Wall -o $@ $<
$(ODIR)/ctree_256_4096: ctree.cc ctree.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=256 -DLEAF_BSIZE=4096 -Wall -o $@ $<


$(ODIR)/ctree2: ctree2.cc ctree2.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -Wall -o $@ $<

$(ODIR)/ctree3: ctree3.cc ctree3.h test.h update.h query.h util.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DINTERNAL_BSIZE=32 -Wall -o $@ $<


.PRECIOUS: $(BDIR)/%.o
.PHONY: clean

clean:
	-@rm -fr $(BDIR)/*
