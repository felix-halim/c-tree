CXX      = g++
ODIR     = ../bin
DEPFILE  = .depends
CXXFLAGS = -O3 -g -std=c++0x# -stdlib=libc++
# CXXFLAGS = -O0 -g -pg -std=c++0x# -stdlib=libc++
FLAGS    = -I.

all: \
	$(ODIR)/comb_noup \
	$(ODIR)/comb_lfhv \
	$(ODIR)/comb2_noup \
	$(ODIR)/comb2_lfhv \
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


$(ODIR)/ctree_noup: ctree.cc ctree.h random.h test_noup.h
	$(CXX) $(CXXFLAGS) $(FLAGS) -DNOUP -Wall -o $@ $<

$(ODIR)/ctree_lfhv: ctree.cc ctree.h random.h test_lfhv.h
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


depend:
	$(CXX) $(CXXFLAGS) -MM *.cc $(FLAGS) | sed -e 's@^\(.*\)\.o:@$(ODIR)/\1.o:@' > $(DEPFILE)

.PHONY: clean

clean:
	-@rm $(ODIR)/* 2> /dev/null || true

sinclude $(DEPFILE)
