CXXFLAGS = -std=c++17 -O3
INCLUDEFLAGS = -I. -Ivendor/include -I../LoadFlowVerifier 
LDDFLAGS = -L../LoadFlowVerifier -Lvendor/linux_lib -Wl,--whole-archive -loppcmdenv -Wl,--no-whole-archive -loppmain -loppenvir -loppenvir -loppsim  -loppnedxml -loppcommon -lloadflow -lxml2 -ldl 

.PHONY: all

debug: CXXFLAGS = -g
debug: alg

all : alg

alg : distAlg.cc distAlg.ned distAlg.h
	$(CXX) $(CXXFLAGS)  distAlg.cc -o $@ $(INCLUDEFLAGS) $(LDDFLAGS) 
   
clean : 
	rm -f alg
