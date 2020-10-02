CXX = g++-9
CC = gcc-9
CXXFLAGS = -std=c++17 -lstdc++fs -O3
CFLAGS = -Ofast

.PHONY: clean

EXEs = bin/alg bin/alg2 bin/gen_fail bin/gen_input bin/libloadflow.a

all : $(EXEs)

debug: CXXFLAGS += -g
debug: CFLAGS += -g
debug: $(EXEs)

############### ALG - original DDLB algorithm ################
INCLUDEFLAGS2 = -I. -Iomnetpp/include
LDDFLAGS2 = -Lbin -Lomnetpp/lib -Wl,--whole-archive -loppcmdenv -Wl,--no-whole-archive -loppmain -loppenvir -loppenvir -loppsim  -loppnedxml -loppcommon -lloadflow -lxml2 -ldl

ALG_SRCDIR = src/original_ddlb_alg

bin/alg : $(ALG_SRCDIR)/distAlg.cc $(ALG_SRCDIR)/distAlg.h bin/libloadflow.a
	$(CXX) $(CXXFLAGS) $(ALG_SRCDIR)/distAlg.cc -o $@ $(INCLUDEFLAGS2) $(LDDFLAGS2)

############## ALG2 - DDLB algorithm with main  ##############
INCLUDEFLAGS = -Iomnetpp/include
LDDFLAGS = -Lbin -Lomnetpp/lib -loppsim -loppnedxml -loppcommon -lxml2 -lloadflow

ALG2_SRC := src/distAlg.cc src/main.cc
ALG2_HEADERS := src/distAlg.h

bin/alg2 : $(ALG2_SRC) $(ALG2_HEADERS) bin/libloadflow.a
	$(CXX) $(CXXFLAGS) $(ALG2_SRC) -o $@ $(INCLUDEFLAGS) $(LDDFLAGS) 
   
############## GENFAIL #################
bin/gen_fail: src/gen_fail.cpp src/bfs.cpp
	$(CXX) $(CXXFLAGS) $^ -lpthread -o $@

############## GENINPUTS ###############
bin/gen_input: src/gen_inps.cpp bin/libloadflow.a
	$(CXX) $(CXXFLAGS) $^ -lpthread -o $@

############## LOADFLOW ################
LOADFLOW_SRCDIR := src/loadflow
LOADFLOW_SRC_FILES := $(wildcard $(LOADFLOW_SRCDIR)/*.c)
LOADFLOW_OBJ_FILES := $(patsubst $(LOADFLOW_SRCDIR)/%.c,obj/%.o,$(LOADFLOW_SRC_FILES))

bin/libloadflow.a: $(LOADFLOW_OBJ_FILES)
	ar rcs $@ $^
	
obj/%.o: $(LOADFLOW_SRCDIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean : 
	rm -f $(EXEs)
	rm -f obj/*.o

	