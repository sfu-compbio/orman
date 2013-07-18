# 786

#CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio124
#CPLEXINC=$(CPLEXDIR)/cplex/include
#CPLEXLIB=$(CPLEXDIR)/cplex/lib/x86-64_sles10_4.1/static_pic
#CONCERTINC=$(CPLEXDIR)/concert/include
#CONCERTLIB=$(CPLEXDIR)/concert/lib/x86-64_sles10_4.1/static_pic
#CPLEXFLAGS=-DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) -L $(CPLEXLIB) -L $(CONCERTLIB) -Wl,--start-group -lconcert -lilocplex -lcplex -lpthread -lm -DNDEBUG
DATE=$(shell date)
CF?=g++
DF?=-O3
LF?=
CFLAGS=-c $(DF) -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -std=gnu++0x -DCOMPILE_TIME='"$(DATE) with $(DF) $(LF)"' -I boost_1_53_0
LDFLAGS= $(LF) -lm -lpthread 

OBJECTS=$(SOURCES:.cc=.o)
SOURCES=$(wildcard *.cc)
EXECUTABLE=orman

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CF) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CF) $(CFLAGS) $< -o $@

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin

clean:
	rm -f *.o $(EXECUTABLE)

#tools:
#	g++ $(LDFLAGS) -o simulator-mea common.cc simulator.cc



