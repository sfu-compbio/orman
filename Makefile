# 786

CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio124
CPLEXINC=$(CPLEXDIR)/cplex/include
CPLEXLIB=$(CPLEXDIR)/cplex/lib/x86-64_sles10_4.1/static_pic
CONCERTINC=$(CPLEXDIR)/concert/include
CONCERTLIB=$(CPLEXDIR)/concert/lib/x86-64_sles10_4.1/static_pic

CC=g++
LDFLAGS=-lm -lpthread -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -I boost150 -std=gnu++0x #-DNDEBUG
SOURCES=stage1_parse.cc stage2_setcover.cc 
OBJECTS=common.o interval.o
CPLEXFLAGS=-DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) -L $(CPLEXLIB) -L $(CONCERTLIB) -Wl,--start-group -lconcert -lilocplex -lcplex -lpthread -lm -DNDEBUG

all:
	g++ $(LDFLAGS) -o stage1 common.cc stage1_parse.cc 
	g++ $(LDFLAGS) $(CPLEXFLAGS) -o stage2 common.cc stage2_setcover.cc 
	g++ $(LDFLAGS) -o stage3 common.cc stage3_interpret.cc 

simulator:
	g++ $(LDFLAGS) -o simulator common.cc simulator.cc

rescue:
	g++ $(LDFLAGS) -o rescue common.cc rescue.cc
#all: $(SOURCES) 
#
#$(EXECUTABLE): $(OBJECTS) 
#	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

#%.o: $(OBJECTS)
#	$(CC) $(CFLAGS) $< -o $@

#clean:
#	rm -f *.o 


