# 786

CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio124
CPLEXINC=$(CPLEXDIR)/cplex/include
CPLEXLIB=$(CPLEXDIR)/cplex/lib/x86-64_sles10_4.1/static_pic
CONCERTINC=$(CPLEXDIR)/concert/include
CONCERTLIB=$(CPLEXDIR)/concert/lib/x86-64_sles10_4.1/static_pic
CPLEXFLAGS=-DIL_STD -I $(CPLEXINC) -I $(CONCERTINC) -L $(CPLEXLIB) -L $(CONCERTLIB) -Wl,--start-group -lconcert -lilocplex -lcplex -lpthread -lm -DNDEBUG

CC=g++

CFLAGS=-c -g -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -I boost150 -std=gnu++0x
LDFLAGS=-g -lm -lpthread 

OBJECTS=$(SOURCES:.cc=.o)
SOURCES=$(wildcard *.cc)
EXECUTABLE=uniqorn

all: $(SOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

patterns.o:
	ld -r -b binary -o patterns.o patterns.bin

clean:
	rm -f *.o $(EXECUTABLE)

#tools:
#	g++ $(LDFLAGS) -o simulator-mea common.cc simulator.cc



