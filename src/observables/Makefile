# compilation flags for the system
LAPACKDIR=~/software/lapack-3.6.1
CFLAGS=  -g -pg    -I $(LAPACKDIR)/LAPACKE/include 
BUILD_DIR=../build

# compiler
CC=mpic++

# linking flags
LFLAGS= -L $(LAPACKDIR) -llapacke -llapack -lblas -lgfortran -lm  -pg

# l directory of libxml2 ( necessary for IO )
XMLDIR=/home/luca/software/libxml

all:	qmc


XML_FLAGS=`${XMLDIR}/bin/xml2-config --cflags`
XML_IDIR=${XMLDIR}/include/libxml2
XML_LIBS=`${XMLDIR}/bin/xml2-config --libs`

all:	qmc

qmc:		
	cd ..;make

accumulator:	accumulator.cpp testAccumulator.cpp ../ptools.cpp 
	$(CC) $(CFLAGS) $(LFLAGS) $(XML_LIBS) -I $(XML_IDIR) accumulator.cpp testAccumulator.cpp ../tools.cpp ../ptools.cpp   -o test

plan: optimizeTest.cpp optimizePlan.cpp optimizePlan.h
	$(CC) -I $(XML_IDIR) $(CFLAGS) $(LFLAGS)  ../xml-input.cpp optimizeTest.cpp optimizePlan.cpp -o test $(XML_LIBS)	
optimize: optimizeTest.cpp optimizeTools.cpp optimizeTools.h
	$(CC) $(CFLAGS)  optimizeTest.cpp optimizeTools.cpp -o test $(LFLAGS)
