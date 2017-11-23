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

orbitals:	 spinOrbital.cpp testOrbital.cpp spinOrbital.h 
	$(CC) $(CFLAGS) $(LFLAGS) $(XML_LIBS) -I $(XML_IDIR) testOrbital.cpp spinOrbital.cpp ../tools.cpp ../ptools.cpp   -o test
