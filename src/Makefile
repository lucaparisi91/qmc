# compilation flags for the system
LAPACKDIR=~/software/lapack-3.6.1
CFLAGS=  -g    -I $(LAPACKDIR)/LAPACKE/include
BUILD_DIR=../build

# compiler
CC=mpic++

# linking flags
LFLAGS= -L $(LAPACKDIR) -llapacke -llapack -lblas -lgfortran -lm  

# l directory of libxml2 ( necessary for IO )
XMLDIR=/home/luca/software/libxml

all:	qmc

SRCS= spline.cpp xml-input.cpp ranlxd.cpp timer.cpp tools.cpp ptools.cpp mesh.cpp random.cpp circular_vector.cpp input.cpp functions.cpp particles.cpp particles/gradientParticles.cpp qmcDriver/qProbability.cpp observables/accumulator.cpp   observables/correlatedEnergyDifference.cpp observables/optimizeTools.cpp observables/optimizePlan.cpp jastrow.cpp jastrow_collection.cpp geometry.cpp potential.cpp factory.cpp  wavefunction.cpp measures.cpp dmc.cpp vmc.cpp buildQMCDriver.cpp main.cpp qmc.cpp

OBJS=$(SRCS:%.cpp=$(BUILD_DIR)/%.o)

XML_FLAGS=`${XMLDIR}/bin/xml2-config --cflags`
XML_IDIR=${XMLDIR}/include/libxml2
XML_LIBS=`${XMLDIR}/bin/xml2-config --libs`

qmc:	$(OBJS)
	$(CC) $^ $(LFLAGS) $(XML_LIBS)   -o ../qmc

$(BUILD_DIR)/%.o: 	%.cpp 	$(BUILD_DIR)/%.d
	$(CC) -c $(CFLAGS) -I $(XML_IDIR) $< -o $@

$(BUILD_DIR)/%.d: 	%.cpp
	@set -e; rm -f $@; \
	$(CC) -M $(CFLAGS) -I $(XML_IDIR) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	find $(BUILD_DIR) -type f -name '*.o' -exec rm {} \;
	find $(BUILD_DIR) -type f -name '*.d' -exec rm {} \;
	find $(BUILD_DIR) -type f -name '*.d.*' -exec rm {} \;
	rm -f ../qmc

-include $(SRCS:%.cpp=$(BUILD_DIR)/%.d)
