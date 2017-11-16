echo "Generating Makefile"
cat qmc.config src/MakefileTemplate > src/Makefile
echo "Building Build Directory..."
mkdir -p build
mkdir -p build/observables
mkdir -p build/particles
mkdir -p build/qmcDriver
echo "Done"
