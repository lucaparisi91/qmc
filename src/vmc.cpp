#include "dmc.h"
#include "vmc.h"
#include "qmc.h"
#include <fstream>
#include <cstdlib>
#include <cassert>
#include "input.h"
#include "tools.h"
#include "random.h"
#include "system.h"
#include "geometry.h"
#include "wavefunction.h"
#include "measures.h"
#include "xml-input.h"
#include <sstream>

using namespace std;
// perform the metropolis step on the single walker

