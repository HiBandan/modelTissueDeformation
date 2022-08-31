#include "CgMinimizerSearchDirectionUpdates.h"

double FletcherReeves::firstResetInterval = 1;
double FletcherReeves::subsequentResetInterval = 10;

double PolakRibiere::firstResetInterval = 0.5;
double PolakRibiere::subsequentResetInterval = 10;
