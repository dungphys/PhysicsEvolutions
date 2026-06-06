#ifndef INIT_H
#define INIT_H

#include "utils.h"

void heaviside(Lattice &lattice);
int exponential(Lattice &lattice, double slope);
int logistic(Lattice &lattice, double slope);

#endif