#include "montecarlo.h"

int fileReader(char *url, mcconfig *config);

void printXYZFile(atom *atoms, int natoms, char* name, double sigma, FILE *file);