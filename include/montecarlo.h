#include "atomlib.h"
#include <stdlib.h>

/*
 * Struct:  MCConfiguration 
 * --------------------
 * Defines a type for the monte carlo simulation
 *
 *  name: string containing the name of the calculation
 *  input: if any, the initial structure provided by the user
 *  hasInitial: 1 if the user provided a xyz structure, otherwise a random structure will be generated
 * 
 *  atoms: list of atoms
 *  natoms: number of atoms
 * 
 *  type: method of calculation for the potential. 0: Lennard-Jones, 1: Stillinger
 *  l: lenght of the periodic box
 *  density: the density of the liquid
 *  mass: the atomic mass of the atom
 * 
 *  epsilon: well depth for the atom
 *  sigma: radius of the "hard" core of the atom
 *  cutoff: the cutoff distance for computing the potential
 *  step: max distance that the atom can move in a monte carlo trial step
 * 
 *  iterations: max number opf montecarlo iterations
 *  printeach: print atoms each n iterations
 *  out: the output file (with xyz extension)
 */
typedef struct MCConfiguration {
    char *name;
    int hasInitial;

    int natoms ;
    atom *atoms ;

    int type;
    double l;
    double density;
    double temp;

    double sigma;
    double squared_cutoff;
    double step;

    int useNei;
    double squared_rskin;

    int equilib;
    int production;
} mcconfig;

typedef struct RadialDistribution {
    int naccum;
    int natoms;
    int *n;
    double dens;

    int resolution;
    double rmax;
    double rmax_squared;
    double deltaR;
} grdist;

/*
 * Function:  readMCFromFile 
 * --------------------
 * Reads the configuration for the simulation from a file
 *
 *  name: string containing the name of the calculation
 *  atoms: list of atoms
 *  natoms: number of atoms
 *  epsilon: well depth for the atom
 *  sigma: radius of the "hard" core of the atom
 *  type: method of calculation for the potential. 0: Lennard-Jones, 1: Stillinger
 */
mcconfig readMCFromFile(char *url);


/// @brief Starts a monte carlo simulation with the given configuration
/// @param config The configuration of the simulation
void performMonteCarlo(mcconfig config);


