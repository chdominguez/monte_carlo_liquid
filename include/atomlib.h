#include "vectormath.h"

// Boltzmann constant in kJ/(mol·K)
//double KB = 8.314462618e-3;

/*
 * Struct:  Atom 
 * --------------------
 * Defines a type for the atoms
 *
 *  name: string containing the name of the atom i.e. hydrogen
 *  element: int that defines the element number i.e. 1 for hydrogen
 *  vector: vector containing the position of the atom
 */
typedef struct Atom {
    char *name;
    int element;
    struct Vector position;

    int nei_capacity;
    int nei_count;
    int *nei;

} atom;

/*
 * Function:  getElementFromName 
 * --------------------
 * Returns the element number of a given atom string
 *
 *  name: string containing the name of the atom
 *
 *  returns: the number corresponging to the atom element
 */
int getElementFromName(char *name);

/*
 * Function:  getElementProperties
 * --------------------
 * Gets the epsilon (well depth) and the sigma (radius of the hard core) for a given element
 *
 *  element: int corresponding to the element number
 *  epsilon: pointer to the epsilon (double) vaule
 *  sigma: pointer to the sigma (double) value
 * 
 * Note: Only two elements (Si and Ar) are allowed.
 */
void getElementProperties(int element, double *epsilon, double *sigma);

/*
 * Function:  getElementProperties
 * --------------------
 * Gets the epsilon (well depth) and the sigma (radius of the hard core) for a given element
 *
 *  element: int corresponding to the element number
 *  epsilon: pointer to the epsilon (double) vaule
 *  sigma: pointer to the sigma (double) value
 * 
 * Note: Only two elements (Si and Ar) are allowed.
 */
void reducedVolume(int element, double *epsilon, double *sigma);

/*
 * Function:  getElementProperties
 * --------------------
 * Gets the epsilon (well depth) and the sigma (radius of the hard core) for a given element
 *
 *  element: int corresponding to the element number
 *  epsilon: pointer to the epsilon (double) vaule
 *  sigma: pointer to the sigma (double) value
 * 
 * Note: Only two elements (Si and Ar) are allowed.
 */
void reducedVolume(int element, double *epsilon, double *sigma);

/*
 * Function:  shortestAtomDistance 
 * --------------------
 * Returns the shortest distance between two atoms taking into account the minimum-image convention
 *
 *  a: first atom
 *  b: second atom
 *  l: the lenght of the box
 *
 *  returns: the shortest distance between the atoms.
 */
double shortestAtomDistance(atom a, atom b, double l, double cutoff);

/*
 * Function:  shortestMIDistance 
 * --------------------
 * Returns the shortest distance between two points taking into account the minimum-image convention
 *
 *  a: first position
 *  b: second position
 *  l: the lenght of the box
 *
 *  returns: the shortest distance between the positions.
 */
double shortestMIDistance(vector a, vector b, double l, double cutoff);

/*
 * Function:  shortestDistance 
 * --------------------
 * computes the longest distance between a list of atoms
 *
 *  n: the number of atoms in the list
 *  atoms: pointer to list of atoms
 *
 *  returns: the longest distance between the atoms.
 */
double longestDistance(int n, atom* atoms);

/*
 * Function:  centerAtomsToBox 
 * --------------------
 * Moves the atoms to all having positive values in the components and within a box of origin 0,0,0
 *
 *  atoms: array of atoms
 *  n: atoms in the list
 */
void centerAtomsToBox(atom *atoms, int n);

/*
 * Function:  boxSize 
 * --------------------
 * Calculates the box size for a 0,0,0 centered box
 *
 *  atoms: array of atoms with >= 0 positions
 *  n: atoms in the list
 *
 *  returns: the box size (l)
 */
double boxSize(atom *atoms, int n);

/// @brief Copy an atom instance into a new one
/// @param from The atom to copy
/// @param to Destination
void copyAtom(atom from, atom *to);

/// @brief Copies an atom list to a new one
/// @param from The atom list to copy
/// @param to Destination
/// @param n Array length
void copyAtomList(atom *from, atom **to, int n);

/// @brief Moves a random atom from the given list
/// @param atoms The list of atoms 
/// @param n The length of the list
/// @param l The length of the box for the mirror image convention
/// @param stepSize The size of the step
/// @return The moved index
int moveRandomAtom(atom **atoms, int n, double l, double stepSize);

/*
 * Function:  lennardJones 
 * --------------------
 * Compute the complete Lennard-Jones potential energy between an array of atoms.
 *
 *  atoms: list of atoms
 *  natoms: atoms on the list
 *  l: lenght of the box for the minimum image convention
 *  cutoff: cutoff distance
 *
 *  returns: the potential energy value for this particular arrangement of atoms in reduced units.
 */
double fullLennardJones(atom *atoms, int natoms, double l, double cutoff);

/*
 * Function:  singleLennardJones 
 * --------------------
 * Compute the Lennard-Jones potential for a single atom.
 *
 *  atoms: list of atoms
 *  natoms: atoms on the list
 *  n: the atom to compute
 *  l: lenght of the box for the minimum image convention
 *  cutoff: cutoff distance
 *
 *  returns: the potential energy value for this particular atom in reduced units.
 */
double singleLennardJones(atom *atoms, int natoms, int n, int m, double l, double cutoff);

double fullNeiLennardJones(atom *atoms, int natoms, double l, double cutoff);

double singleNeiLennardJones(atom *atoms, int natoms, int n, double l, double cutoff);

/// @brief Updates the atom neighbours for a list of atoms
/// @param atomlist The list of atoms to compute its neighbours
/// @param natoms The number of atoms in the list
/// @param l The box lenght for the MI convention
/// @param cutoff The cutoff for the neighbour distance
void updateNeighbours(atom **atomlist, int natoms, double l, double cutoff);

/// @brief Initializes an atom instance
/// @param element The element number
/// @param position The vector with the atom position
/// @return The initialized atom
atom atomInit(int element, vector position);