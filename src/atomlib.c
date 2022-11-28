#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atomlib.h"

int getElementFromName(char *name)
{
    // Only two types of atoms are supported on this software, Silicon and Argon.
    if (strcmp(name, "Si") == 0)
    {
        return 14;
    }

    if (strcmp(name, "Ar") == 0)
    {
        return 18;
    }

    return 0; // Element not supported
}

char *getNameFromElement(int element)
{
    if (element == 14)
    {
        return "Si";
    }

    if (element == 18)
    {
        return "Ar";
    }

    return "X"; // Element not supported
}

double shortestAtomDistance(atom a, atom b, double l, double cutoff)
{
    return shortestMIDistance(a.position, b.position, l, cutoff);
}

double shortestMIDistance(vector a, vector b, double l, double cutoff)
{
    double rij = sumSquaredComponents(substractVector(b, a));
    double rmij = sumSquaredComponents(substractVector(b, addToComponents(l, a)));

    double max_distance = cutoff * cutoff;

    if (rij > max_distance && rmij > max_distance) // Cancel the computation of the square root because its outside the cutoff
    {
        return -1;
    }

    if (rij > rmij)
    {
        return sqrt(rmij);
    }
    else
    {
        return sqrt(rij);
    }
}

double longestDistance(int n, atom *atoms)
{
    double longest = 0;
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            double current = norm(substractVector(atoms[j].position, atoms[i].position));
            if (current > longest)
            {
                longest = current;
            }
        }
    }
    return longest;
}

void centerAtomsToBox(atom *atoms, int n)
{
    double longestX = 0; // Find the largest negative distance
    double longestY = 0; // Find the largest negative distance
    double longestZ = 0; // Find the largest negative distance
    for (int i = 0; i < n; i++)
    {
        double distX = atoms[i].position.x;
        if (distX < longestX)
        {
            longestX = distX;
        }
        double distY = atoms[i].position.y;
        if (distY < longestY)
        {
            longestY = distY;
        }
        double distZ = atoms[i].position.z;
        if (distZ < longestZ)
        {
            longestZ = distZ;
        }
    }

    for (int i = 0; i < n; i++)
    {
        atoms[i].position.x -= longestX;
        atoms[i].position.y -= longestY;
        atoms[i].position.z -= longestZ;
    }
}

double boxSize(atom *atoms, int n)
{
    double maxX = 0;
    double maxY = 0;
    double maxZ = 0;
    for (int i = 0; i < n; i++)
    {
        vector pos = atoms[i].position;
        if (pos.x > maxX)
        {
            maxX = pos.x;
        }
        if (pos.y > maxY)
        {
            maxY = pos.y;
        }
        if (pos.z > maxZ)
        {
            maxZ = pos.z;
        }
    }

    double lenght = maxX;

    if (lenght < maxY)
    {
        lenght = maxY;
    }
    if (lenght < maxZ)
    {
        lenght = maxZ;
    }

    return lenght;
}

void copyAtom(atom from, atom *to)
{
    atom newAtom;
    newAtom.element = from.element;
    newAtom.name = "Ar";
    //newAtom.name = (char *)calloc(3, sizeof(char)); // Allocating memory for the name
    //strcpy(newAtom.name, from.name);
    newAtom.nei_count = from.nei_count;
    newAtom.nei = from.nei;
    vector newVector;
    newVector.x = from.position.x;
    newVector.y = from.position.y;
    newVector.z = from.position.z;
    newAtom.position = newVector;
    *to = newAtom;
}

void copyAtomList(atom *from, atom **to, int n)
{
    atom *newlist = (atom *)calloc(n, sizeof(atom));
    for (int i = 0; i < n; i++)
    {
        copyAtom(from[i], &newlist[i]);
    }
    *to = newlist;
}

int moveRandomAtom(atom **atoms, int n, double l, double stepSize)
{
    int r = getRandomInt(0, n - 1);

    (*atoms)[r].position.x += getRandomDouble(-1, 1) * stepSize;
    (*atoms)[r].position.y += getRandomDouble(-1, 1) * stepSize;
    (*atoms)[r].position.z += getRandomDouble(-1, 1) * stepSize;

    // X
    if ((*atoms)[r].position.x > l)
    {
        (*atoms)[r].position.x -= l;
    }
    if ((*atoms)[r].position.x < 0)
    {
        (*atoms)[r].position.x += l;
    }

    // Y
    if ((*atoms)[r].position.y > l)
    {
        (*atoms)[r].position.y -= l;
    }
    if ((*atoms)[r].position.y < 0)
    {
        (*atoms)[r].position.y += l;
    }

    // Z
    if ((*atoms)[r].position.z > l)
    {
        (*atoms)[r].position.z -= l;
    }
    if ((*atoms)[r].position.z < 0)
    {
        (*atoms)[r].position.z += l;
    }

    return r;
}

double fullLennardJones(atom *atoms, int natoms, double l, double cutoff)
{
    double etot = 0;
    int iter = 0;
    for (int i = 0; i < natoms - 1; i++) // The last atom will have been counted by all the other ones
    {
        etot += singleLennardJones(atoms, natoms, i, i + 1, l, cutoff);
        iter++;
    }
    return etot;
};

double singleLennardJones(atom *atoms, int natoms, int i, int m, double l, double cutoff)
{
    double etot = 0;
    atom a = atoms[i];
    for (int k = m; k < natoms; k++)
    {
        atom b = atoms[k];
        double dist = shortestAtomDistance(a, b, l, cutoff);
        if (dist > 0) // Skip because of cutoff (dist would be -1 // should not have bc of neighbour list)
        {
            double term = 1.0 / dist;
            double fterm = pow(term, 12);
            double sterm = pow(term, 6);
            double enow = fterm - sterm;
            etot += enow;
        }
    }
    return 4 * etot;
};

double singleNeiLennardJones(atom *atoms, int natoms, int i, double l, double cutoff)
{
    double etot = 0;
    atom a = atoms[i];
    for (int k = 0; k < atoms[i].nei_count; k++)
    {
        atom b = atoms[a.nei[k]];
        double dist = shortestAtomDistance(a, b, l, cutoff);
        if (dist > 0) // Skip because of cutoff (dist would be -1 // should not have bc of neighbour list)
        {
            double term = 1.0 / dist;
            double fterm = pow(term, 12);
            double sterm = pow(term, 6);
            double enow = fterm - sterm;
            etot += enow;
        }
    }
    return 4 * etot;
};


void addNeighbour(atom *a, int neighbour)
{
    a->nei_count++;
    if (a->nei_count == a->nei_capacity)
    {
        a->nei_capacity *= 2;
        a->nei = (int *)realloc(a->nei, a->nei_capacity*sizeof(int));
    }
    a->nei[a->nei_count - 1] = neighbour;
}

void updateNeighbours(atom **atomlist, int natoms, double l, double rskin)
{
    // Empty the current list
    for (int i = 0; i < natoms; i++)
    {
        free((*atomlist)[i].nei);
        (*atomlist)[i].nei_count = 0;
        (*atomlist)[i].nei_capacity = 100;
        (*atomlist)[i].nei = (int *)calloc(100, sizeof(int));
    }

    // Recalculate neighbours
    for (int i = 0; i < natoms; i++)
    {
        for (int j = i + 1; j < natoms; j++)
        {
            double dist = shortestAtomDistance((*atomlist)[i], (*atomlist)[j], l, rskin);
            if (dist > 0) // The rskin is already accounted in the distance function, (-1 if outside the rskin)
            {
                addNeighbour(&(*atomlist)[i], j);
                addNeighbour(&(*atomlist)[j], i);
            }
        }
    }
}

atom atomInit(int element, vector position)
{
    // Define the atom
    atom a;
    a.name = (char *)calloc(3, sizeof(char)); // Allocating memory for the name

    strcpy(a.name, getNameFromElement(element));
    a.element = element;

    // Start with no neighbours
    a.nei_capacity = 100;
    a.nei_count = 0;
    a.nei = (int *)calloc(100, sizeof(int));

    // Position of the atom
    a.position = position;

    return a;
}