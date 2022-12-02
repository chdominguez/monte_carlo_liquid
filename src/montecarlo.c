#include <stdio.h>
#include <string.h>
#include "tools.h"
#include <time.h>

double monteCarloProbability(double change, double reducedT)
{
    return exp(-(change / reducedT));
}

mcconfig readMCFromFile(char *url)
{
    mcconfig config;
    int i = fileReader(url, &config);
    if (i == 1)
    {
        exit(1);
    }
    return config;
}

/// @brief Generates a distribution of atoms based on the density and the box lenght
/// @param mass The mass of the atom
/// @param density The density of the fluid
/// @param lenght  The lenght of the box
/// @param atoms List of atoms
/// @param name The atom name ("Ar" or "Si")
/// @return The quantity of atoms placed
int generateAtoms(double density, double lenght, atom **atoms, int element)
{
    printf("Generating structure...\n");

    // Find the number of atoms in the box
    int n = (int)ceil(density * pow(lenght, 3));

    // Find the cubic grid M
    int M = (int)ceil(pow(n, 1.0 / 3.0));
    // Populate the cells with atoms
    // Cuantity of atoms to be placed
    int mCube = (int)ceil(pow(M, 3));
    // Allocate n space for the temporal atoms list
    atom *atemp = (atom *)calloc(mCube, sizeof(atom));

    int aCount = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int k = 0; k < M; k++)
            {
                vector u;
                u.x = lenght / M * (i + 1);
                u.y = lenght / M * (j + 1);
                u.z = lenght / M * (k + 1);
                atemp[aCount] = atomInit(element, u, aCount + 1);
                aCount++;
            }
        }
    }

    // Now random delete unneeded atoms
    // Allocate n space for the final atoms list
    *atoms = (atom *)calloc(n, sizeof(atom));
    // Atoms to be deleted
    int del = mCube - n;

    // Positions to be deleted
    int *randomDel = (int *)calloc(del, sizeof(int));
    int d = 0;

    while (d < del)
    {
        int rndInt = getRandomInt(0, mCube);
        // Check if there are repeated indices
        int repeated = 0;
        for (int i = 0; i < d; i++)
        {
            if (randomDel[i] == rndInt)
            {
                repeated = 1;
            }
        }
        if (repeated == 0)
        {
            randomDel[d] = rndInt;
            d++;
        }
    }

    // Restart the counter
    aCount = 0;
    // Copy atoms to the new final array
    for (int i = 0; i < mCube; i++)
    {
        int continueLoop = 0;
        for (int j = 0; j < d; j++)
        {
            if (randomDel[j] == i)
            {
                // Skip one atom (therefore deleting it from the new array)
                continueLoop = 1;
                break;
            }
        }
        if (continueLoop == 1)
        {
            continue;
        }
        copyAtom(atemp[i], &(*atoms)[aCount]);
        //(*atoms)[aCount] = atemp[i];
        aCount++;
    }
    free(atemp);
    printf("Total atoms placed: %d\n", n);
    return n;
}

void printAtom(atom a, int n)
{
    printf("%d %s %lf %lf %lf\n", n, a.name, a.position.x, a.position.y, a.position.z);
}

void printAtomList(atom *list, int size)
{
    for (int i = 0; i < size; i++)
    {
        printAtom(list[i], i);
    }
}

void trialIteration(mcconfig *config, double *energy, int *changes)
{
    int r = getRandomInt(0, config->natoms - 1);
    // copyAtomList(config->atoms, &newPositions, config->natoms);
    double prevEnergy = singleLennardJones(config->atoms, config->natoms, r, 0, config->l, config->squared_cutoff);
    vector oldpos = moveRandomAtom(&config->atoms, r, config->l, config->step);
    double postEnergy = singleLennardJones(config->atoms, config->natoms, r, 0, config->l, config->squared_cutoff);
    double change = postEnergy - prevEnergy;

    int printAtoms = 0;

    if (change <= 0) // Accept change
    {
        printAtoms = 1;
    }
    else // Distribution
    {
        double p = monteCarloProbability(change, config->temp);
        double w = getRandomDouble(0, 1);
        if (p >= w)
        {
            printAtoms = 1;
        }
    }
    if (printAtoms)
    {
        (*energy) += change;
        (*changes)++;
    }
    else
    {
        config->atoms[r].position = oldpos;
    }
}
// changes count to update the neighbours accordingly
// squared_max_dist (rskin - cutoff) ^ 2
void trialNeiIteration(mcconfig *config, double *energy, int *changes, double squared_max_dist)
{
    int r = getRandomInt(0, config->natoms - 1);
    double prevEnergy = singleNeiLennardJones(config->atoms, config->natoms, r, config->l, config->squared_cutoff);
    vector oldpos = moveRandomAtom(&config->atoms, r, config->l, config->step);
    double postEnergy = singleNeiLennardJones(config->atoms, config->natoms, r, config->l, config->squared_cutoff);
    double change = postEnergy - prevEnergy;

    int printAtoms = 0;

    if (change <= 0) // Accept change
    {
        printAtoms = 1;
    }
    else // Distribution
    {
        double p = monteCarloProbability(change, config->temp);
        double w = getRandomDouble(0, 1);
        if (p >= w)
        {
            printAtoms = 1;
        }
    }
    if (printAtoms)
    {
        // free(config->atoms);
        // copyAtomList(newPositions, &config->atoms, config->natoms);
        (*energy) += change;
        (*changes)++;

        atom movedAtom = config->atoms[r];
        vector v_moved = substractVector(movedAtom.position, movedAtom.nei_position);
        v_moved = MIVector(v_moved, config->l);
        double squared_moved = sumSquaredComponents(v_moved);
        if (squared_moved >= squared_max_dist)
        {
            // printf("Updating neighbours after %d changes\n", *changes);
            updateNeighbours(&config->atoms, config->natoms, config->l, config->squared_rskin);
        }
    }
    else
    {
        config->atoms[r].position = oldpos;
    }
    // free(newPositions);
}

grdist gr_init(mcconfig config)
{
    grdist gr;
    gr.naccum = 0;
    gr.natoms = config.natoms;
    gr.dens = config.density;

    gr.resolution = 400;
    gr.rmax = 4;
    gr.rmax_squared = gr.rmax * gr.rmax;
    gr.deltaR = gr.rmax / gr.resolution;

    gr.n = (int *)calloc(gr.resolution, sizeof(int));
    for (int i = 0; i < gr.resolution; i++)
    {
        gr.n[i] = 0;
    }

    return gr;
}

void save_gr(char *name, grdist gr)
{
    FILE *fptr;
    fptr = fopen(name, "w");
    if (fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }

    double volPI = atan(1.0) * 4;
    for (int i = 0; i < gr.resolution; i++)
    {
        int ni = gr.n[i];
        double ri = i * gr.deltaR;
        double rterm = pow(ri + gr.deltaR, 3) - pow(ri, 3);
        double result = (ni * 3) / (gr.naccum * (gr.natoms / 2) * 4 * volPI * rterm * gr.dens);
        fprintf(fptr, "%lf %lf\n", ri, result);
    }

    fclose(fptr);
}

void update_gr(grdist *gr, mcconfig config)
{
    gr->naccum++;

    for (int j = 0; j < config.natoms; j++)
    {
        for (int k = j + 1; k < config.natoms; k++)
        {
            double dist = shortestAtomDistance(config.atoms[j], config.atoms[k], config.l, gr->rmax_squared);
            if (dist < 0)
            {
                continue; // Skip because distance is over the max radius of g(r)
            }
            for (int l = 0; l < gr->resolution; l++)
            {
                double currentR = l * gr->deltaR;
                double nextR = currentR + gr->deltaR;
                if (dist < nextR && dist >= currentR)
                {
                    gr->n[l]++;
                }
            }
        }
    }
}

void performMonteCarlo(mcconfig config)
{
    printf("%s\n", config.name);

    // Check if the configuration has initial structure in order to generate a new one
    if (config.hasInitial == 0)
    {
        int element = 18;
        if (config.type == 2)
        {
            element = 14;
        }
        config.natoms = generateAtoms(config.density, config.l, &config.atoms, element);
    }

    // Prepare the output file
    FILE *fptr;
    char *outputxyz = (char *)calloc(20, sizeof(char));
    strcpy(outputxyz, config.name);
    strcat(outputxyz, ".xyz");
    fptr = fopen(outputxyz, "w");
    if (fptr == NULL)
    {
        printf("Error!");
        exit(1);
    }

    // Initialize the distribution function class
    grdist gr = gr_init(config);
    // Calculate Lennard-Jones energy
    double energy = fullLennardJones(config.atoms, config.natoms, config.l, config.squared_cutoff);
    printf("Initial energy: %lf\n", energy);

    int iterations = (config.equilib + config.production) * config.natoms;
    int changes = 0; // Count accepted moves

    // Print in terminal the total number of iterations
    printf("Total iterations: %d\nEquilibrium: %d\nProduction: %d\n", iterations, config.equilib, config.production);

    // Save initial positions
    printXYZFile(config.atoms, config.natoms, config.name, fptr);

    // The actual monte-carlo iterations
    if (config.useNei == 0)
    {
        for (int i = 0; i < iterations; i++)
        {
            trialIteration(&config, &energy, &changes);
            if (i > config.equilib * config.natoms)
            {
                if (i % config.natoms == 0) // Save each sweep
                {
                    //printXYZFile(config.atoms, config.natoms, config.name, fptr);
                    update_gr(&gr, config);
                }
            }
        }
    }
    else if (config.useNei == 1)
    {
        printf("Using neighbours\n");

        double max_squared_dist = pow(sqrt(config.squared_rskin) - sqrt(config.squared_cutoff), 2);
        updateNeighbours(&config.atoms, config.natoms, config.l, config.squared_rskin);
        for (int i = 0; i < iterations; i++)
        {
            trialNeiIteration(&config, &energy, &changes, max_squared_dist);
            if (i > config.equilib * config.natoms)
            {
                if (i % config.natoms == 0) // Save each sweep
                {
                    //printXYZFile(config.atoms, config.natoms, config.name, fptr);
                    update_gr(&gr, config);
                }
            }
        }
    }

    printf("Accepted %.2lf%% of the changes\n", ((1.0 * changes) / iterations) * 100);
    printf("Final energy with changes: %lf\n", energy);
    printf("Saved arrangement has full %lf energy\n", fullLennardJones(config.atoms, config.natoms, config.l, config.squared_cutoff));

    // Close the output file
    fclose(fptr);

    // Save the distribution function
    char *outputgidst = config.name;
    strcat(outputgidst, ".gdist");
    save_gr(outputgidst, gr);
}