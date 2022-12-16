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
    mcconfig config = initMCConfig();
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
    return n;
}

void trialIteration(mcconfig *config, double *energy, int *changes, double squared_max_dist, double (*potential)())
{
    // Get a random integer between 0 and the total number of atoms - 1
    int r = getRandomInt(0, config->natoms - 1);

    // Compute the energy corresponding to atom r
    double prevEnergy = potential(config->atoms, r, config->l, config->squared_cutoff);

    // Change the position of atom r and saving the old position
    vector oldpos = moveRandomAtom(&config->atoms, r, config->l, config->stepSize);

    // Compute the new energy after moving atom r
    double postEnergy = potential(config->atoms, r, config->l, config->squared_cutoff);

    // The change of energy
    double change = postEnergy - prevEnergy;

    // To control whether to save or not the new arrangement
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

    // Save the new position
    if (printAtoms)
    {
        (*energy) += change;
        (*changes)++;
        double moved = shortestMIDistance(config->atoms[r].nei_position, config->atoms[r].position, config->l, squared_max_dist, NULL);
        if (moved == -1) // Update the neighbours if the distance moved by this atoms is superior to rskin / 2 (squared_max_dist)
        {
            updateNeighbours(&config->atoms, config->natoms, config->l, config->squared_rskin_plus_cutoff);
        }
    }
    else // Erase the new position
    {
        config->atoms[r].position = oldpos;
    }
}

grdist gr_init(mcconfig config)
{
    grdist gr;
    gr.naccum = 0;
    gr.natoms = config.natoms;
    gr.dens = config.density;

    gr.resolution = config.gr_resolution;
    gr.rmax_squared = config.gr_rmax_squared;
    gr.deltaR = sqrt(gr.rmax_squared) / gr.resolution;

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

    double volPI = atan(1.0) * 4.0;
    for (int i = 0; i < gr.resolution; i++)
    {
        int ni = gr.n[i];
        double ri = i * gr.deltaR;
        double rterm = pow(ri + gr.deltaR, 3) - pow(ri, 3);
        double result = (ni * 3.0) / (gr.naccum * (gr.natoms / 2.0) * 4.0 * volPI * rterm * gr.dens);
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
            double dist = shortestAtomDistance(config.atoms[j], config.atoms[k], config.l, gr->rmax_squared, NULL);
            if (dist > 0)
            {
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
}

void performMonteCarlo(mcconfig config)
{
    printf("============Monte Carlo===========\n");
    printf("Starting '%s' simulation\n", config.name);

    // Check if the configuration has initial structure in order to generate a new one
    if (config.hasInitial == 0)
    {
        int element = 18;
        if (config.type == 2)
        {
            element = 14;
        }
        printf("Generating structure... ");
        config.natoms = generateAtoms(config.density, config.l, &config.atoms, element);
        printf("Total atoms placed: %d\n", config.natoms);
    }
    else
    {
        printf("Read structure from '%s'\n", config.input);
    }

    // Max atom movement before updating neighbours again rskin / 2, but squared.
    double max_squared_dist = config.squared_rskin / 4.0;
    printf("Cutoff: %f\n", sqrt(config.squared_cutoff));
    printf("rskin: %f\n", sqrt(config.squared_rskin));
    printf("Neighbours up to: %f\n", sqrt(config.squared_rskin_plus_cutoff));
    printf("Update after: %f\n", sqrt(max_squared_dist));

    // Prepare the output file
    FILE *fptr;
    char *outputxyz = (char *)calloc(20, sizeof(char));
    strcpy(outputxyz, config.name);
    strcat(outputxyz, ".xyz");
    fptr = fopen(outputxyz, "w");
    if (fptr == NULL)
    {
        printf("Error generating output xyz...\n");
        exit(1);
    }

    // Initialize the distribution function class
    grdist gr = gr_init(config);

    int equilib_it = config.equilib * config.natoms;
    int prod_it = config.production * config.natoms;
    int total_it = prod_it + equilib_it;
    int changes = 0; // Count accepted moves

    // Print in terminal the total number of iterations
    printf("============Iterations============\n");
    printf("Equilibrium: %d\nProduction: %d\nTotal: %d\n", equilib_it, prod_it, total_it);

    updateNeighbours(&config.atoms, config.natoms, config.l, config.squared_rskin_plus_cutoff);

    double energy = 0;    // Initialize energy variable
    double completed = 0; // Keep track of the simulation status

    printf("===============Model==============\n");
    // Assign function pointer so code can be reused
    double (*potential)();
    double (*fullPotential)();
    if (config.type == 1)
    {
        printf("Lennard-Jones model\n");
        potential = &singleNeiLennardJones;
        fullPotential = &fullLennardJones;
    }
    else if (config.type == 2)
    {
        printf("Stillinger-Weber model\n");
        potential = &stillingerModel;
        fullPotential = &fullStillinger;
    }
    else
    {
        printf("Error reading simulation type...\n");
        exit(1);
    }

    // Calculate full energy
    energy = fullPotential(config.atoms, config.natoms, config.l, config.squared_cutoff);
    printf("Initial energy: %lf\n", energy);

    printf("============Simulating============\n");

    // Perform the Monte Carlo equlibration iterations
    for (int i = 0; i < equilib_it; i++)
    {
        trialIteration(&config, &energy, &changes, max_squared_dist, potential);
        if (i % config.natoms == 0) // Save each sweep
        {
            completed = ((double)i / equilib_it) * 100.0;
            printf("\rEquilibration: %.1f%% ", completed);
            fflush(stdout);
        }
    }
    printf("\n");
    // Perform the Monte Carlo production iterations
    for (int i = 0; i < prod_it; i++)
    {
        trialIteration(&config, &energy, &changes, max_squared_dist, potential);
        if (i % config.natoms == 0) // Save each sweep
        {
            completed = ((double)i / prod_it) * 100.0;
            printf("\rProduction: %.1f%% ", completed);
            fflush(stdout);
            update_gr(&gr, config);
            // Uncoment to save a xyz file with a snapshot of each sweep
            //printXYZFile(config.atoms, config.natoms, config.name, config.sigma, fptr);
        }
    }

    printf("\n=============Results==============\n");
    printf("Accepted %.2lf%% of the changes\n", ((1.0 * changes) / total_it) * 100);
    printf("Final energy: %lf\n", energy);

    updateNeighbours(&config.atoms, config.natoms, config.l, config.squared_rskin_plus_cutoff);

    // Close the output file
    fclose(fptr);

    // Save the distribution function
    char *outputgidst = config.name;
    strcat(outputgidst, ".gdist");
    save_gr(outputgidst, gr);
    printf("Saved %s\n", outputgidst);
    printf("=========Normal Termination=======\n");
}

mcconfig initMCConfig()
{
    mcconfig config;
    config.name = "Montecarlo";
    config.hasInitial = 0;
    config.natoms = 0;
    // config.atoms = (atom *)calloc(1, sizeof(atom));
    config.type = 0;
    config.l = 1;
    config.density = 1;
    config.temp = 1;
    config.sigma = 1;
    config.squared_cutoff = 1;
    config.squared_rskin_plus_cutoff = sqrt(2);
    config.stepSize = 1;
    config.useNei = 0;
    config.squared_rskin = 1;
    config.equilib = 0;
    config.production = 0;

    return config;
}