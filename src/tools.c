#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tools.h"

/// @brief Checks if the given string starts with the given prefix
/// @param a full string
/// @param b prefix
/// @return 0 if false, 1 if true
int startsWith(const char *a, const char *b)
{
    if (strncmp(a, b, strlen(b)) == 0)
        return 1;
    return 0;
}

atom readAtomFromXYZ(char *line, double sigma, int natoms)
{
    atom a;
    vector u;
    char delim[] = " ";

    char *ptr = strtok(line, delim);
    int n = 0;

    // Needed for the strtod function
    char *extra;

    while (ptr != NULL)
    {
        // Get the element number of the atom
        if (n == 0)
        {
            a.name = ptr;
            a.element = getElementFromName(ptr);
        }
        // x, y ,z values
        if (n == 1)
        {
            u.x = strtod(ptr, &extra) / sigma;
        }
        if (n == 2)
        {
            u.y = strtod(ptr, &extra) / sigma;
        }
        if (n == 3)
        {
            u.z = strtod(ptr, &extra) / sigma;
        }

        ptr = strtok(NULL, delim);
        n++;
    }
    a.position = u;
    a = atomInit(a.element, a.position, natoms);
    return a;
}

void readXYZ(mcconfig *config)
{
    FILE *file;
    file = fopen(config->input, "r");
    if (NULL == file)
    {
        printf("File '%s' can't be opened \n", config->input);
        exit(0);
    }
    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    int ln = 0;

    atom *atoms;
    int n = 0;

    while ((read = getline(&line, &len, file)) != -1)
    {
        // Increment line number
        ln++;

        // Reading first line, which contains number of atoms
        // Allocates the memory for the number of atoms
        if (ln == 1)
        {
            atoms = (atom *)calloc(atoi(line), sizeof(atom));
        }
        // Title is on the second line
        else if (ln == 2)
        {
            continue;
        }
        else
        {
            atoms[n] = readAtomFromXYZ(line, config->sigma, n);
            n++;
        }
    }

    config->natoms = n;
    config->atoms = atoms;

    centerAtomsToBox(config->atoms, config->natoms);
    config->l = boxSize(config->atoms, config->natoms);
}

void readInput(FILE *file, mcconfig *config)
{
    char *line = NULL;
    size_t len = 0;
    ssize_t read;

    char *extra; // Nedded for the strod function

    // Initialize hasInitial
    config->hasInitial = 0;

    int hasBoxSize = 0;
    double boxsize = 0;

    while ((read = getline(&line, &len, file)) != -1)
    {

        line[strcspn(line, "\n")] = 0; // Removing trailing "\n" character
        if (startsWith(line, "\%name="))
        {
            config->name = (char *)calloc(len, sizeof(char)); // Allocating memory for the name
            strcpy(config->name, line + 6);
        }
        if (startsWith(line, "\%type="))
        {
            char *name = line + 6;
            if (strcmp(name, "lennard") == 0)
            {
                config->type = 1;
            }
            else if (strcmp(name, "stillinger") == 0)
            {
                config->type = 2;
            }
            else
            {
                printf("Unsupported type. Allowed types are 'lennard' and 'stillinger'\n");
            }
        }
        if (startsWith(line, "\%neighbours="))
        {
            char *nei = line + 12;
            if (strcmp(nei, "true") == 0)
            {
                config->useNei = 1;
            }
            else
            {
                config->useNei = 0;
            }
        }
        if (startsWith(line, "\%sigma="))
        {
            config->sigma = strtod(line + 7, &extra);
        }
        if (startsWith(line, "\%density="))
        {
            config->density = strtod(line + 9, &extra);
        }
        if (startsWith(line, "\%cutoff="))
        {
            config->squared_cutoff = strtod(line + 8, &extra);
        }
        if (startsWith(line, "\%equilibration="))
        {
            config->equilib = atoi(line + 15);
        }
        if (startsWith(line, "\%production="))
        {
            config->production = atoi(line + 12);
        }
        if (startsWith(line, "\%stepsize="))
        {
            config->step = strtod(line + 10, &extra);
        }
        if (startsWith(line, "\%temp="))
        {
            config->temp = strtod(line + 6, &extra);
        }
        if (startsWith(line, "\%input="))
        {
            config->hasInitial = 1;
            config->input = (char *)calloc(len, sizeof(char)); // Allocating memory for the input file name
            strcpy(config->input, line + 7);
        }
        if (startsWith(line, "\%boxsize="))
        {
            hasBoxSize = 1;
            boxsize = strtod(line + 9, &extra);
            config->l = boxsize;
        }
        if (startsWith(line, "\%grmax="))
        {
            config->gr_rmax_squared = strtod(line + 7, &extra);
            config->gr_rmax_squared *= config->gr_rmax_squared;
        }
        if (startsWith(line, "\%res="))
        {
            config->gr_resolution = atoi(line + 5);
        }
    }

    // Check initial structure
    if (config->hasInitial == 1)
    {
        readXYZ(config);
    }

    if (hasBoxSize)
    {
        config->l = boxsize;
    }

    // Finalize cutoff setup
    config->squared_rskin = config->step * 6; // 6 times was the fastest as the tests shown
    config->squared_rskin_plus_cutoff = config->squared_cutoff + config->squared_rskin;

    // Square all the variables
    config->squared_rskin *= config->squared_rskin;
    config->squared_rskin_plus_cutoff *= config->squared_rskin_plus_cutoff;
    config->squared_cutoff *= config->squared_cutoff;
}

const char *get_filename_ext(const char *filename)
{
    const char *dot = strrchr(filename, '.');
    if (!dot || dot == filename)
        return "";
    return dot + 1;
}

int fileReader(char *url, mcconfig *config)
{
    FILE *file;
    const char *extension;

    file = fopen(url, "r");
    if (NULL == file)
    {
        printf("file can't be opened \n");
        return 1;
    }

    extension = get_filename_ext(url);

    // if (strcmp(extension, "xyz") == 0) {
    //     printf("Reading xyz file:\n");
    //     readXYZ(file, config);
    // }
    if (strcmp(extension, "mc") == 0)
    {
        printf("Reading montecarlo file...\n\n");
        readInput(file, config);
    }
    else
    {
        printf("Incompatible file\n");
        return 1;
    }

    fclose(file);
    return 0;
}

void printXYZFile(atom *atoms, int natoms, char *name, double sigma, FILE *file)
{
    fprintf(file, "%d\n", natoms);
    fprintf(file, "%s\n", name);
    for (int j = 0; j < natoms; j++)
    {
        fprintf(file, "%s %f %f %f\n", atoms[j].name, atoms[j].position.x * sigma, atoms[j].position.y * sigma, atoms[j].position.z * sigma);
    }
}