#include <stdio.h>
#include "montecarlo.h"
#include <time.h>

int main(int argc, char **argv)
{
    // Setting seed for random rumber generation
    srand(time(NULL));
    if (argc > 1)
    {
        mcconfig config = readMCFromFile(argv[1]);
        performMonteCarlo(config);
    }
    else 
    {
        printf("Usage: mc.exe <input.mc>\n");
    }

    printf("Program finished.\n");

    return 0;
}