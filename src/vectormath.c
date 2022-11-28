#include "vectormath.h"
#include <stdlib.h>

vector addVector(vector u, vector v)
{
    vector result;

    result.x = u.x + v.x;
    result.y = u.y + v.y;
    result.z = u.z + v.z;

    return result;
}

vector substractVector(vector u, vector v)
{
    vector result;

    result.x = u.x - v.x;
    result.y = u.y - v.y;
    result.z = u.z - v.z;

    return result;
}

vector scaledVector(double a, vector u)
{
    vector result;

    result.x = a * u.x;
    result.y = a * u.y;
    result.z = a * u.z;

    return result;
}

vector addToComponents(double a, vector u)
{
    vector result;

    result.x = a + u.x;
    result.y = a + u.y;
    result.z = a + u.z;

    return result;
}

double sumSquaredComponents(vector u)
{
    return u.x * u.x + u.y * u.y + u.z * u.z;
}

double norm(vector u)
{
    return sqrt(sumSquaredComponents(u));
}

void normalize(vector *u)
{
    double n = norm(*u);
    u->x /= n;
    u->y /= n;
    u->z /= n;
}

double getRandomDouble(double lower, double upper)
{
    double range = (upper - lower);
    double div = RAND_MAX / range;
    return lower + (rand() / div);
}

int getRandomInt(int lower, int upper)
{
    return (rand() % (upper - lower + 1)) + lower;
}