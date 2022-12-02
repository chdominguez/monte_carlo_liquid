#include<math.h>
/*
 * Struct:  Vector 
 * --------------------
 * Store data about the x, y and z positions in space
 *
 *  x, y ,z: doubles containing the position
 */
typedef struct Vector {
    double x;
    double y;
    double z;
} vector;

/*
 * Function:  addVector 
 * --------------------
 * Adds a pair of vectors
 *
 *  u: the first vector
 *  v: the second vector
 * 
 *  returns the vector resulting from u + v
 */
vector addVector(vector u, vector v);

/*
 * Function:  substractVector 
 * --------------------
 * Substracts a pair of vectors
 *
 *  u: the first vector
 *  v: the second vector
 * 
 *  returns the vector resulting from u - v
 */
vector substractVector(vector u, vector v);

/*
 * Function:  scaledVector 
 * --------------------
 * Scales a vector
 *
 *  a: the scalar value
 *  v: the vector
 * 
 *  returns the vector resulting from a*v
 */
vector scaledVector(double a, vector u);

/*
 * Function:  addToComponents 
 * --------------------
 * Adds the scalar to each of the components of the vector
 *
 *  a: the scalar value
 *  v: the vector
 * 
 *  returns the vector resulting from a+v.x, a+v.y, a+v.z
 */
vector addToComponents(double a, vector u);

/// @brief Returns the squared sum of the components
/// @param u The vector
/// @return u.x + u.y + u.z
double sumSquaredComponents(vector u);

/*
 * Function:  norm 
 * --------------------
 * Computes the norm of the vector
 *  
 *  u: the vector
 *  
 *  returns sqrt(x^2+y^2+z^2)
 */
double norm(vector u);

/*
 * Function:  normalize 
 * --------------------
 * Modifies a vector in order to normalize it
 *
 *  u: the vector to normalize
 */
void normalize(vector *u);

/*
 * Function:  getRandomDouble
 * --------------------
 * Returns a random double number between a range
 *
 *  lower: the lower number in the range
 *  upper: the upper number in the range
 */
double getRandomDouble(double lower, double upper);

/*
 * Function:  getRandomInt 
 * --------------------
 * Returns a random int number between a range
 *
 *  lower: the lower number in the range
 *  upper: the upper number in the range
 */
int getRandomInt(int lower, int upper);

/// @brief Compute the angle between two vectors
/// @param a First vector
/// @param b Second vector
/// @return The angle in radians
double angleBetweenVectors(vector a, vector b);