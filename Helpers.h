#include <iostream>            // for cout and cin
#include <stdlib.h>     /* abs */
#include <vector>
#include <array>
#include <math.h>       /* sqrt */
#include <string>
 
double sum_coo(std::array<double, 3> a);

std::array<double, 3> colMax(std::vector<std::array<double, 3>> a);

std::array<double, 3> colMin(std::vector<std::array<double, 3>> a);

double dot(std::array<double, 3> a, std::array<double, 3> b);

std::array<double, 3> cross(std::array<double, 3> a,std::array<double, 3> b);

std::array<double, 3> minus(std::array<double, 3> a, std::array<double, 3> b);

void normalize(std::array<double, 3>& v);

double randZO();

std::array<double, 3> getRandOrth(std::array<double, 3> a);

double isLeft( std::array<double, 3> P0,std::array<double, 3> P1, std::array<double, 3> P2 );

std::array<double, 3> faceNormal( std::array<double, 3> V0,std::array<double, 3> V1, std::array<double, 3> V2 );

void printArray(std::array<double, 3> v);




