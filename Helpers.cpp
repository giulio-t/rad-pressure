#include <iostream>            // for cout and cin
#include <stdlib.h>     /* abs */
#include <vector>
#include <array>
#include <math.h>       /* sqrt */
#include <cmath>        // std::abs
#include <string>
#include <Helpers.h>
 
 double sum_coo(std::array<double, 3> a)
 {
 return a[0] + a[1] + a[2];
 }

std::array<double, 3> colMax(std::vector<std::array<double, 3>> a)
{ 
std::array<double, 3> m = {0, 0, 0};
int l = a.size();
for (int i=0; i < l; i++){
	m[0] =std::max(m[0],a[i][0]);
	m[1] =std::max(m[1],a[i][1]);
	m[2] =std::max(m[2],a[i][2]);
	}
return m;
}

std::array<double, 3> colMin(std::vector<std::array<double, 3>> a)
{ 
std::array<double, 3> m = {0, 0, 0,};
int l = a.size();
for (int i=0; i < l; i++){
	m[0] =std::min(m[0],a[i][0]);
	m[1] =std::min(m[1],a[i][1]);
	m[2] =std::min(m[2],a[i][2]);
	}
return m;
}

double dot(std::array<double, 3> a, std::array<double, 3> b)
{ 
double d = 0;
for (int i=0; i < 3; i++){
	d = d + a[i]*b[i];
}
return d;

}

std::array<double, 3> cross(std::array<double, 3> a,std::array<double, 3> b)
{
std::array<double, 3> c;
c[0] = a[1]*b[2] - a[2]*b[1];
c[1] = a[2]*b[0] - a[0]*b[2];
c[2] = a[0]*b[1] - a[1]*b[0];
return c;
}

std::array<double, 3> minus(std::array<double, 3> a, std::array<double, 3> b)
{
std::array<double, 3> c;
c[0] = a[0] - b[0];
c[1] = a[1] - b[1];
c[2] = a[2] - b[2];
return c;
}

void normalize(std::array<double, 3>& v)
{
double norm = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
//std::cout << "norm: " << norm<< "\n";
v[0] = v[0] / norm;
v[1] = v[1] / norm;
v[2] = v[2] / norm;
}

double randZO()
{
    return rand() / (RAND_MAX + 1.);
}

std::array<double, 3> getRandOrth(std::array<double, 3> a)
{
std::array<double, 3> o;
double proj = 1;
while (std::abs(proj) > 0.99) // dunno if needed anyway it does not hurt
{ 
	o[0] = (randZO()-0.5)*2;
	o[1] = (randZO()-0.5)*2;
	o[2] = (randZO()-0.5)*2;
	normalize(o);
	proj = dot(o,a);
}
o[0] = o[0] - proj*a[0];
o[1] = o[1] - proj*a[1];
o[2] = o[2] - proj*a[2];
normalize(o);
return o;
}



// isLeft(): test if a point is Left|On|Right of an infinite 2D line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 to P1
//          =0 for P2 on the line
//          <0 for P2 right of the line
double isLeft( std::array<double, 3> P0,std::array<double, 3> P1, std::array<double, 3> P2 )
{
    return ( (P1[0] - P0[0]) * (P2[1] - P0[1])
           - (P2[0] - P0[0]) * (P1[1] - P0[1]) );
}
std::array<double, 3> faceNormal( std::array<double, 3> V0,std::array<double, 3> V1, std::array<double, 3> V2 )
{
std::array<double, 3> n;
n = cross(V0,V1);
}

void printArray(std::array<double, 3> v)
{
int l = v.size();
for (int i =0; i < l ; i++){
	std::cout << v[i] <<" ";
}
std::cout << "\n";
}





