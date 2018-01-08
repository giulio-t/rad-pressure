#include <iostream>
#include <string>
#include <Mesh.h>

int main(int argc, char** argv)
{	
    double res;
    if (argc > 1)
    {
    	res = atof(argv[2]);
    }else{
    	res = 0.01;
    }
    /**std::array<double,3> a = { 0.4,,0};
    std::array<double,3> b = { 0,1,0};
    std::cout << "dot:" << dot(a,b);
    std::cout << "cross:";
    printArray(cross(a,b));**/
    std::string fname(argv[1]);
    std::array<double,3> prj ={ 0,0,-1};
    normalize(prj);
    Mesh mesh;
    mesh.loadX3D(fname);
    mesh.calcProjectedPts(prj);
    mesh.computeFaceNorm();
    std::array<double, 3> F = mesh.computeSurfRad(res);
    printArray(F);
    mesh.dumpZVal("zvals.txt");
    return 0;
    
}
