#include <iostream>            // for cout and cin
#include <stdlib.h>     /* abs */
#include <vector>
#include <array>
#include <math.h>       /* sqrt */
#include <string>
#include <sstream>
#include <fstream>
#include <Helpers.h>

class Mesh                      // begin declaration of the class
{
  public:                      // begin public section

    void calcProjectedPts(std::array<double, 3> proj);        // accessor function
    void computeFaceNorm();
    std::array<double, 3> computeSurfRad(double reso);
    void getContainedIdx(int face_idx, std::vector<int>& idx, std::vector<double>& z_val);
    void loadAsciiPly(std::string filename);
    void loadX3D(std::string fname);
    void dumpZVal(std::string fname);
    void dumpF(std::string fname);
    int n_pts;
    int n_faces = -2;
    std::vector<std::array<double, 3>> pts;
    std::vector<std::array<double, 3>> proj_pts;
    std::vector<std::array<int, 3>> faces; // only triangular faces 
    std::vector<std::array<double, 3>> faces_norm;
    std::vector<std::array<double, 2>> props; // rho and delta for all faces
    std::vector<std::array<double, 4>> colors; // RGB alpha
 private:                      // begin private section
    int n_x;
    int n_y;
    double res;
    std::array<double, 3> proj_x;
    std::array<double, 3> proj_y;
    std::array<double, 3> proj_z;
    std::vector<double> z_pix;
    std::vector<double> fx;
    std::vector<double> fy;
    std::vector<double> fz;
};
