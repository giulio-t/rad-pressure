#include <iostream>            // for cout and cin
#include <stdlib.h>     /* abs */
#include <vector>
#include <array>
#include <math.h>       /* sqrt */
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <Mesh.h>

#define EPS 1e-9
#define EPS2 1e-4
void Mesh::loadAsciiPly(std::string fname)
{
std::ifstream input(fname);
bool header = true;
int v_count = 0;
int p_count = 0;
int t_count = 0;
for( std::string line; getline( input, line ); )
{	std::stringstream stream(line);
	if (header)
	{
		if (line.find("end_header"))
		{
			header = false;
		}
		else if (line.find("element vertex"))
		{
			n_pts = std::stoi(line.substr(14));
			pts.resize(n_pts);
		}
	} else if(v_count < n_pts){
		//sscanf (line.c_str(),"%f %f %f",&pts[v_count][0],&pts[v_count][1],&pts[v_count][2]);
		stream >> pts[v_count][0];
		stream >> pts[v_count][1];
		stream >> pts[v_count][2];
	} else {
		std::vector<int> vertex;
		int n_v;
		stream >> n_v;
		int n; 
		for(int i=0; i < n_v ; i++) {
		stream >> n;
		vertex.push_back(n);}
		// convert polygon in traingle
		for (int i =0; i < (n_v - 2); i++){
			std::array<int, 3> tri;
			tri[0] = vertex[0];
			tri[1] = vertex[i+1];
			tri[2] = vertex[i+2];
			faces.push_back(tri);
			t_count++;
		}
	}
}
n_faces = t_count;
}
void Mesh::loadX3D(std::string fname) // warnig does only supports 1 shape with color per faces 
{
	std::ifstream input(fname);
	bool main = true;
	bool read = false;
	int v_count = 0;
	int p_count = 0;
	int t_count = 0;
	std::string obj;
	std::string coo;
	for( std::string line; getline( input, line ); )
	{	std::stringstream stream(line);
		if (main)
		{
			if (line.find("<IndexedFaceSet") != std::string::npos )
			{
				main = false;
				read = true;
			}
		
		} 
		if (read){
			obj.append(line);
			if  (line.find("</IndexedFaceSet>") != std::string::npos ){
			break;
			}
		
		}
	
	}
	// probably extremly inefficient
	int idxs = obj.find("colorPerVertex=\"")+16;
	int idxe = obj.find("\"",idxs);
	std::string cpv  = obj.substr(idxs,idxe - idxs);
	
	if (not(cpv.find("false" )!= std::string::npos)){
		std::cout <<  "ERROR!! color on vertex instead of faces\n";
		return;
	}
	
	idxs = obj.find("coordIndex=\"") + 12;
	idxe = obj.find("\"",idxs);

	std::string face  = obj.substr(idxs,idxe - idxs);
	
	idxs = obj.find("<Coordinate point=\"")+19;
	idxe = obj.find("\"",idxs);
	std::string vertex  = obj.substr(idxs,idxe - idxs);
	
	idxs = obj.find("<ColorRGBA color=\"")+18;
	idxe = obj.find("\"",idxs);
	std:: string color  = obj.substr(idxs,idxe - idxs);
	
	
	
	// reading vertex into pts
	std::stringstream stream;
	stream.str(vertex);
	double c;
	std::array<double, 3> p;
	int n_c = 0;
	int n_p = 0;
	while (stream >> c){
		p[n_c] = c;
		if (n_c == 2){
			//std::cout << p[0] << " " << p[1] << " " << p[2]<<"\n";
			pts.push_back(p);
			n_p++;
			n_c = -1;
		}
		n_c++;
		
	}
	n_pts = n_p;
	// reading idx into faces
	stream.str("");
 	stream.clear();
	stream.str(face);
	int idxp;
	std::array<int, 3> tf;
	std::array<int, 3> f;
	n_p = 0;
	int n_f = 0;
	while (stream >> idxp){
		tf[n_p] = idxp;
		if (n_p  == 2){
			//std::cout << f[0] << " " << f[1] << " " << f[2]<<"\n";
			
			stream >> idxp;
			if (idxp >0){
				f[0] = tf[0];
				f[1] = tf[2];
				f[2] = tf[1];
			}else{
				f = tf;
			
			}
			faces.push_back(f);
			n_f++;
			n_p = -1;
		}
		n_p++;
		
	}
	n_faces = n_f;
	// reading color into colors
	stream.str("");
 	stream.clear();
	stream.str(color);
	n_c = 0;
	n_f = 0;
	std::array<double, 4> cf;
	while (stream >> c){
		cf[n_c ] = c;
		
		if (n_c  == 3){
			colors.push_back(cf);
			n_f++;
			n_c = -1;
		}
		n_c++;
		
	}
	
	// reading normals into normlas
	idxs = obj.find("<Normal vector=\"")+16;
	if (idxs != std::string::npos ){
	idxe = obj.find("\"",idxs);
	std:: string normals  = obj.substr(idxs,idxe - idxs);
	stream.str("");
 	stream.clear();
	stream.str(normals);
	n_c = 0;
	n_f = 0;
	while (stream >> c){
		p[n_c ] = c;
		
		if (n_c  == 2){
			faces_norm.push_back(p);
			n_f++;
			n_c = -1;
		}
		n_c++;
		
	}
	}
	std::cout << "Mesh loaded, number of vertex: " << n_pts << " number of triangles: " << n_faces <<"\n";
	
}

void Mesh::calcProjectedPts(std::array<double, 3> proj) // project point along proj directions, azimuth is random
{
proj_pts.resize(n_pts);
 proj_y = getRandOrth(proj);
 proj_x = cross(proj_y,proj);
 proj_z = proj;
//printArray(proj_x);
//printArray(proj_y);
//printArray(proj_z);
	
	for (int i=0; i < n_pts; i++)
	{
	proj_pts[i][0] = dot(pts[i],proj_x);
	proj_pts[i][1] = dot(pts[i],proj_y);
	proj_pts[i][2] = dot(pts[i],proj);
	
	}
	if (!faces_norm.empty()){
	for (int i=0; i < n_faces; i++)
	{
	std::array<double, 3> o_n = faces_norm[i];
	faces_norm[i][0] = dot(o_n,proj_x);
	faces_norm[i][1] = dot(o_n,proj_y);
	faces_norm[i][2] = dot(o_n,proj);
	//printArray(faces_norm[i]);
	}
	}
	
}


std::array<double, 3> Mesh::computeSurfRad(double reso){
	res = reso;
	std::array<double, 3> min_coo = colMin(proj_pts);
	std::array<double, 3> max_coo = colMax(proj_pts);
	// set min_coo to 0 and divide by res
	for(int i=0; i<n_pts; i++)
	{
	proj_pts[i][0] = (proj_pts[i][0]-min_coo[0]);
	proj_pts[i][1] = (proj_pts[i][1]-min_coo[1]);
	}
	n_x = ceil((max_coo[0] -min_coo[0])/res);
	n_y = ceil((max_coo[1] -min_coo[1])/res);
	//std::vector<int> idx_face(n_px); // uncomment in case we want to go back
	//std::fill(idx_face.begin(), idx_face.end(), -1);
	int n_px = n_x*n_y;
	z_pix.resize(n_px);
	std::fill(z_pix.begin(), z_pix.end(), -std::numeric_limits<double>::infinity());
	fx.resize(n_px);
	std::fill(fx.begin(), fx.end(), 0);
	fy.resize(n_px);
	std::fill(fy.begin(), fy.end(), 0);
	fz.resize(n_px);
	std::fill(fz.begin(), fz.end(), 0);
	
	for ( int i = 0 ; i < n_faces ; i++ ) 
	{
		if (faces_norm[i][2] < 0) //face not pointing the emitter
		{continue;}
		// get idx inside and z
		std::array<double, 3> pt1 = proj_pts[faces[i][0]];
		std::array<double, 3> pt2 = proj_pts[faces[i][1]];
		std::array<double, 3> pt3 = proj_pts[faces[i][2]];
		
	
		// for the idx if the value is more than z  compute the emission force and put it in the force
		std::vector<int> tt_idx;
		std::vector<double> z_val;
		getContainedIdx(i, tt_idx, z_val);
		
		
		
		// compute the force value for the face
		double rho =  colors[i][0];0.7;// props[i][0];
		double d = colors[i][0];;//props[i][1];
		std::array<double, 3> f_norm = faces_norm[i];
		double cosb = f_norm[2];
		std::array<double, 3> df;
		df[0] = (2/3*d + 2*rho*cosb) * cosb * f_norm[0];
		df[1] = (2/3*d + 2*rho*cosb) * cosb * f_norm[1];
		df[2] = (2/3*d + 2*rho*cosb) * cosb * f_norm[2]+ (1-rho) * cosb;
		//printArray(df);
		int n_tt = tt_idx.size();
		for (int j =0; j < n_tt; j++){
			int c_idx = tt_idx[j];
			//std::cout << z_val[c_idx] << " ";
			if (z_pix[c_idx] < z_val[j]) // if it is in front
			{
				//std::cout << c_idx<< " ";
				fx[c_idx] = df[0];
				//std::cout << fx[c_idx]<< " "<<df[0]<< " ";
				fy[c_idx] = df[1];
				fz[c_idx] = df[2];
				z_pix[c_idx] = z_val[j];
			
			}
		}
	}
	std::array<double, 3> F = {0,0,0};
	/*
	// sum all df
	
	for (int i = 0; i < n_px; i++){
		//std::cout << fx[i]<<" ";
		F[0] = F[0] + fx[i];
		F[1] = F[1] + fy[i];
		F[2] = F[2] + fz[i];
	} 
	*/

    F[0] = std::accumulate(fx.begin(), fx.end(), 0.0);
    F[1] = std::accumulate(fy.begin(), fy.end(), 0.0);
    F[2] = std::accumulate(fz.begin(), fz.end(), 0.0);
	//std::cout<<"n_px" << n_px <<"    ";
	F[0] = F[0]*pow(res,2);
	F[1] = F[1]*pow(res,2);
	F[2] = F[2]*pow(res,2);
	std::array<double, 3> tproj_x = {proj_x[0], proj_y[0], proj_z[0]};
	std::array<double, 3> tproj_y = {proj_x[1], proj_y[1], proj_z[1]};
	std::array<double, 3> tproj_z = {proj_x[2], proj_y[2], proj_z[2]};
	F[0] = dot(F,tproj_x);
	F[1] = dot(F,tproj_y);
	F[2] = dot(F,tproj_z);
	return F;
	}

void Mesh::getContainedIdx(int face_idx, std::vector<int>& idx, std::vector<double>& z_val){
	std::array<double, 3> pt1 = proj_pts[faces[face_idx][0]];
	std::array<double, 3> pt2 = proj_pts[faces[face_idx][1]];
	std::array<double, 3> pt3 = proj_pts[faces[face_idx][2]];
	// gettng idx of the bbbox
	double min_x = std::min(std::min(pt1[0],pt2[0]),pt3[0]);
	double max_x = std::max(std::max(pt1[0],pt2[0]),pt3[0]);
	
	double min_y = std::min(std::min(pt1[1],pt2[1]),pt3[1]);
	double max_y = std::max(std::max(pt1[1],pt2[1]),pt3[1]);
	int n_tt_x = ceil((max_x - min_x)/res);
	int n_tt_y = ceil((max_y - min_y)/res);
	int min_x_i = floor((min_x)/res);
	int min_y_i = floor((min_y)/res);
	// getting slope and intercept of the three segment
	std::array<double, 3> f_norm = faces_norm[face_idx];
	
	double a = -f_norm[0] / f_norm[2];
	double b = -f_norm[1] / f_norm[2];
	double c = dot(f_norm,pt1)/f_norm[2];// z = a*x + b*y +c 
	//printf("z = %f x + %f y + %f\n",a,b,c);
	for (int i = 0; i < n_tt_x; i++)
	{
		for (int j = 0; j < n_tt_y; j++){
			int ii = min_x_i+i;
			int jj = min_y_i+j;
			double x = static_cast<double>(ii)*res;
			double y = static_cast<double>(jj)*res;
			std::array<double, 3> tt = {x,y, 0};
			if (isLeft(pt1,pt2,tt) > -EPS && isLeft(pt2,pt3,tt) > -EPS && isLeft(pt3,pt1,tt) > -EPS){ // if is contained
				idx.push_back(jj*n_x+ii);
				// get the z value
				z_val.push_back(a*x +b*y +c);
			}
		}
	}
	
	}

void Mesh::computeFaceNorm(){ //WARNING compute the normal based on the prjocted points
	//std::cout << n_faces <<"\n";
	if (faces_norm.empty()){
	faces_norm.resize(n_faces);
	for (int i = 0; i < n_faces; i++){
		std::array<double, 3> pt1 = proj_pts[faces[i][0]];
		std::array<double, 3> pt2 = proj_pts[faces[i][1]];
		std::array<double, 3> pt3 = proj_pts[faces[i][2]];
		std::array<double, 3> a = minus(pt2,pt1);
		std::array<double, 3> b = minus(pt3,pt1);
		std::array<double, 3> n = cross(a,b);
		normalize(n);
		faces_norm[i] = n;
		//std::cout << dot(n,faces_norm[i]) <<"   ";
	}
	}
	
}
    void Mesh::dumpZVal(std::string fname){
    std::ofstream myfile;
  myfile.open(fname);
  for (int i = 0; i < n_x; i++){
  	for (int j = 0; j < n_y; j++){
  	myfile <<std::setprecision(6) << z_pix[j*n_x+i] << " ";
  	}
  	myfile << "\n";
  	}
  myfile.close();
  };
    void Mesh::dumpF(std::string fname){};




