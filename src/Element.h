#pragma once
#include<iostream>
#include<array>
#define zbX 0
#define zbY 1
#define zbZ 2
class Node {

public:
	Node();

	std::array<double, 3> zb ;
	int unknown_uvw[3];
	int unknown_abc[3];
	double von_Mises;//

	int unknown_T;
	int unknown_q;
	int ver;
	friend class Element;
};

class Face {
public:
	Face();
	int opp[2]; //opp[0] Stores the global number of a common coplanar tetrahedron    opp[1] Stores the local number of an element
	int opp_stru[2];
	int opp_T_q[2];
	double sigma_face;
	double N_Vector[3];
	double Area;
	friend  class Element;
	int boundInfo;
	int oppo_material;  //material of opposite face
	bool whether_boundary;// the face is on the boundary(1) or not(0)
};

class Element
{
	
public:
	int Global_num;                //Global Number of the element
	int domain;
	Node node[4];
	Face face[4];
	//general parameter
	double a[4];
	double b[4]; 
	double c[4];
	double d[4];
	double Ve;
	double length[12];

	double nu;
	double E;
	double alpha;
	double Q;


	double Cp;
	double Rho;
	double Kh[3];


	double sigma_temp;
	double sigma;
	double epsil;
	double eta;
	int Eedge_GBNO[24];           //Global numbering of edges
	int Hedge_GBNO[24];

	int Eedge_T[10];
	int Eedge_q[10];

	int Material;
	void InitalElement();
	//void ElementInit();
	~Element();
	Element();
private:

};


