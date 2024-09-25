#include <fstream>
#include <sstream>
#include"Global_Data.h"
#include "Element.h"
#include <Eigen>
#include <cstdlib>
#include <ctime>
#include"Partion.h"
#include"ReadFile.h"
using namespace std;
void* ReadFile() {
	Element* el = metis();
	//cout << "num_domain = " << num_domain << endl;
	return el;
}
 void BoundaryRead() {
	string filename = FaceInfo;
	string filename1 = filename + ".mphtxt";
	
	ifstream inFile;

	inFile.open(filename1);
	if (!inFile) {
		cout << "Cannot open File correctly" << endl;
		exit(EXIT_FAILURE);
	}
	inFile >> num_node_boundary;
	//cout << "num_node_boundary = " << num_node_boundary << endl;
	zb_boundary.resize(num_node_boundary, 3);
	//string str_zb;
	for (int i = 0; i < num_node_boundary; i++) {
		for (int j = 0; j < 3; j++) {
			inFile >> zb_boundary(i, j);//Store coordinates(x y z) of boundary nodes
		}
	}
	inFile >> num_element_boundary;
	//cout << "num_element_boundary=" << num_element_boundary << endl;
	Bondver = new int*[num_element_boundary];
	for (int i = 0; i < num_element_boundary; i++) {
		Bondver[i] = new int[3];
	}
	for (int i = 0; i < num_element_boundary; i++) {
		for (int j = 0; j < 3; j++) {
			inFile >> Bondver[i][j];     //zb_boundary(Bondver[i][j]-1, j),j=0,1,2
			++Bondver[i][j];
		}
	}
	boundaryInfo = new int[num_element_boundary];
	for (int i = 0; i < num_element_boundary; i++) {
		inFile >> boundaryInfo[i];   // boundary conditions  Boundary information
		++boundaryInfo[i];
	}
	inFile.close();
}