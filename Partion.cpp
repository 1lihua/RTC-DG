#include<Eigen>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include"Element.h"
#include"Global_Data.h"
#include "ReadFile.h"
using namespace Eigen;
using namespace std;
Element* metis() {
	int i, j, cc;
	string filename = BulkInfo;
	string filename1 = ".mphtxt";
	string epartname = filename + "_epart.txt";
	ifstream infile;
	infile.open(filename + filename1);
	//infile >> num_domain;
	//cout << "num_domain=" << num_domain << endl;
	//cout << "Input domain number: ";
	infile >> num_node;
	//cout << "num_node=" << num_node << endl;
	/***************************************************************************/
	/*                                                                         */
	/*             read the mesh file from the hard disk	                   */
	/*                                                                         */
	/***************************************************************************/
	zb.resize(num_node, 3);
	for (cc = 0; cc < num_node; cc++) {
		for (j = 0; j < 3; j++) {
			infile >> zb(cc, j);//gain the node's coordinate x, y ,z
		}
	}
	infile >> num_element;
	//cout << "num_element=" << num_element << endl;
	Element* el = (Element*)malloc(num_element * sizeof(Element));
	//for (size_t i = 0; i < num_element; ++i) {
	//	el[i].InitalElement();
	//}
	for (int i = 0; i < num_element; i++) {
		for (int j = 0; j < 4; j++) {
			infile >> el[i].node[j].ver;//el[i].node[j].ver   Store global node number
			++el[i].node[j].ver;//?
		}
	}
	for (int ith = 0; ith < num_element; ith++) {
		infile >> el[ith].Material;
	}
	infile.close();
	infile.open(epartname);
	for (int i = 0; i < num_element; i++) {
		infile >> el[i].domain;
	}
	infile.close();
#   pragma omp parallel for 	
	for (int i = 0; i < num_element; i++) {
		//el[i].Global_num = i + 1;
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 3; k++) {
				el[i].node[j].zb[k] = zb(el[i].node[j].ver - 1, k);     //zb   coordinate
			}
		}
	}

	for (int i = 0; i < num_element; i++) {
		el[i].Global_num = i + 1;
	}


	return el;
}