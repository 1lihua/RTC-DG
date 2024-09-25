

#include"Face_Vector.h"
using namespace std;


//Find the normal vector on each surface


void Face_Vector(Element* el,int num_element_subdomain) {

	array<double, 3> e1 = { 0,0,0 };
	array<double, 3> e2 = { 0,0,0 };
	array<double, 3> e3 = { 0,0,0 };
	array<double, 3> vtemp = { 0,0,0 };

	for (size_t i = 0; i != num_element_subdomain; ++i) {

		//face 1
		for (size_t j = 0; j != 3; ++j) {
			e1[j] = el[i].node[1].zb[j] - el[i].node[0].zb[j];
			e2[j] = el[i].node[2].zb[j] - el[i].node[0].zb[j];
			e3[j] = el[i].node[3].zb[j] - el[i].node[0].zb[j];
		}

		//cross product
		vtemp[0] = -(e1[1] * e2[2] - e1[2] * e2[1]);
		vtemp[1] = -(e1[2] * e2[0] - e1[0] * e2[2]);
		vtemp[2] = -(e1[0] * e2[1] - e1[1] * e2[0]);
		//if ((vtemp[0] * e3[0] + vtemp[1] * e3[1] + vtemp[2] * e3[2]) > 0) {
		//	vtemp[0] *= -1.0; vtemp[1] *= -1.0; vtemp[2] *= -1.0;
		//	cout << "1" << ' ';
		//}
		el[i].face[0].N_Vector[0] = vtemp[0] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[0].N_Vector[1] = vtemp[1] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[0].N_Vector[2] = vtemp[2] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));


		

		// Face 2
		for (size_t j = 0; j != 3; ++j) {
			e1[j] = el[i].node[2].zb[j] - el[i].node[0].zb[j];
			e2[j] = el[i].node[3].zb[j] - el[i].node[0].zb[j];
			e3[j] = el[i].node[1].zb[j] - el[i].node[0].zb[j];
		}

		//cross product
		vtemp[0] = -(e1[1] * e2[2] - e1[2] * e2[1]);
		vtemp[1] = -(e1[2] * e2[0] - e1[0] * e2[2]);
		vtemp[2] = -(e1[0] * e2[1] - e1[1] * e2[0]);
		//if ((vtemp[0] * e3[0] + vtemp[1] * e3[1] + vtemp[2] * e3[2]) > 0) {
		//	vtemp[0] *= -1.0; vtemp[1] *= -1.0; vtemp[2] *= -1.0;
		//	cout << "2"<< ' ';
		//}
		el[i].face[1].N_Vector[0] = vtemp[0] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[1].N_Vector[1] = vtemp[1] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[1].N_Vector[2] = vtemp[2] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));

		

		// Face 3
		for (size_t j = 0; j != 3; ++j) {
			e1[j] = el[i].node[3].zb[j] - el[i].node[0].zb[j];
			e2[j] = el[i].node[1].zb[j] - el[i].node[0].zb[j];
			e3[j] = el[i].node[2].zb[j] - el[i].node[0].zb[j];
		}

		//cross product
		vtemp[0] = -(e1[1] * e2[2] - e1[2] * e2[1]);
		vtemp[1] = -(e1[2] * e2[0] - e1[0] * e2[2]);
		vtemp[2] = -(e1[0] * e2[1] - e1[1] * e2[0]);
		//if ((vtemp[0] * e3[0] + vtemp[1] * e3[1] + vtemp[2] * e3[2]) > 0) {
		//	vtemp[0] *= -1.0; vtemp[1] *= -1.0; vtemp[2] *= -1.0;
		//	cout << "3" << ' ';
		//}
		el[i].face[2].N_Vector[0] = vtemp[0] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[2].N_Vector[1] = vtemp[1] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[2].N_Vector[2] = vtemp[2] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));

	

		
		//face 4
		for (size_t j = 0; j != 3; ++j) {
			e1[j] = el[i].node[2].zb[j] - el[i].node[1].zb[j];
			e2[j] = el[i].node[1].zb[j] - el[i].node[3].zb[j];
			e3[j] = el[i].node[0].zb[j] - el[i].node[1].zb[j];
		}

		//cross product
		vtemp[0] = -(e1[1] * e2[2] - e1[2] * e2[1]);
		vtemp[1] = -(e1[2] * e2[0] - e1[0] * e2[2]);
		vtemp[2] = -(e1[0] * e2[1] - e1[1] * e2[0]);
		//if ((vtemp[0] * e3[0] + vtemp[1] * e3[1] + vtemp[2] * e3[2]) > 0) {
		//	vtemp[0] *= -1.0; vtemp[1] *= -1.0; vtemp[2] *= -1.0;
		//	cout << "4" << ' ';
		//}
		el[i].face[3].N_Vector[0] = vtemp[0] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[3].N_Vector[1] = vtemp[1] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
		el[i].face[3].N_Vector[2] = vtemp[2] / sqrt(pow(vtemp[0], 2) + pow(vtemp[1], 2) + pow(vtemp[2], 2));
	}
}