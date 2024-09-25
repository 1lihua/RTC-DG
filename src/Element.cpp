#include "Element.h"

//node constructor for class
Node::Node() {

	//std::cout << "Node::Node()" << std::endl;
	ver = 0;
	zb[0] = 0.0;
	zb[1] = 0.0;
	zb[2] = 0.0;
	unknown_T = 0;
	unknown_q = 0;
	unknown_uvw[0] = 0;
	unknown_uvw[1] = 0;
	unknown_uvw[2] = 0;

	unknown_abc[0] = 0;
	unknown_abc[1] = 0;
	unknown_abc[2] = 0;
	von_Mises = 0;//
	//unknown_T = 0;

}

Element::Element() {
	for (size_t i = 0; i < 4; ++i) {
		node[i]; face[i];
	}
	Global_num = 0;
	for (size_t i = 0; i < 4; ++i) {
		a[i] = 0; b[i] = 0; c[i] = 0; d[i] = 0;

	}
	domain = 0;    Ve = 0;	epsil = 0; eta = 0; Material = 0; sigma = 0; sigma_temp = 0;
	for (size_t i = 0; i < 12; ++i) { length[i] = 0; }
	for (size_t i = 0; i < 24; ++i) {
		Eedge_GBNO[i] = 0; Hedge_GBNO[i] = 0;
	}
	for (size_t i = 0; i < 10; ++i) {
		Eedge_T[i] = 0; Eedge_q[i] = 0;
	}
	Q = 0.0;
	Cp = 0;
	Rho = 0;
	Kh[0] = 0;
	Kh[1] = 0;
	Kh[2] = 0;
	nu = 0;
	E = 0;
	alpha = 0;

}


void Element::InitalElement() {

	Global_num = 0;
	for (size_t i = 0; i < 4; ++i) {
		node[i].unknown_T = 0;
		node[i].unknown_q = 0;
		node[i].ver = 0;
		node[i].zb[0] = 0;
		node[i].zb[1] = 0;
		node[i].zb[2] = 0;
		face[i].oppo_material = 0;
		face[i].opp[0] = 0;
		face[i].opp[1] = 0;
		face[i].N_Vector[0] = 0;
		face[i].N_Vector[1] = 0;
		face[i].N_Vector[2] = 0;
		face[i].Area = 0;
		face[i].whether_boundary = 0;
		a[i] = 0; b[i] = 0; c[i] = 0; d[i] = 0;

	}
	domain = 0;   Ve = 0; 	epsil = 0; eta = 0; Material = 0; sigma=0; sigma_temp = 0;
	for (size_t i = 0; i < 12; ++i) { length[i] = 0; }
	for (size_t i = 0; i < 24; ++i) {
		Eedge_GBNO[i] = 0; Hedge_GBNO[i] = 0;
	}
	for (size_t i = 0; i < 10; ++i) {
		Eedge_T[i] = 0; Eedge_q[i] = 0;
	}


}
Element::~Element() {
	std::cout << "~Element()" << std::endl;

}
Face::Face() {

	for (size_t i = 0; i != 2; ++i) {
		opp[i] = 0;
		opp_stru[i] = 0;
		opp_T_q[i] = 0;
	}
	sigma_face = 0;
	boundInfo = 0;
	oppo_material = 0; opp[0] = 0; opp[1] = 0; N_Vector[0] = 0; N_Vector[1] = 0; N_Vector[2] = 0; Area = 0; whether_boundary = 0;


}

