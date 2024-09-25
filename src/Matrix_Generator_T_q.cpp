

#include "Matrix_Generator_T_q.h"
using namespace Eigen;
using namespace std;
unordered_map<int, double> mat_full1;
unordered_map<int, double> mat_full2;
unordered_map<int, double> ma_full;

unordered_map<long long, double> ma_full_T_q;



double Flabel_T_q(int vnode1, int vnode2);

double sign_judge(int edgeNo1, int edgeNo2);


vector<Triplet<double>> TriList_T_q;

void Direct_solution_T_q(Element* element, int myid);
void Direct_solution_T_q_2(Element* element, int myid);

//int Thermal_Solver_T_q(Element* element, Element* oppo_element, int myid) {
//	ofstream ofs;
//	cout.precision(16);
//	int material, el, ii, node_ii_loc;
//	double Q;
//	for (el = 0; el < num_element_subdomain; el++) {
//		material = element[el].Material;
//		for (ii = 0; ii < 4; ii++) {
//			node_ii_loc = element[el].node[ii].unknown_T;
//			if (material == 9) {
//				Q = 2e9*4/7500;
//				//Q = 0;
//				fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
//			}
//			else if (material == 11) {
//				Q = 2e9/1800;
//				//Q = 0;
//				fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
//			}
//			else if (material == 12) {
//				Q = 2e9 / 1800;
//				//Q = 0;
//				fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
//			}
//			else if((material<=8&& material >=2) || material==10) {
//				Q = 1e9/900;
//				//Q = 0;
//				fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
//			}
//		}
//	}
//	Direct_solution_T_q(element, myid);
//	return 0;
//}

int Thermal_Solver_T_q(Element* element, Element* oppo_element, int myid) {
	int material, el, ii, node_ii_loc;
	double Q;
	double Xgs, Ygs, Zgs;
	Vector3d  zb_node1, zb_node2, zb_node3, zb_node4;
	double rs[4];
	double L_node;
	double E_real[3], E_imag[3];
	double ncross_E_real[3], ncross_E_imag[3];
	double E_square;
	int node_ii1, node_ii2, signE_Ni, edgeii_loc, edgeii_E;
	double sigma;
	int count_sigma=0;
	double zb1[3], zb2[3], zb3[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	double Ftemp[6][3];
	double nx, ny, nz;
	double E_max = 0;





	for (el = 0; el < num_element_subdomain; el++) {
		material = element[el].Material;
		sigma = element[el].sigma;
		zb_node1(0) = element[el].node[0].zb[0]; zb_node1(1) = element[el].node[0].zb[1]; zb_node1(2) = element[el].node[0].zb[2];
		zb_node2(0) = element[el].node[1].zb[0]; zb_node2(1) = element[el].node[1].zb[1]; zb_node2(2) = element[el].node[1].zb[2];
		zb_node3(0) = element[el].node[2].zb[0]; zb_node3(1) = element[el].node[2].zb[1]; zb_node3(2) = element[el].node[2].zb[2];
		zb_node4(0) = element[el].node[3].zb[0]; zb_node4(1) = element[el].node[3].zb[1]; zb_node4(2) = element[el].node[3].zb[2];
		for (ii = 0; ii < 4; ii++) {
			node_ii_loc = element[el].node[ii].unknown_T;
			for (int pp = 0; pp < 5; pp++) {
				//Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0) + delta0[pp] * zb_node4(0);
				//Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1) + delta0[pp] * zb_node4(1);
				//Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2) + delta0[pp] * zb_node4(2);
				rs[0] = alpha0_3D[pp] * zb_node1(0) + beta0_3D[pp] * zb_node2(0) + gamma0_3D[pp] * zb_node3(0) + delta0_3D[pp] * zb_node4(0);
				rs[1] = alpha0_3D[pp] * zb_node1(1) + beta0_3D[pp] * zb_node2(1) + gamma0_3D[pp] * zb_node3(1) + delta0_3D[pp] * zb_node4(1);
				rs[2] = alpha0_3D[pp] * zb_node1(2) + beta0_3D[pp] * zb_node2(2) + gamma0_3D[pp] * zb_node3(2) + delta0_3D[pp] * zb_node4(2);
				L_node = (element[el].a[ii] + element[el].b[ii] * rs[0] + element[el].c[ii] * rs[1] + element[el].d[ii] * rs[2])/(6* element[el].Ve);
				E_real[0] = 0; E_real[1] = 0; E_real[2] = 0;
				E_imag[0] = 0; E_imag[1] = 0; E_imag[2] = 0;
				//E_square = 0;
				for (int edgeii_loc = 1; edgeii_loc < 13; edgeii_loc++) {
					node_ii1 = edge_node_local[edgeii_loc-1][0]; node_ii2 = edge_node_local[edgeii_loc-1][1];
					edgeii_E = element[el].Eedge_GBNO[edgeii_loc-1];
					if (edgeii_E != 0) {
						signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc - 1 + 12];
						Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
						Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
						Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
						Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
						Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
						Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
						E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
						E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					}
					//signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc-1 + 12];
					//Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
					//Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
					//Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
					//Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
					//Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
					//Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
					//E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					//E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
				}
				E_square = E_real[0] * E_real[0] + E_real[1] * E_real[1] + E_real[2] * E_real[2] + E_imag[0] * E_imag[0] + E_imag[1] * E_imag[1] + E_imag[2] * E_imag[2];
				if (node_ii_loc != 0) {
					fbr_T_q1[node_ii_loc - 1] += 1.0 * element[el].Ve * L_node * E_square * sigma * weight03_3D[pp]*0.5;
				}
				//fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve * L_node * E_square * sigma * weight0[pp];
			}
		}
	}
	ofstream ofs;
	cout.precision(16);
	for (el = 0; el < num_element_subdomain; el++) {
		material = element[el].Material;
		for (ii = 0; ii < 4; ii++) {
			node_ii_loc = element[el].node[ii].unknown_T;
			Q = element[el].Q;
			//Q = 0;
			fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;


			//if (material == 2) {
			//	Q = 1e5;
			//	//Q = 0;
			//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
			//}
			//if (material == 111) {
			//	Q = 5e7;
			//	//Q = 0;
			//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
			//}
			//if (material == 3) {
			//	Q = 2.5e7;
			//	//Q = 0;
			//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
			//}
			//if (material == 103) {
			//	Q = 1.25e7;
			//	//Q = 0;
			//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
			//}
		}
	}
	//cout << "QQQ" << endl;
	Direct_solution_T_q(element, myid);
	return 0;
}

int Thermal_Solver_T_q2(Element* element, Element* oppo_element, int myid) {
	int material, el, ii, node_ii_loc;
	double Q;
	double Xgs, Ygs, Zgs;
	Vector3d  zb_node1, zb_node2, zb_node3, zb_node4;
	double rs[4];
	double L_node;
	double E_real[3], E_imag[3];
	double ncross_E_real[3], ncross_E_imag[3];
	double E_square;
	int node_ii1, node_ii2, signE_Ni, edgeii_loc, edgeii_E;
	double sigma;
	int count_sigma = 0;
	double zb1[3], zb2[3], zb3[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	double Ftemp[6][3];
	double nx, ny, nz;
	double E_max = 0;


	ofstream ofs;
	cout.precision(16);
	for (el = 0; el < num_element_subdomain; el++) {
		material = element[el].Material;
		for (ii = 0; ii < 4; ii++) {
			node_ii_loc = element[el].node[ii].unknown_T;
			Q = element[el].Q;
			//Q = 0;
			fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;

		}
	}
	//cout << "QQQ" << endl;
	Direct_solution_T_q(element, myid);
	return 0;
}


int Thermal_Solver_T_q_2(Element* element, Element* oppo_element, int myid) {

	double weight0[5] = { {-0.8},{0.45},{0.45},{0.45},{0.45} };
	double alpha0[5] = { {0.25},{0.50000000000000},{0.166666666666667},{0.166666666666667},{0.166666666666667} };
	double beta0[5] = { {0.25},{0.166666666666667},{0.50000000000000},{0.166666666666667},{0.166666666666667} };
	double gamma0[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.50000000000000},{0.166666666666667} };
	double delta0[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.166666666666667},{0.50000000000000} };


	int material, el, ii, node_ii_loc;
	double Q;
	double Xgs, Ygs, Zgs;
	Vector3d  zb_node1, zb_node2, zb_node3, zb_node4;
	double rs[4];
	double L_node;
	double E_real[3], E_imag[3];
	double E_square;
	int node_ii1, node_ii2, signE_Ni, edgeii_loc, edgeii_E;
	double sigma;
	int count_sigma = 0;
	double zb1[3], zb2[3], zb3[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];


	for (el = 0; el < num_element_subdomain; el++) {
		material = element[el].Material;
		sigma = element[el].sigma;
		zb_node1(0) = element[el].node[0].zb[0]; zb_node1(1) = element[el].node[0].zb[1]; zb_node1(2) = element[el].node[0].zb[2];
		zb_node2(0) = element[el].node[1].zb[0]; zb_node2(1) = element[el].node[1].zb[1]; zb_node2(2) = element[el].node[1].zb[2];
		zb_node3(0) = element[el].node[2].zb[0]; zb_node3(1) = element[el].node[2].zb[1]; zb_node3(2) = element[el].node[2].zb[2];
		zb_node4(0) = element[el].node[3].zb[0]; zb_node4(1) = element[el].node[3].zb[1]; zb_node4(2) = element[el].node[3].zb[2];
		for (ii = 0; ii < 4; ii++) {
			node_ii_loc = element[el].Eedge_T[ii];
			for (int pp = 0; pp < 5; pp++) {
				//Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0) + delta0[pp] * zb_node4(0);
				//Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1) + delta0[pp] * zb_node4(1);
				//Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2) + delta0[pp] * zb_node4(2);
				rs[0] = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0) + delta0[pp] * zb_node4(0);
				rs[1] = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1) + delta0[pp] * zb_node4(1);
				rs[2] = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2) + delta0[pp] * zb_node4(2);
				L_node = (2.0*(element[el].a[ii] + element[el].b[ii] * rs[0] + element[el].c[ii] * rs[1] + element[el].d[ii] * rs[2])* (element[el].a[ii] + element[el].b[ii] * rs[0] + element[el].c[ii] * rs[1] + element[el].d[ii] * rs[2])/(6 * element[el].Ve) - (element[el].a[ii] + element[el].b[ii] * rs[0] + element[el].c[ii] * rs[1] + element[el].d[ii] * rs[2])) / (6 * element[el].Ve);
				E_real[0] = 0; E_real[1] = 0; E_real[2] = 0;
				E_imag[0] = 0; E_imag[1] = 0; E_imag[2] = 0;
				//E_square = 0;
				for (int edgeii_loc = 1; edgeii_loc < 13; edgeii_loc++) {
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					edgeii_E = element[el].Eedge_GBNO[edgeii_loc - 1];
					if (edgeii_E != 0) {
						signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc - 1 + 12];
						Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
						Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
						Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
						Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
						Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
						Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
						E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
						E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					}
					//signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc-1 + 12];
					//Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
					//Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
					//Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
					//Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
					//Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
					//Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
					//E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					//E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
				}
				E_square = E_real[0] * E_real[0] + E_real[1] * E_real[1] + E_real[2] * E_real[2] + E_imag[0] * E_imag[0] + E_imag[1] * E_imag[1] + E_imag[2] * E_imag[2];
				if (node_ii_loc != 0) {
					fbr_T_q1[node_ii_loc - 1] += 1.0 * element[el].Ve * L_node * E_square * sigma * weight0[pp] * 0.5;
				}
				//fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve * L_node * E_square * sigma * weight0[pp];
			}
		}
		for (ii = 0; ii < 6; ii++) {
			node_ii_loc = element[el].Eedge_T[ii+4];
			for (int pp = 0; pp < 5; pp++) {
				//Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0) + delta0[pp] * zb_node4(0);
				//Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1) + delta0[pp] * zb_node4(1);
				//Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2) + delta0[pp] * zb_node4(2);
				rs[0] = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0) + delta0[pp] * zb_node4(0);
				rs[1] = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1) + delta0[pp] * zb_node4(1);
				rs[2] = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2) + delta0[pp] * zb_node4(2);
				int node_ii11 = edge_node_local[ii][0] - 1; int node_ii22 = edge_node_local[ii][1] - 1;
				L_node = (element[el].a[node_ii11] + element[el].b[node_ii11] * rs[0] + element[el].c[node_ii11] * rs[1] + element[el].d[node_ii11] * rs[2]) / (6 * element[el].Ve)*(element[el].a[node_ii22] + element[el].b[node_ii22] * rs[0] + element[el].c[node_ii22] * rs[1] + element[el].d[node_ii22] * rs[2]) / (6 * element[el].Ve);
				E_real[0] = 0; E_real[1] = 0; E_real[2] = 0;
				E_imag[0] = 0; E_imag[1] = 0; E_imag[2] = 0;
				//E_square = 0;
				for (int edgeii_loc = 1; edgeii_loc < 13; edgeii_loc++) {
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					edgeii_E = element[el].Eedge_GBNO[edgeii_loc - 1];
					if (edgeii_E != 0) {
						signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc - 1 + 12];
						Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
						Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
						Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
						Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
						Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
						Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
						E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
						E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
						E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
						E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					}
					//signE_Ni = 1.0 * element[el].Eedge_GBNO[edgeii_loc-1 + 12];
					//Li1[0] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].b[node_ii2 - 1];
					//Li1[1] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].c[node_ii2 - 1];
					//Li1[2] = element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii1 - 1] + element[el].b[node_ii1 - 1] * rs[0] + element[el].c[node_ii1 - 1] * rs[1] + element[el].d[node_ii1 - 1] * rs[2]) * element[el].d[node_ii2 - 1];
					//Li2[0] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].b[node_ii1 - 1];
					//Li2[1] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].c[node_ii1 - 1];
					//Li2[2] = sign_judge(edgeii_loc, 6) * element[el].length[edgeii_loc - 1] / (36.0 * element[el].Ve * element[el].Ve) * (element[el].a[node_ii2 - 1] + element[el].b[node_ii2 - 1] * rs[0] + element[el].c[node_ii2 - 1] * rs[1] + element[el].d[node_ii2 - 1] * rs[2]) * element[el].d[node_ii1 - 1];
					//E_real[0] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_real[1] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_real[2] += Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
					//E_imag[0] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[0] + Li2[0]);
					//E_imag[1] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[1] + Li2[1]);
					//E_imag[2] += Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]);
				}
				E_square = E_real[0] * E_real[0] + E_real[1] * E_real[1] + E_real[2] * E_real[2] + E_imag[0] * E_imag[0] + E_imag[1] * E_imag[1] + E_imag[2] * E_imag[2];
				if (node_ii_loc != 0) {
					fbr_T_q1[node_ii_loc - 1] += 1.0 * element[el].Ve * L_node * E_square * sigma * weight0[pp] * 0.5;
				}
				//fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve * L_node * E_square * sigma * weight0[pp];
			}
		}
	}


	////cout << "count_sigma = " << count_sigma << endl;
	//cout << "num_element_subdomain = " << num_element_subdomain << endl;
	////ofstream ofs;
	////cout.precision(16);
	////int material, el, ii, node_ii_loc;
	////double Q;
	//for (el = 0; el < num_element_subdomain; el++) {
	//	material = element[el].Material;
	//	for (ii = 0; ii < 10; ii++) {
	//		node_ii_loc = element[el].Eedge_T[ii];
	//		if (material == 2) {
	//			Q = 1e5;
	//			//Q = 0;
	//			if (ii < 4) {
	//				fbr_T_q1[node_ii_loc - 1] -= 1.0 / 20 * element[el].Ve  * Q;
	//			}
	//			else {
	//				fbr_T_q1[node_ii_loc - 1] += 1.0 / 5 * element[el].Ve  * Q;
	//			}
	//			//fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
	//		}
	//		//if (material == 111) {
	//		//	Q = 5e7;
	//		//	//Q = 0;
	//		//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
	//		//}
	//		//if (material == 3) {
	//		//	Q = 2.5e7;
	//		//	//Q = 0;
	//		//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
	//		//}
	//		//if (material == 103) {
	//		//	Q = 1.25e7;
	//		//	//Q = 0;
	//		//	fbr_T_q1[node_ii_loc - 1] += 1 * element[el].Ve / 4 * Q;
	//		//}
	//	}
	//}
	//cout << "QQQ111" << endl;


	Direct_solution_T_q_2(element, myid);
	return 0;
}



void* Matrix_Generator1_T_q(int myid) {
	int unknown_dm, row1, nzero1, nn1, nn2, nn3, final_row, nnz_dm, tot_dm, nrhs, cnt_n1, end1, start1, cnt_n2, end2, start2, final_col, mm2;
	double* bb, * x, * a_pro, * b_pro, * pro;
	int maxfct, mnum, mtype, phase1, phase2, phase3, n, msglvl, error, comm, info;
	int iparm[64]{ 0 };
	int pt[64]{ 0 };
	int* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
	iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;

	for (int n1 = myid; n1 < myid+1; n1++) {
		//unknown_dm = num_unKnown_subdomain(n1, solver_type1);
		unknown_dm = num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		nnz_dm = num_nzero_Pmatrix_T_q[n1];
		row1 = 0;
		nzero1 = 0;
		//nn3 = n3_dm[n1];
		//nn2 = nn3_dm[n1];
		nn2 = 0;
		final_row = r_dm_T_q[n1];
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 };
		mkl_dcsrcoo(job, &unknown_dm, p_pardiso_T_q , jp_T_q, ip_T_q, &nnz_dm, acoo_T_q , rowind1 , colind1, &info);


		cnt_n1 = Vertice_Boundary_T_q[n1];

		bb = new double[cnt_n1 * unknown_dm]();
		x = new double[cnt_n1 * unknown_dm];
		nrhs = cnt_n1;

		for (int i = 0; i < cnt_n1; i++) {
			bb[i * unknown_dm + i + unknown_dm - cnt_n1] = 1.0;//Identity matrix
		}

		tot_dm = unknown_dm * unknown_dm;
		//bb = new double[tot_dm]();
		//x = new double[tot_dm];
		//nrhs = unknown_dm;
		//for (int ith = 0; ith < unknown_dm; ith++) {
		//	bb[ith * (unknown_dm + 1)] = 1;
		//}
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_T_q , ip_T_q, jp_T_q, perm2, &nrhs, iparm, &msglvl, bb, x, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_T_q , ip_T_q, jp_T_q, perm2, &nrhs, iparm, &msglvl, bb, x, &error);
		delete[] bb;
		//if (myid == 0) {
		//	ofstream ofs12153("inv_T_q1.csv");
		//	for (int j = 0; j < tot_dm; j++) {
		//		ofs12153 << x[j] << endl;
		//	}

		//}
		//cnt_n1 = Vertice_Boundary_T_q[n1];
		
		//end1 = num_unKnown_subdomain(n1, solver_type1);
		end1 = num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];

		start1 = end1 - cnt_n1 + 1;
		a_pro = new double[cnt_n1 * cnt_n1]();
//#		pragma omp parallel for num_threads(8)

//#		pragma omp parallel for 
//		for (int i = 0; i < tot_dm; i++) {
//			if (i % unknown_dm + 1 >= start1 && i % unknown_dm + 1 <= end1 && i / unknown_dm + 1 >= start1 && i / unknown_dm + 1 <= end1) {
//				a_pro[cnt_n1 * (i % unknown_dm + 1 - start1) + i / unknown_dm + 1 - start1] = x[i];
//			}
//		}
		
#		pragma omp parallel for 
		for (int i = 0; i < cnt_n1 * unknown_dm; i++) {
			if (i % unknown_dm + 1 >= start1 && i % unknown_dm + 1 <= end1 && i / unknown_dm + 1 >= 1 && i / unknown_dm + 1 <= cnt_n1) {
				a_pro[cnt_n1 * (i % unknown_dm + 1 - start1) + i / unknown_dm + 1 - 1] = x[i];
			}
		}
		delete[] x;

		double* fbr_temp = new double[unknown_dm];
		nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_T_q, ip_T_q, jp_T_q, perm2, &nrhs, iparm, &msglvl, fbr_T_q1, fbr_temp, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_T_q, ip_T_q, jp_T_q, perm2, &nrhs, iparm, &msglvl, fbr_T_q1, fbr_temp, &error);

		for (int i = 0; i < cnt_n1; i++) {
			fbb_T_q_G[i] = fbr_temp[i + unknown_dm - cnt_n1];
		}
		delete[] fbr_temp;



		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_C(n1, n2) != 0) {
				mm2 = nn2 + nnz_C(n1, n2);
				if (n1 != n2) {
					cnt_n2 = Vertice_Boundary_T_q[n2];
					end2 = num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];

					//cnt_n2 = Node_Boundary_Thermal[n2];
					//end2 = num_unKnown_subdomain(n2, solver_type1);
					start2 = end2 - cnt_n2 + 1;
					b_pro = new double[cnt_n1 * cnt_n2]();
					pro = new double[cnt_n1 * cnt_n2]();
					for (int i = nn2; i < mm2; i++) {
						if (mrow_T_q[i] >= start1 && mrow_T_q[i] <= end1 && mcol_T_q[i] >= start2 && mcol_T_q[i] <= end2) {
							b_pro[cnt_n2 * (mrow_T_q[i] - start1) + mcol_T_q[i] - start2] = m_T_q[i];
						}
					}

					//if (myid == 0) {
					//	ofstream ofs12154("b_pro.csv");
					//	for (int j = 0; j < cnt_n1 * cnt_n2; j++) {
					//		ofs12154 << b_pro[j] << endl;
					//	}

					//}



					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, 1, a_pro, cnt_n1, b_pro, cnt_n2, 0, pro, cnt_n2);

					//if (myid == 0) {
					//	ofstream ofs12155("bpro.csv");
					//	for (int j = 0; j < cnt_n1 * cnt_n2; j++) {
					//		ofs12155 << pro[j] << endl;
					//	}

					//}
					vector<Triplet<double>> Tri_tmp(cnt_n1 * cnt_n2);
					int List_step = 0;
					for (int i = 0; i < cnt_n1 * cnt_n2; i++) {
						if (pro[i] != 0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<double>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					TriList_T_q.insert(TriList_T_q.end(), Tri_tmp.begin(), Tri_tmp.end());
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Vertice_Boundary_T_q[n2];
			//final_col += Node_Boundary_Thermal[n2];
		}
		delete[] a_pro;


//		a_pro = new double[unknown_dm * cnt_n1];
////#		pragma omp parallel for num_threads(8)
//#		pragma omp parallel for 
//		for (int i = 0; i < tot_dm; i++) {
//			if (i % unknown_dm + 1 >= start1 && i % unknown_dm + 1 <= end1) {
//				a_pro[unknown_dm * (i % unknown_dm + 1 - start1) + i / unknown_dm + 1 - 1] = x[i];
//			}
//		}
//		pro = new double[cnt_n1];
//		double* fbr_tmp = fbr_T_q.data();
//		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, unknown_dm, 1, a_pro, unknown_dm, fbr_T_q1, 1, 0, pro, 1);
//		for (int i = 0; i < cnt_n1; i++) {
//			fbb[i] = pro[i];
//		}
//
//		delete[] a_pro;
//		delete[] pro;

		//if (myid == 0) {
		//	ofstream ofs12155("fbr_T_q1.csv");
		//	for (int j = 0; j < unknown_dm; j++) {
		//		ofs12155 << fbr_T_q1[j] << endl;
		//	}

		//}



		//delete[] invA;
		//delete[] rowA;
		//delete[] colA;
	}
}

void* Solver_T_q(int myid) {
	int maxfct, mnum, mtype, phase1, phase2, n, nrhs, msglvl, error;
	int* perm = nullptr;;
	int iparm[64];
	int pt[64];
	maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; nrhs = 1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 2;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;
	for (int n1 = myid; n1 < myid+1; n1++) {
		//int unknown_dm = num_unKnown_subdomain(n1, solver_type1);
		int unknown_dm = num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		int nnz_dm = num_nzero_Pmatrix_T_q[n1];
		//int tot_dm = unknown_dm * unknown_dm;
		int cnt_n2 = 0, mm3, nn3 = nn3_dm_T_q[n1], unknown_dm2;
		int nn4 = 0, mm4;
		int nn2 = 0;
		int nzero1 = 0, row1 = 0;
		//double* b_pro, * pro, * fb_tmp, * fbb = fbr_T_q.data();
		int cnt_n1 = 0;
		double* b_pro, * pro, * fb_tmp, * fbb ;
		fbb = new double[unknown_dm];
		for (int i = 0; i < unknown_dm; i++) {
			fbb[nn2 + i] = fbr_T_q1[i];
		}
		cnt_n1 = Vertice_Boundary_T_q[n1];
		for (int n2 = 0; n2 < num_domain; n2++) {
			//cnt_n2 = Node_Boundary_Thermal(n2);
			cnt_n2 = Vertice_Boundary_T_q[n2];
			
			mm4 = nn4 + cnt_n2;
			if (nnz_C(n1, n2) != 0) {
				mm3 = nn3 + nnz_C(n1, n2);
				//unknown_dm2 = num_unKnown_subdomain(n2, solver_type1);
				unknown_dm2 = num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];

				int end1, start1, end2, start2;
				start1 = unknown_dm - cnt_n1 + 1;

				start2 = unknown_dm2 - cnt_n2 + 1;
				end1 = unknown_dm;

				end2 = unknown_dm2;


				if (n1 != n2) {
					b_pro = new double[cnt_n1 * cnt_n2]();
					pro = new double[cnt_n1];
					fb_tmp = new double[cnt_n2]();
					for (int i = 0; i < cnt_n1 * cnt_n2; i++) {
						b_pro[i] = 0;
					}
					//b_pro = new double[unknown_dm * unknown_dm2]();
					//pro = new double[unknown_dm];
					//fb_tmp = new double[unknown_dm2]();


					for (int i = 0; i < cnt_n2; i++) {
						fb_tmp[i] = xx[nn4 + i];
					}
					for (int i = nn3; i < mm3; i++) {
						if (mrow_T_q[i] >= start1 && mrow_T_q[i] <= end1 && mcol_T_q[i] >= start2 && mcol_T_q[i] <= end2) {

							b_pro[cnt_n2 * (mrow_T_q[i] - start1) + mcol_T_q[i] - start2] = m_T_q[i];
						}


						//b_pro[unknown_dm2 * (mrow[i] - 1) + mcol[i] - 1] = m_T_q[i];
					}
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, 1, b_pro, cnt_n2, fb_tmp, 1, 0, pro, 1);
					for (int i = 0; i < cnt_n1; i++) {
						fbb[nn2 + i+ unknown_dm- cnt_n1] -= pro[i];
					}
					delete[] b_pro;
					delete[] pro;
					delete[] fb_tmp;

				}
				nn3 = mm3;
			}
			nn4 = mm4;
		}
		int nrhs = 1;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_T_q, ip_T_q, jp_T_q, perm, &nrhs, iparm, &msglvl, fbb + nn2, Xh + nn2, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_T_q, ip_T_q, jp_T_q, perm, &nrhs, iparm, &msglvl, fbb + nn2, Xh + nn2, &error);
	}
	return 0;
}

void Direct_solution_T_q(Element* element,int myid) {
	//int info, nzero1, row1, nn1, nn2, nn3, unknown_dm, nnz_dm;

	fbb_T_q_G = new double[Vertice_Boundary_T_q[myid]]();


	//cout << "Vertice_Boundary_T_q[" << myid << "] is " << Vertice_Boundary_T_q[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);
	
	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	Matrix_Generator1_T_q(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	//cout << "myid is " << myid << "  Volumn solve time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);




	int mat_size = TriList_T_q.size();   //  nnz of Zi*Cij (j=1,2,3...,n)

	//cout << "myid is " << myid << "  mat_size is " << mat_size << endl;

	double* m_final = new double[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i] = (TriList_T_q[i].value());
		
		rowm[i] = TriList_T_q[i].row() + 1;
		colm[i] = TriList_T_q[i].col() + 1;
	}
	TriList_T_q.clear();
	//cout << " my id is " << myid << endl;
	cout << " mat_size is " << mat_size << endl;
	int* nnz_ZC_process = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		nnz_ZC_process[j] = 0;
	}
	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	//cout << "after,  my id is " << myid << endl;
	int face_mat_size = 0;
	for (int j = 0; j < num_domain; j++) {
		//cout << nnz_ZC_process[j] << endl;
		face_mat_size += nnz_ZC_process[j] + Vertice_Boundary_T_q[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		//cout << "face_mat_size_T_q = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Vertice_Boundary_T_q[j];
	}
	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	double* m_final_total = nullptr;
	double* fbb_total = nullptr;

	if (myid == 0) {
		m_final_total = new double[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new double[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j] = 0.0;
		}
	}
	//m_final_total = new double[face_mat_size];
	//rowm_total = new int[face_mat_size];
	//colm_total = new int[face_mat_size];
	//fbb_total = new double[num_unk_boundary];
	//for (int j = 0; j < num_unk_boundary; j++) {
	//	fbb_total[j] = 0.0;
	//}

	double* fbb_temp = new double[Vertice_Boundary_T_q[myid]];
	for (int j = 0; j < Vertice_Boundary_T_q[myid]; j++) {
		fbb_temp[j]= fbb_T_q_G[j];
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Vertice_Boundary_T_q[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Vertice_Boundary_T_q[myid], MPI_DOUBLE, fbb_total, Vertice_Boundary_T_q, address_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;

	time_t start, end;

	start = time(NULL);
	xx = new double[num_unk_boundary];
	for (int i = 0; i < num_unk_boundary; i++) {
		xx[i] = 0.0;
	}


	//omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j]=1.0;
		}
		ofstream ofs11141("m_final_total.txt");
		ofstream ofs11142("rowm_total.txt");
		ofstream ofs11143("colm_total.txt");
		for (int j = 0; j < 100; j++) {
			ofs11141 << m_final_total[j]<< endl;
			ofs11142 << rowm_total[j] << endl;
			ofs11143 << colm_total[j] << endl;
		}
		//delete[] m_final_total;
		double* m_pardiso = new double[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_dcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_total, rowm_total, colm_total, &info);
		//mkl_dcsrcoo(job, &num_unKnown_eb, m_pardiso1, jm1, im1, &nnz, m_final, rowm, colm, &info);
		//delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;
		delete[] m_final_total;


		//cout << "test5555" << endl;
		MKL_INT maxfct, mnum, mtype, phase1, phase2, phase3, nrhs, msglvl, error;
		MKL_INT iparm[64]{ 0 };
		int pt[64]{ 0 };
		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
		iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;

		//iparm[0] = 0;
		//iparm[1] = 3;
		//iparm[2] = 16;
		//iparm[9] = 13;
		//iparm[10] = 1;
		//iparm[12] = 1;
		//iparm[17] = -1;
		//iparm[18] = -1;
		//iparm[1] = 3; iparm[24] = 2;


		//nzero1 = 0; row1 = 0;
		//int cnt_n1, cnt_n2, start1, start2, end1, end2, final_row = 0, final_col = 0;
		nrhs = 1;
		cout << "test8888" << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		cout << "test6666" << endl;
		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		 delete[] m_pardiso; delete[] im; delete[] jm; delete[]  fbb_total;
		//delete[]fbb_total_mkl;
		// calulate full solution Um

	}




	MPI_Barrier(MPI_COMM_WORLD);
	end = time(NULL);
	double time1 = (double)(end - start);



	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}



	Xh = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	MPI_Bcast(xx, num_unk_boundary, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	Solver_T_q(myid);
	delete[] xx;
	xx = nullptr;
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);


	if (myid == 0) {
		cout << "time is " << time1 << endl;
		ofstream ofs1119("data_DG_RTC_patch.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain_T_q[0][0] + num_unknown_subdomain_T_q[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time1 << "s" << ',' << num_nzero_Pmatrix_T_q[myid] << ',' << face_mat_size << ',' << Vertice_Boundary_T_q[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << mat_size << endl;
	}



	int un_T_temp=0;
	int* Accumulated_T_temp = new int[num_domain];
	int* unknown_T_temp = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		Accumulated_T_temp[j]=0;
		unknown_T_temp[j] = num_unknown_subdomain_T_q[j][0];

	}

	for (int j = 1; j < num_domain; j++) {
		Accumulated_T_temp[j] = Accumulated_T_temp[j-1]+ num_unknown_subdomain_T_q[j-1][0];
	}


	for (int j = 0; j < num_domain; j++) {
		un_T_temp += num_unknown_subdomain_T_q[j][0];
	}

	double* X_T = new double[un_T_temp];
	double* X_h_temp = new double[num_unknown_subdomain_T_q[myid][0]];
	int* node_global_temp = new int[num_unknown_subdomain_T_q[myid][0]];

	int* V_T = new int[un_T_temp];

	for (int j = 0; j < num_element_subdomain; j++) {
		for (int ii = 0; ii < 4; ii++) {
			if (element[j].node[ii].unknown_T != 0) {

				node_global_temp[element[j].node[ii].unknown_T - 1] = element[j].node[ii].ver;

			}
		}
	}

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0]; j++) {
		X_h_temp[j] = Xh[j];
	}

	MPI_Gatherv(X_h_temp, num_unknown_subdomain_T_q[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp, num_unknown_subdomain_T_q[myid][0], MPI_INT, V_T, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);

	delete[]Accumulated_T_temp; delete[]unknown_T_temp;  delete[]node_global_temp; delete[]X_h_temp; //delete[]Xh;




	double* XX1 = new double[num_node];
	for (int i = 0; i < num_node; i++) {
		XX1[i] = T0;
	}
	if (myid == 0) {
		for (int i = 0; i < un_T_temp; i++) {
			int temp_ver = V_T[i];
			XX1[temp_ver - 1] = X_T[i];
			
		}

		ofstream ofs12157("Patch_T_ht.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12157 << XX1[i] << endl;
		}


	}


	delete[]X_T; delete[]XX1; delete[]V_T;



}

void Direct_solution_T_q_2(Element* element, int myid) {
	int info, nzero1, row1, nn1, nn2, nn3, unknown_dm, nnz_dm;

	fbb_T_q_G = new double[Vertice_Boundary_T_q[myid]]();


	cout << "Vertice_Boundary_T_q[" << myid << "] is " << Vertice_Boundary_T_q[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	Matrix_Generator1_T_q(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	//cout << "myid is " << myid << "  Volumn solve time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);

	double time_out = (double)(end_out - start_out);



	int mat_size = TriList_T_q.size();   //  nnz of Zi*Cij (j=1,2,3...,n)

	//cout << "myid is " << myid << "  mat_size is " << mat_size << endl;

	double* m_final = new double[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i] = (TriList_T_q[i].value());

		rowm[i] = TriList_T_q[i].row() + 1;
		colm[i] = TriList_T_q[i].col() + 1;
	}
	TriList_T_q.clear();
	//cout << " my id is " << myid << endl;
	//cout << " mat_size is " << mat_size << endl;
	int* nnz_ZC_process = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		nnz_ZC_process[j] = 0;
	}
	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	//cout << "after,  my id is " << myid << endl;
	int face_mat_size = 0;
	for (int j = 0; j < num_domain; j++) {
		//cout << nnz_ZC_process[j] << endl;
		face_mat_size += nnz_ZC_process[j] + Vertice_Boundary_T_q[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size_T_q = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Vertice_Boundary_T_q[j];
	}
	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	double* m_final_total = nullptr;
	double* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new double[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new double[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j] = 0.0;
		}



	}
	double* fbb_temp = new double[Vertice_Boundary_T_q[myid]];
	for (int j = 0; j < Vertice_Boundary_T_q[myid]; j++) {
		fbb_temp[j] = fbb_T_q_G[j];
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Vertice_Boundary_T_q[j - 1];
	}

	cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Vertice_Boundary_T_q[myid], MPI_DOUBLE, fbb_total, Vertice_Boundary_T_q, address_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;

	time_t start, end;

	start = time(NULL);
	xx = new double[num_unk_boundary];
	for (int i = 0; i < num_unk_boundary; i++) {
		xx[i] = 0.0;
	}


	//omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j] = 1.0;
		}
		//m_final_mkl = new MKL_Complex16[face_mat_size];
		////#pragma omp  parallel for
		//for (int j = 0; j < face_mat_size; j++) {
		//	m_final_mkl[j].real = m_final_total[j].real();
		//	m_final_mkl[j].imag = m_final_total[j].imag();
		//}

		//delete[] m_final_total;
		double* m_pardiso = new double[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_dcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_total, rowm_total, colm_total, &info);
		//mkl_dcsrcoo(job, &num_unKnown_eb, m_pardiso1, jm1, im1, &nnz, m_final, rowm, colm, &info);
		//delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;
		delete[] m_final_total;


		//cout << "test5555" << endl;
		MKL_INT maxfct, mnum, mtype, phase1, phase2, phase3, nrhs, msglvl, error;
		MKL_INT iparm[64]{ 0 };
		int pt[64]{ 0 };
		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
		iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;
		nzero1 = 0; row1 = 0;
		int cnt_n1, cnt_n2, start1, start2, end1, end2, final_row = 0, final_col = 0;
		nrhs = 1;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		//cout << "test6666" << endl;
		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		delete[] m_pardiso; delete[] im; delete[] jm; delete[]  fbb_total;
		//delete[]fbb_total_mkl;
		// calulate full solution Um

	}


	end = time(NULL);

	double time = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	if (myid == 0) {
		cout << "time is " << time << endl;
		ofstream ofs1119("data_DG_RTC_BGA.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain_T_q[0][0] + num_unknown_subdomain_T_q[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix_T_q[myid] << ',' << face_mat_size << ',' << Vertice_Boundary_T_q[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << mat_size << endl;

		//ofstream ofs12150("face_T_q.csv");
		//for (int j = 0; j < num_unk_boundary; j++) {
		//	ofs12150 << xx[j] << endl;
		//}



	}


	Xh = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	MPI_Bcast(xx, num_unk_boundary, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	cout << "test7777" << endl;
	Solver_T_q(myid);
	delete[] xx;


	if (myid == 0) {
		ofstream ofs20240304("microstrip_line_T_ht111.txt");
		for (int i = 0; i < num_unknown_subdomain_T_q[myid][0]; i++) {
			ofs20240304 << Xh[i] << endl;
		}
	}


	cout << "test8888" << endl;
	int un_T_temp = 0;
	int* Accumulated_T_temp = new int[num_domain];
	int* unknown_T_temp = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		Accumulated_T_temp[j] = 0;
		unknown_T_temp[j] = num_unknown_subdomain_T_q[j][0];

	}

	for (int j = 1; j < num_domain; j++) {
		Accumulated_T_temp[j] = Accumulated_T_temp[j - 1] + num_unknown_subdomain_T_q[j - 1][0];
	}


	for (int j = 0; j < num_domain; j++) {
		un_T_temp += num_unknown_subdomain_T_q[j][0];
	}

	double* X_T = new double[un_T_temp];
	double* X_h_temp = new double[num_unknown_subdomain_T_q[myid][0]];
	int* node_global_temp = new int[num_unknown_subdomain_T_q[myid][0]];

	int* V_T = new int[un_T_temp];

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0]; j++) {
		node_global_temp[j] = 0;
	}

	for (int j = 0; j < num_element_subdomain; j++) {
		for (int ii = 0; ii < 4; ii++) {
			if (element[j].Eedge_T[ii] != 0) {

				node_global_temp[element[j].Eedge_T[ii] - 1] = element[j].node[ii].ver;

			}
		}
	}
	cout << "test2222" << endl;
	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0]; j++) {
		X_h_temp[j] = Xh[j];
	}

	MPI_Gatherv(X_h_temp, num_unknown_subdomain_T_q[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp, num_unknown_subdomain_T_q[myid][0], MPI_INT, V_T, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);

	delete[]Accumulated_T_temp; delete[]unknown_T_temp; delete[]Xh; delete[]node_global_temp; delete[]X_h_temp;



	cout << "test9999" << endl;
	double* XX1 = new double[num_node];
	for (int i = 0; i < num_node; i++) {
		XX1[i] = T0;
	}
	if (myid == 0) {
		for (int i = 0; i < un_T_temp; i++) {
			int temp_ver = V_T[i];
			if (temp_ver != 0) {
				XX1[temp_ver - 1] = X_T[i];
			}
			

		}

		ofstream ofs12157("microstrip_line_T_ht.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12157 << XX1[i] << endl;
		}




	}


	delete[]X_T; delete[]XX1; delete[]V_T;



}


int Matrix_Generator_T_q(Element* element, Element* oppo_element, int myid) {
	ofstream ofs;
	cout.precision(16);
	//cout << "test44" << endl;
	long long temp;
	int ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int  nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;
	
	fbr_T_q.resize(num_unKnown_Thermal);
	//fbr_T_q.resize(num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]);
	//cout << "test77" << endl;

	fbr_T_q.setZero();
	fbr_T_q1 = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]; j++) {
		fbr_T_q1[j] = 0;
	}


	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii_glb = element[el].node[ii].ver; node_ii_loc = element[el].node[ii].unknown_T;
			ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
			for (jj = 0; jj < 4; jj++) {
				node_jj_glb = element[el].node[jj].ver; node_jj_loc = element[el].node[jj].unknown_T;
				aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					ma_full_T_q[temp] += (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
				if ((node_ii_loc != 0) && (node_jj_loc == 0)) {
					fbr_T_q1[(node_ii_loc - 1)] -= Td * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
			}
		}

		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp_T_q[0]; ofn1 = element[el].face[nn].opp_T_q[1];
			if ((op1 == 0) && (ofn1 == -7)) {
				for (ii = 0; ii < 3; ii++) {
					node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
					for (jj = 0; jj < 3; jj++) {
						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 1.0 / 3 * element[el].face[nn].Area;
					
				}
			}
		}

		//double Hcvt1 = 10.0;
		//for (nn = 0; nn < 4; nn++) {
		//	op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
		//	if ((op1 == 0) && ( ofn1 == -6 || ofn1 == -66)) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
		//			for (jj = 0; jj < 3; jj++) {
		//				node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
		//				temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
		//				if (node_ii_loc != 0 && node_jj_loc != 0) {
		//					ma_full_T_q[temp] += Hcvt1 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//				if (node_ii_loc != 0 && node_jj_loc == 0) {
		//					fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt1 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//			}
		//			fbr_T_q1[(node_ii_loc - 1)] += Hcvt1 * Ta * 1.0 / 3 * element[el].face[nn].Area;
		//		}
		//	}
		//}



		//for (nn = 0; nn < 4; nn++) {
		//	op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
		//	if ((op1 == 0) && (ofn1 == -7 )) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
		//			for (jj = 0; jj < 3; jj++) {
		//				node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
		//				temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
		//				if (node_ii_loc != 0 && node_jj_loc != 0) {
		//					ma_full_T_q[temp] += 0 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//				if (node_ii_loc != 0 && node_jj_loc == 0) {
		//					fbr_T_q1[(node_ii_loc - 1)] -= Td * 0 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//			}
		//			fbr_T_q1[(node_ii_loc - 1)] += 0 * Ta * 1.0 / 3 * element[el].face[nn].Area;
		//		}
		//	}
		//}



		//for (nn = 0; nn < 4; nn++) {
		//	op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
		//	if ((op1 == 0) && (ofn1 == -5)) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
		//			for (jj = 0; jj < 3; jj++) {
		//				node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
		//				temp = (node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
		//				if (node_ii_loc != 0 && node_jj_loc != 0) {
		//					ma_full_T_q[temp] += 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel(node_ii, node_jj));
		//				}
		//				if (node_ii_loc != 0 && node_jj_loc == 0) {
		//					fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel(node_ii, node_jj));
		//					
		//				}
		//			}
		//			fbr_T_q1[(node_ii_loc - 1)] +=   T0 * 1.0 / 3 * element[el].face[nn].Area;
		//			
		//		}
		//	}
		//}


	}
	int matsize_T_q_temp = ma_full_T_q.size();


	//cout << "matsize_T_q_temp = " << matsize_T_q_temp << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if (op1 != 0) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					if (count == 10) {
						cout << "yess" << endl;
					}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}
					//k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2] + oppo_element[count_oppo_element-1].Kh[0] + oppo_element[count_oppo_element-1].Kh[1] + oppo_element[count_oppo_element-1].Kh[2])/6;
					k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2]) / 3;
					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid]+ node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}
					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver; 
						node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element-1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element-1].node[face_node[ofn1 - 1][jj] - 1].unknown_T;
							node_jj2_loc = oppo_element[count_oppo_element-1].node[face_node[ofn1 - 1][jj] - 1].unknown_q;

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element-1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element-1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element-1].Kh[0] + oppo_element[count_oppo_element-1].Kh[1] + oppo_element[count_oppo_element-1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] +=- 0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element-1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 *k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element-1].domain-1]+ node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}
				}
			}
		}
	}
	//fbt.resize(num_unKnown_Thermal);
	//fbt = fbr_T_q;
	//cout << "test77" << endl;
	int nnz_T_q = 0;
	int matsize_T_q = ma_full_T_q.size();

	int* iRow_T_q = new int[matsize_T_q];
	int* jCol_T_q = new int[matsize_T_q];

	double* Me_full_T_q = new double[matsize_T_q];

	for (auto full : ma_full_T_q) {
		Me_full_T_q[nnz_T_q] = full.second;
		iRow_T_q[nnz_T_q] = full.first / num_unKnown_Thermal + 1;
		jCol_T_q[nnz_T_q++] = full.first % num_unKnown_Thermal + 1;
	}
	ma_full_T_q.clear();
	int cc = 0; nnz_tot = 0;
	//cout << "matsize_T_q = " << matsize_T_q << endl;
	cout <<"nnz_T_q = "<< nnz_T_q << endl;
	int mm1_T_q = 0; int nn1_T_q = 0; int mm2_T_q = 0; int nn2_T_q = 0;

	m_T_q = new double[matsize_T_q];
	mrow_T_q = new int[matsize_T_q];
	mcol_T_q = new int[matsize_T_q];
	acoo_T_q = new double[matsize_T_q];
	rowind1 = new int[matsize_T_q];
	colind1 = new int[matsize_T_q];

	nnz_C.resize(num_domain, num_domain);
	nnz_C.setZero();

	num_nzero_Pmatrix_T_q = new int[num_domain];
	for (int n1 = myid; n1 < myid+1; n1++) {
		mm1_T_q = nn1_T_q + num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		nn2_T_q = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2_T_q = nn2_T_q + num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz_T_q; j++) {
				if (iRow_T_q[j] > nn1_T_q && iRow_T_q[j] <= mm1_T_q && jCol_T_q[j] > nn2_T_q && jCol_T_q[j] <= mm2_T_q) {
#						pragma omp critical
					{
						++nnz_dm;
						m_T_q[cc] = Me_full_T_q[j];
						mrow_T_q[cc] = iRow_T_q[j] - nn1_T_q;
						mcol_T_q[cc] = jCol_T_q[j] - nn2_T_q;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo_T_q[nnz_tot] = Me_full_T_q[j];
							rowind1[nnz_tot] = iRow_T_q[j] - nn1_T_q;
							colind1[nnz_tot] = jCol_T_q[j] - nn2_T_q;
							++nnz_tot;
						}
					}

				}
			}
			//cout << "myid is " << myid << endl;
			//cout << "nnz_dm = " << nnz_dm << endl;
			nn2_T_q = mm2_T_q;
			nnz_C(n1, n2) = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_T_q[n1] = nnz_dm;
		}
		nn1_T_q = mm1_T_q;
	}

	delete[] Me_full_T_q;
	delete[] iRow_T_q;
	delete[] jCol_T_q;

	p_pardiso_T_q = new double[nnz_tot]();
	ip_T_q = new int[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1] + 1]();
	jp_T_q = new int[nnz_tot]();

	r_dm_T_q = new int[num_domain];
	r_dm_T_q[0] = 0;
	nn3_dm_T_q = new int[num_domain];
	nn3_dm_T_q[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_T_q[i] = r_dm_T_q[i - 1] + Vertice_Boundary_T_q[i - 1];
		nn3_dm_T_q[i] = nn3_dm_T_q[i - 1];
		for (int j = 0; j < num_domain; j++) {
			nn3_dm_T_q[i] += nnz_C(i - 1, j);
		}
	}


	//cout << "num_nzero_Pmatrix_T_q["<<myid<<"] = " << num_nzero_Pmatrix_T_q[myid] << endl;




	return 0;
}
int Matrix_Generator_T_q_21(Element* element, Element* oppo_element, int myid) {
	ofstream ofs;
	cout.precision(16);
	long long temp;
	int ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int  nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;

	fbr_T_q.resize(num_unKnown_Thermal);
	fbr_T_q.setZero();


	fbr_T_q1 = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]; j++) {
		fbr_T_q1[j] = 0;
	}


	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 10; ii++) {
			//node_ii_glb = element[el].node[ii].ver;
			node_ii_loc = element[el].Eedge_T[ii];
			
			for (jj = 0; jj < 10; jj++) {
				//node_jj_glb = element[el].node[jj].ver; 
				node_jj_loc = element[el].Eedge_T[jj];
				
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					if (ii < 4 && jj < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						ma_full_T_q[temp] += (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
					}
					else if (jj < 4) {
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						int node_ii1 = edge_node_local[ii - 4][0]-1; int node_ii2 = edge_node_local[ii - 4][1]-1;
						//int node_jj1 = edge_node_local[jj - 4][0]-1; int node_jj2 = edge_node_local[jj - 4][1]-1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						//int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						
						ma_full_T_q[temp] += ((element[el].Kh[0] * element[el].b[node_ii1] * bj + element[el].Kh[1] * element[el].c[node_ii1] * cj + element[el].Kh[2] * element[el].d[node_ii1] * dj)\
							+(element[el].Kh[0] * element[el].b[node_ii2] * bj + element[el].Kh[1] * element[el].c[node_ii2] * cj + element[el].Kh[2] * element[el].d[node_ii2] * dj))/ (36.0 * element[el].Ve );
					}
					else if (ii < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

						//int node_ii1 = edge_node_local[ii - 4][0]-1; int node_ii2 = edge_node_local[ii - 4][1]-1;
						int node_jj1 = edge_node_local[jj - 4][0]-1; int node_jj2 = edge_node_local[jj - 4][1]-1;
						//int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						ma_full_T_q[temp] += ((element[el].Kh[0] * bi * element[el].b[node_jj1] + element[el].Kh[1] * ci * element[el].c[node_jj1] + element[el].Kh[2] * di * element[el].d[node_jj1])\
							+ (element[el].Kh[0] * bi * element[el].b[node_jj2] + element[el].Kh[1] * ci * element[el].c[node_jj2] + element[el].Kh[2] * di * element[el].d[node_jj2]))/ (36.0 * element[el].Ve);
					}
					else {
						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1= element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						ma_full_T_q[temp] += ((element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj1])*(1+ Flabel_T_q(node_ii_glb2, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb1))\
							+(element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb2))\
							+(element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb1)))/ (36.0 * element[el].Ve)*2.0/15.0*6.0;
					}

					
				}
				if ((node_ii_loc != 0) && (node_jj_loc == 0)) {
					cout << "erong" << endl;
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					if (ii < 4 && jj < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						fbr_T_q1[(node_ii_loc - 1)] -= Td * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
					}
					else if (jj < 4) {
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						//int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						//int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;


						fbr_T_q1[(node_ii_loc - 1)] -= Td * ((element[el].Kh[0] * element[el].b[node_ii1] * bj + element[el].Kh[1] * element[el].c[node_ii1] * cj + element[el].Kh[2] * element[el].d[node_ii1] * dj)\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * bj + element[el].Kh[1] * element[el].c[node_ii2] * cj + element[el].Kh[2] * element[el].d[node_ii2] * dj)) / (36.0 * element[el].Ve);
					}
					else if (ii < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						//int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						//int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						fbr_T_q1[(node_ii_loc - 1)] -= Td * ((element[el].Kh[0] * bi * element[el].b[node_jj1] + element[el].Kh[1] * ci * element[el].c[node_jj1] + element[el].Kh[2] * di * element[el].d[node_jj1])\
							+ (element[el].Kh[0] * bi * element[el].b[node_jj2] + element[el].Kh[1] * ci * element[el].c[node_jj2] + element[el].Kh[2] * di * element[el].d[node_jj2])) / (36.0 * element[el].Ve);
					}
					else {
						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						fbr_T_q1[(node_ii_loc - 1)] -= Td * ((element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb1))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb1))) / (36.0 * element[el].Ve) * 2.0 / 15.0*6.0;
					}
					//fbr_T_q1[(node_ii_loc - 1)] -= Td * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
			}
		}
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if ((op1 == 0) && (ofn1 == -7)) {
				for (ii = 0; ii < 3; ii++) {
					node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].Eedge_T[node_ii];
					for (jj = 0; jj < 3; jj++) {
						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
					}

					for (jj = 0; jj < 3; jj++) {
						int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0] - 1; int node_jj2 = edge_node_local[edge_jj - 1][1] - 1;

						//node_jj = face_node[nn][jj] - 1;
						node_jj_loc = element[el].Eedge_T[edge_jj+3];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1))*(1 + Flabel_T_q(node_ii, node_jj2));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
						}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 1.0 / 3 * element[el].face[nn].Area;
				}

				for (ii = 0; ii < 3; ii++) {

					int edge_ii = face_edge[nn][jj]; int node_ii1 = edge_node_local[edge_ii - 1][0] - 1; int node_ii2 = edge_node_local[edge_ii - 1][1] - 1;

					//node_ii = face_node[nn][ii] - 1;
					node_ii_loc = element[el].Eedge_T[edge_ii+3];
					for (jj = 0; jj < 3; jj++) {
						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj)) * (1 + Flabel_T_q(node_ii2, node_jj));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj)) * (1 + Flabel_T_q(node_ii2, node_jj));
						}
					}

					for (jj = 0; jj < 3; jj++) {
						int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0] - 1; int node_jj2 = edge_node_local[edge_jj - 1][1] - 1;

						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[edge_jj + 3];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 360 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj1)) * (1 + Flabel_T_q(node_ii1, node_jj2))*(1 + Flabel_T_q(node_ii2, node_jj1)) * (1 + Flabel_T_q(node_ii2, node_jj2));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 360 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj1)) * (1 + Flabel_T_q(node_ii1, node_jj2)) * (1 + Flabel_T_q(node_ii2, node_jj1)) * (1 + Flabel_T_q(node_ii2, node_jj2));
						}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 1.0 / 12 * element[el].face[nn].Area;

				}

			}
		}

	}
	int matsize_T_q_temp = ma_full_T_q.size();


	cout << "matsize_T_q_temp = " << matsize_T_q_temp << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if (op1 != 0) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					if (count == 10) {
						cout << "yess" << endl;
					}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}
					//k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2] + oppo_element[count_oppo_element-1].Kh[0] + oppo_element[count_oppo_element-1].Kh[1] + oppo_element[count_oppo_element-1].Kh[2])/6;
					k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2]) / 3;
					for (ii = 0; ii < 3; ii++) {


						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];
						//node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						//node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].Eedge_T[face_node[nn][jj] - 1];
							node_jj2_loc = element[el].Eedge_q[face_node[nn][jj] - 1];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}

						for (jj = 0; jj < 3; jj++) {

							int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = element[el].node[node_jj1 - 1].ver;
							int node_jj_glb2 = element[el].node[node_jj2 - 1].ver;
							node_jj1_loc = element[el].Eedge_T[edge_jj +3];
							node_jj2_loc = element[el].Eedge_q[edge_jj +3];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}

					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						//node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];

						int edge_ii = face_edge[nn][ii]; int node_ii1 = edge_node_local[edge_ii - 1][0]; int node_ii2 = edge_node_local[edge_ii - 1][1];

						int node_ii_glb1 = element[el].node[node_ii1 - 1].ver;
						int node_ii_glb2 = element[el].node[node_ii2 - 1].ver;
						node_ii1_loc = element[el].Eedge_T[edge_ii + 3];
						node_ii2_loc = element[el].Eedge_q[edge_ii + 3];

						//node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						//node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].Eedge_T[face_node[nn][jj] - 1];
							node_jj2_loc = element[el].Eedge_q[face_node[nn][jj] - 1];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb, node_ii_glb2));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}

						for (jj = 0; jj < 3; jj++) {

							int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = element[el].node[node_jj1 - 1].ver;
							int node_jj_glb2 = element[el].node[node_jj2 - 1].ver;
							node_jj1_loc = element[el].Eedge_T[edge_jj + 3];
							node_jj2_loc = element[el].Eedge_q[edge_jj + 3];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 360 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb2)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb2));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}

					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];


						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];



							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							//node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							//node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];


							int edge_jj = face_edge[ofn1 - 1][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = oppo_element[count_oppo_element - 1].node[node_jj1 - 1].ver;
							int node_jj_glb2 = oppo_element[count_oppo_element - 1].node[node_jj2 - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[edge_jj + 3];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[edge_jj + 3];


							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}




					}

					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						//node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];

						int edge_ii = face_edge[nn][ii]; int node_ii1 = edge_node_local[edge_ii - 1][0]; int node_ii2 = edge_node_local[edge_ii - 1][1];

						int node_ii_glb1 = element[el].node[node_ii1 - 1].ver;
						int node_ii_glb2 = element[el].node[node_ii2 - 1].ver;
						node_ii1_loc = element[el].Eedge_T[edge_ii + 3];
						node_ii2_loc = element[el].Eedge_q[edge_ii + 3];



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];



							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb, node_ii_glb2));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							//node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							//node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];


							int edge_jj = face_edge[ofn1 - 1][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = oppo_element[count_oppo_element - 1].node[node_jj1 - 1].ver;
							int node_jj_glb2 = oppo_element[count_oppo_element - 1].node[node_jj2 - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[edge_jj + 3];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[edge_jj + 3];


							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 360 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb1)) *(1 + Flabel_T_q(node_jj_glb1, node_ii_glb2)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb2));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}




					}


				}
			}
		}
	}
	//fbt.resize(num_unKnown_Thermal);
	//fbt = fbr_T_q;
	//cout << "test77" << endl;
	int nnz_T_q = 0;
	int matsize_T_q = ma_full_T_q.size();

	int* iRow_T_q = new int[matsize_T_q];
	int* jCol_T_q = new int[matsize_T_q];

	double* Me_full_T_q = new double[matsize_T_q];

	for (auto full : ma_full_T_q) {
		Me_full_T_q[nnz_T_q] = full.second;
		iRow_T_q[nnz_T_q] = full.first / num_unKnown_Thermal + 1;
		jCol_T_q[nnz_T_q++] = full.first % num_unKnown_Thermal + 1;
	}
	ma_full_T_q.clear();
	int cc = 0; nnz_tot = 0;
	cout << "matsize_T_q = " << matsize_T_q << endl;
	cout <<"nnz_T_q = "<< nnz_T_q << endl;
	int mm1_T_q = 0; int nn1_T_q = 0; int mm2_T_q = 0; int nn2_T_q = 0;

	m_T_q = new double[matsize_T_q];
	mrow_T_q = new int[matsize_T_q];
	mcol_T_q = new int[matsize_T_q];
	acoo_T_q = new double[matsize_T_q];
	rowind1 = new int[matsize_T_q];
	colind1 = new int[matsize_T_q];

	nnz_C.resize(num_domain, num_domain);
	nnz_C.setZero();

	num_nzero_Pmatrix_T_q = new int[num_domain];
	for (int n1 = myid; n1 < myid + 1; n1++) {
		mm1_T_q = nn1_T_q + num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		nn2_T_q = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2_T_q = nn2_T_q + num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz_T_q; j++) {
				if (iRow_T_q[j] > nn1_T_q && iRow_T_q[j] <= mm1_T_q && jCol_T_q[j] > nn2_T_q && jCol_T_q[j] <= mm2_T_q) {
#						pragma omp critical
					{
						++nnz_dm;
						m_T_q[cc] = Me_full_T_q[j];
						mrow_T_q[cc] = iRow_T_q[j] - nn1_T_q;
						mcol_T_q[cc] = jCol_T_q[j] - nn2_T_q;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo_T_q[nnz_tot] = Me_full_T_q[j];
							rowind1[nnz_tot] = iRow_T_q[j] - nn1_T_q;
							colind1[nnz_tot] = jCol_T_q[j] - nn2_T_q;
							++nnz_tot;
						}
					}

				}
			}
			//cout << "myid is " << myid << endl;
			//cout << "nnz_dm = " << nnz_dm << endl;
			nn2_T_q = mm2_T_q;
			nnz_C(n1, n2) = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_T_q[n1] = nnz_dm;
		}
		nn1_T_q = mm1_T_q;
	}

	delete[] Me_full_T_q;
	delete[] iRow_T_q;
	delete[] jCol_T_q;

	p_pardiso_T_q = new double[nnz_tot]();
	ip_T_q = new int[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1] + 1]();
	jp_T_q = new int[nnz_tot]();

	r_dm_T_q = new int[num_domain];
	r_dm_T_q[0] = 0;
	nn3_dm_T_q = new int[num_domain];
	nn3_dm_T_q[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_T_q[i] = r_dm_T_q[i - 1] + Vertice_Boundary_T_q[i - 1];
		nn3_dm_T_q[i] = nn3_dm_T_q[i - 1];
		for (int j = 0; j < num_domain; j++) {
			nn3_dm_T_q[i] += nnz_C(i - 1, j);
		}
	}


	cout << "num_nzero_Pmatrix_T_q["<<myid<<"] = " << num_nzero_Pmatrix_T_q[myid] << endl;




	return 0;
}

int Matrix_Generator_T_q_2(Element* element, Element* oppo_element, int myid) {
	ofstream ofs;
	cout.precision(16);
	long long temp;
	int ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int  nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;

	fbr_T_q.resize(num_unKnown_Thermal);
	fbr_T_q.setZero();


	fbr_T_q1 = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]; j++) {
		fbr_T_q1[j] = 0;
	}


	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 10; ii++) {
			//node_ii_glb = element[el].node[ii].ver;
			node_ii_loc = element[el].Eedge_T[ii];

			for (jj = 0; jj < 10; jj++) {
				//node_jj_glb = element[el].node[jj].ver; 
				node_jj_loc = element[el].Eedge_T[jj];

				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					if (ii < 4 && jj < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						ma_full_T_q[temp] += (-1.0 + 4.0 / 5.0 * (1 + Flabel_T_q(jj, ii))) * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
					}
					else if (jj < 4) {
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						//int node_jj1 = edge_node_local[jj - 4][0]-1; int node_jj2 = edge_node_local[jj - 4][1]-1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						//int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;


						ma_full_T_q[temp] += ((4.0 / 5.0 * (1 + Flabel_T_q(jj, node_ii2)) - 1.0) * (element[el].Kh[0] * element[el].b[node_ii1] * bj + element[el].Kh[1] * element[el].c[node_ii1] * cj + element[el].Kh[2] * element[el].d[node_ii1] * dj)\
							+ (4.0 / 5.0 * (1 + Flabel_T_q(jj, node_ii1)) - 1.0) * (element[el].Kh[0] * element[el].b[node_ii2] * bj + element[el].Kh[1] * element[el].c[node_ii2] * cj + element[el].Kh[2] * element[el].d[node_ii2] * dj)) / (36.0 * element[el].Ve);
					}
					else if (ii < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

						//int node_ii1 = edge_node_local[ii - 4][0]-1; int node_ii2 = edge_node_local[ii - 4][1]-1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						//int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						ma_full_T_q[temp] += ((4.0 / 5.0 * (1 + Flabel_T_q(ii, node_jj2)) - 1.0) * (element[el].Kh[0] * bi * element[el].b[node_jj1] + element[el].Kh[1] * ci * element[el].c[node_jj1] + element[el].Kh[2] * di * element[el].d[node_jj1])\
							+ (4.0 / 5.0 * (1 + Flabel_T_q(ii, node_jj1)) - 1.0) * (element[el].Kh[0] * bi * element[el].b[node_jj2] + element[el].Kh[1] * ci * element[el].c[node_jj2] + element[el].Kh[2] * di * element[el].d[node_jj2])) / (36.0 * element[el].Ve);
					}
					else {
						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						ma_full_T_q[temp] += ((element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb1))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb1))) / (36.0 * element[el].Ve) * 4.0 / 5.0;
					}


				}
				if ((node_ii_loc != 0) && (node_jj_loc == 0)) {
					//cout << "erong" << endl;
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					if (ii < 4 && jj < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						fbr_T_q1[(node_ii_loc - 1)] -= Td * (-1.0 + 4.0 / 5.0 * (1 + Flabel_T_q(jj, ii))) * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
					}
					else if (jj < 4) {
						aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];

						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						//int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						//int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;


						fbr_T_q1[(node_ii_loc - 1)] -= Td *  ((4.0 / 5.0 * (1 + Flabel_T_q(jj, node_ii2)) - 1.0) * (element[el].Kh[0] * element[el].b[node_ii1] * bj + element[el].Kh[1] * element[el].c[node_ii1] * cj + element[el].Kh[2] * element[el].d[node_ii1] * dj)\
							+ (4.0 / 5.0 * (1 + Flabel_T_q(jj, node_ii1)) - 1.0) * (element[el].Kh[0] * element[el].b[node_ii2] * bj + element[el].Kh[1] * element[el].c[node_ii2] * cj + element[el].Kh[2] * element[el].d[node_ii2] * dj)) / (36.0 * element[el].Ve);
					}
					else if (ii < 4) {
						ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						//int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						//int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						fbr_T_q1[(node_ii_loc - 1)] -= Td * ((4.0 / 5.0 * (1 + Flabel_T_q(ii, node_jj2)) - 1.0)*(element[el].Kh[0] * bi * element[el].b[node_jj1] + element[el].Kh[1] * ci * element[el].c[node_jj1] + element[el].Kh[2] * di * element[el].d[node_jj1])\
							+ (4.0 / 5.0 * (1 + Flabel_T_q(ii, node_jj1)) - 1.0) * (element[el].Kh[0] * bi * element[el].b[node_jj2] + element[el].Kh[1] * ci * element[el].c[node_jj2] + element[el].Kh[2] * di * element[el].d[node_jj2])) / (36.0 * element[el].Ve);
					}
					else {
						//int node_ii1 = edge_node_local[ii - 4][0]; int node_ii2 = edge_node_local[ii - 4][1];
						//int node_jj1 = edge_node_local[jj - 4][0]; int node_jj2 = edge_node_local[jj - 4][1];
						int node_ii1 = edge_node_local[ii - 4][0] - 1; int node_ii2 = edge_node_local[ii - 4][1] - 1;
						int node_jj1 = edge_node_local[jj - 4][0] - 1; int node_jj2 = edge_node_local[jj - 4][1] - 1;
						int node_ii_glb1 = element[el].node[node_ii1].ver; int node_ii_glb2 = element[el].node[node_ii2].ver;
						int node_jj_glb1 = element[el].node[node_jj1].ver; int node_jj_glb2 = element[el].node[node_jj2].ver;

						fbr_T_q1[(node_ii_loc - 1)] -= Td * ((element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii1] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii1] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii1] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb2, node_jj_glb1))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj1] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj1] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj1]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb2))\
							+ (element[el].Kh[0] * element[el].b[node_ii2] * element[el].b[node_jj2] + element[el].Kh[1] * element[el].c[node_ii2] * element[el].c[node_jj2] + element[el].Kh[2] * element[el].d[node_ii2] * element[el].d[node_jj2]) * (1 + Flabel_T_q(node_ii_glb1, node_jj_glb1))) / (36.0 * element[el].Ve) * 2.0 / 15.0 * 6.0;
					}
					//fbr_T_q1[(node_ii_loc - 1)] -= Td * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
			}
		}
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if ((op1 == 0) && (ofn1 == -7)) {        //convection  boundary condition
				for (ii = 0; ii < 3; ii++) {   //three nodes
					node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].Eedge_T[node_ii];
					for (jj = 0; jj < 3; jj++) {   //three nodes
						if (ii == jj) {
							node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
							temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
							if (node_ii_loc != 0 && node_jj_loc != 0) {
								ma_full_T_q[temp] += Hcvt * 1.0 / 30 * element[el].face[nn].Area ;
							}


						}
						else {
							node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
							temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
							if (node_ii_loc != 0 && node_jj_loc != 0) {
								ma_full_T_q[temp] += Hcvt *  element[el].face[nn].Area * (2.0 / 45.0 - 2.0 / 15.0 + 1.0 / 12.0);
							}
							if (node_ii_loc != 0 && node_jj_loc == 0) {
								fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt *  element[el].face[nn].Area * (2.0 / 45.0 - 2.0 / 15.0 + 1.0 / 12.0);
								//fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (2.0 / 45.0 - 2.0 / 15.0 + 1.0 / 12.0);
							}
						}

						//node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
						//temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						//if (node_ii_loc != 0 && node_jj_loc != 0) {
						//	ma_full_T_q[temp] += Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						//}
						//if (node_ii_loc != 0 && node_jj_loc == 0) {
						//	fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						//}
					}

					for (jj = 0; jj < 3; jj++) {    //three edges
						int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0] - 1; int node_jj2 = edge_node_local[edge_jj - 1][1] - 1;

						//node_jj = face_node[nn][jj] - 1;
						node_jj_loc = element[el].Eedge_T[edge_jj + 3];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (ii == node_jj1|| ii == node_jj2) {
							ma_full_T_q[temp] += 0.0;
							//if (node_ii_loc != 0 && node_jj_loc != 0) {
							//	ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
							//}
							//if (node_ii_loc != 0 && node_jj_loc == 0) {
							//	fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
							//}

						}
						else {
							if (node_ii_loc != 0 && node_jj_loc != 0) {
								ma_full_T_q[temp] += Hcvt * -4.0 / 180 * element[el].face[nn].Area ;
							}
							if (node_ii_loc != 0 && node_jj_loc == 0) {
								fbr_T_q1[(node_ii_loc - 1)] += Td * Hcvt * 4.0 / 180 * element[el].face[nn].Area ;
							}
						}


						
						//if (node_ii_loc != 0 && node_jj_loc != 0) {
						//	ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
						//}
						//if (node_ii_loc != 0 && node_jj_loc == 0) {
						//	fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
						//}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 0.0 * element[el].face[nn].Area;
				}

				for (ii = 0; ii < 3; ii++) { //three edges

					int edge_ii = face_edge[nn][jj]; int node_ii1 = edge_node_local[edge_ii - 1][0] - 1; int node_ii2 = edge_node_local[edge_ii - 1][1] - 1;

					//node_ii = face_node[nn][ii] - 1;
					node_ii_loc = element[el].Eedge_T[edge_ii + 3];
					for (jj = 0; jj < 3; jj++) { //three nodes
						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[node_jj];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;

						if (jj == node_ii1 || jj == node_ii2) {
							ma_full_T_q[temp] += 0.0;
							//if (node_ii_loc != 0 && node_jj_loc != 0) {
							//	ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
							//}
							//if (node_ii_loc != 0 && node_jj_loc == 0) {
							//	fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj1)) * (1 + Flabel_T_q(node_ii, node_jj2));
							//}

						}
						else {
							if (node_ii_loc != 0 && node_jj_loc != 0) {
								ma_full_T_q[temp] += Hcvt * -4.0 / 180 * element[el].face[nn].Area;
							}
							if (node_ii_loc != 0 && node_jj_loc == 0) {
								fbr_T_q1[(node_ii_loc - 1)] += Td * Hcvt * 4.0 / 180 * element[el].face[nn].Area;
							}
						}

						//if (node_ii_loc != 0 && node_jj_loc != 0) {
						//	ma_full_T_q[temp] += Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj)) * (1 + Flabel_T_q(node_ii2, node_jj));
						//}
						//if (node_ii_loc != 0 && node_jj_loc == 0) {
						//	fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 60 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj)) * (1 + Flabel_T_q(node_ii2, node_jj));
						//}
					}

					for (jj = 0; jj < 3; jj++) { //three edges
						int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0] - 1; int node_jj2 = edge_node_local[edge_jj - 1][1] - 1;

						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].Eedge_T[edge_jj + 3];
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 16.0 / 360 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj1)) * (1 + Flabel_T_q(node_ii1, node_jj2)) * (1 + Flabel_T_q(node_ii2, node_jj1)) * (1 + Flabel_T_q(node_ii2, node_jj2));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 16.0 / 360 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii1, node_jj1)) * (1 + Flabel_T_q(node_ii1, node_jj2)) * (1 + Flabel_T_q(node_ii2, node_jj1)) * (1 + Flabel_T_q(node_ii2, node_jj2));
						}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 1.0 / 3 * element[el].face[nn].Area;

				}

			}
		}

	}
	int matsize_T_q_temp = ma_full_T_q.size();


	cout << "matsize_T_q_temp = " << matsize_T_q_temp << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if (op1 != 0) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					if (count == 10) {
						cout << "yess" << endl;
					}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}
					k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2] + oppo_element[count_oppo_element-1].Kh[0] + oppo_element[count_oppo_element-1].Kh[1] + oppo_element[count_oppo_element-1].Kh[2])/6;
					//k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2]) / 3;


					//numerical flux
					for (ii = 0; ii < 3; ii++) {//three nodes


						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];
						//node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						//node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {//three nodes
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].Eedge_T[face_node[nn][jj] - 1];
							node_jj2_loc = element[el].Eedge_q[face_node[nn][jj] - 1];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;
							if (ii == jj) {
								factor = element[el].face[nn].Area * (4.0 / 15 - 2.0 / 5 + 1.0 / 6);
							}
							else {
								factor = element[el].face[nn].Area * (2.0 / 45 - 2.0 / 15 + 1.0 / 12);
							}


							//factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}

						for (jj = 0; jj < 3; jj++) { //three edges

							int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = element[el].node[node_jj1 - 1].ver;
							int node_jj_glb2 = element[el].node[node_jj2 - 1].ver;
							node_jj1_loc = element[el].Eedge_T[edge_jj + 3];
							node_jj2_loc = element[el].Eedge_q[edge_jj + 3];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;
							if (node_ii_glb == node_jj_glb1 || node_ii_glb == node_jj_glb2) {
								factor = 0;
							}
							else {
								factor = element[el].face[nn].Area * -4.0 / 180;
							}
							//factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}

					for (ii = 0; ii < 3; ii++) { //three edges
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						//node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];

						int edge_ii = face_edge[nn][ii]; int node_ii1 = edge_node_local[edge_ii - 1][0]; int node_ii2 = edge_node_local[edge_ii - 1][1];

						int node_ii_glb1 = element[el].node[node_ii1 - 1].ver;
						int node_ii_glb2 = element[el].node[node_ii2 - 1].ver;
						node_ii1_loc = element[el].Eedge_T[edge_ii + 3];
						node_ii2_loc = element[el].Eedge_q[edge_ii + 3];

						//node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						//node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) { //three nodes
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].Eedge_T[face_node[nn][jj] - 1];
							node_jj2_loc = element[el].Eedge_q[face_node[nn][jj] - 1];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;
							
							if (node_ii_glb1 == node_jj_glb || node_ii_glb2 == node_jj_glb) {
								factor = 0;
							}
							else {
								factor = element[el].face[nn].Area*-4.0 / 180;
							}

							//factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb, node_ii_glb2));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}

						for (jj = 0; jj < 3; jj++) { //three edges

							int edge_jj = face_edge[nn][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = element[el].node[node_jj1 - 1].ver;
							int node_jj_glb2 = element[el].node[node_jj2 - 1].ver;
							node_jj1_loc = element[el].Eedge_T[edge_jj + 3];
							node_jj2_loc = element[el].Eedge_q[edge_jj + 3];

							//node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							//node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area *16.0/ 360 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb2)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb2));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}




					//RTC
					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];


						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];

							if (node_ii_glb == node_jj_glb) {
								factor = element[el].face[nn].Area * (4.0 / 15 - 2.0 / 5 + 1.0 / 6);
							}
							else {
								factor = element[el].face[nn].Area * (2.0 / 45 - 2.0 / 15 + 1.0 / 12);
							}

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							//factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}

							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							//node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							//node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];


							int edge_jj = face_edge[ofn1 - 1][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = oppo_element[count_oppo_element - 1].node[node_jj1 - 1].ver;
							int node_jj_glb2 = oppo_element[count_oppo_element - 1].node[node_jj2 - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[edge_jj + 3];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[edge_jj + 3];

							if (node_ii_glb == node_jj_glb1 || node_ii_glb == node_jj_glb2) {
								factor = 0;
							}
							else {
								factor = element[el].face[nn].Area*-4.0 / 180;
							}


							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							//factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}




					}

					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = element[el].Eedge_T[face_node[nn][ii] - 1];
						//node_ii2_loc = element[el].Eedge_q[face_node[nn][ii] - 1];

						int edge_ii = face_edge[nn][ii]; int node_ii1 = edge_node_local[edge_ii - 1][0]; int node_ii2 = edge_node_local[edge_ii - 1][1];

						int node_ii_glb1 = element[el].node[node_ii1 - 1].ver;
						int node_ii_glb2 = element[el].node[node_ii2 - 1].ver;
						node_ii1_loc = element[el].Eedge_T[edge_ii + 3];
						node_ii2_loc = element[el].Eedge_q[edge_ii + 3];



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];

							if (node_ii_glb1 == node_jj_glb || node_ii_glb2 == node_jj_glb) {
								factor = 0;
							}
							else {
								factor = element[el].face[nn].Area *-4.0/ 180;
							}

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							//factor = element[el].face[nn].Area / 60 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb, node_ii_glb2));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}



						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							//node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[face_node[ofn1 - 1][jj] - 1];
							//node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[face_node[ofn1 - 1][jj] - 1];


							int edge_jj = face_edge[ofn1 - 1][jj]; int node_jj1 = edge_node_local[edge_jj - 1][0]; int node_jj2 = edge_node_local[edge_jj - 1][1];

							int node_jj_glb1 = oppo_element[count_oppo_element - 1].node[node_jj1 - 1].ver;
							int node_jj_glb2 = oppo_element[count_oppo_element - 1].node[node_jj2 - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].Eedge_T[edge_jj + 3];
							node_jj2_loc = oppo_element[count_oppo_element - 1].Eedge_q[edge_jj + 3];


							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area*16.0 / 360 * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb1)) * (1 + Flabel_T_q(node_jj_glb1, node_ii_glb2)) * (1 + Flabel_T_q(node_jj_glb2, node_ii_glb2));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq

							//double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];
							//temp_Kh = temp_Kh / 3;
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}




					}


				}
			}
		}
	}
	//fbt.resize(num_unKnown_Thermal);
	//fbt = fbr_T_q;
	//cout << "test77" << endl;
	int nnz_T_q = 0;
	int matsize_T_q = ma_full_T_q.size();

	int* iRow_T_q = new int[matsize_T_q];
	int* jCol_T_q = new int[matsize_T_q];

	double* Me_full_T_q = new double[matsize_T_q];

	for (auto full : ma_full_T_q) {
		Me_full_T_q[nnz_T_q] = full.second;
		iRow_T_q[nnz_T_q] = full.first / num_unKnown_Thermal + 1;
		jCol_T_q[nnz_T_q++] = full.first % num_unKnown_Thermal + 1;
	}
	ma_full_T_q.clear();
	int cc = 0; nnz_tot = 0;
	cout << "matsize_T_q = " << matsize_T_q << endl;
	cout << "nnz_T_q = " << nnz_T_q << endl;
	int mm1_T_q = 0; int nn1_T_q = 0; int mm2_T_q = 0; int nn2_T_q = 0;

	m_T_q = new double[matsize_T_q];
	mrow_T_q = new int[matsize_T_q];
	mcol_T_q = new int[matsize_T_q];
	acoo_T_q = new double[matsize_T_q];
	rowind1 = new int[matsize_T_q];
	colind1 = new int[matsize_T_q];

	nnz_C.resize(num_domain, num_domain);
	nnz_C.setZero();

	num_nzero_Pmatrix_T_q = new int[num_domain];
	for (int n1 = myid; n1 < myid + 1; n1++) {
		mm1_T_q = nn1_T_q + num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		nn2_T_q = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2_T_q = nn2_T_q + num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz_T_q; j++) {
				if (iRow_T_q[j] > nn1_T_q && iRow_T_q[j] <= mm1_T_q && jCol_T_q[j] > nn2_T_q && jCol_T_q[j] <= mm2_T_q) {
#						pragma omp critical
					{
						++nnz_dm;
						m_T_q[cc] = Me_full_T_q[j];
						mrow_T_q[cc] = iRow_T_q[j] - nn1_T_q;
						mcol_T_q[cc] = jCol_T_q[j] - nn2_T_q;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo_T_q[nnz_tot] = Me_full_T_q[j];
							rowind1[nnz_tot] = iRow_T_q[j] - nn1_T_q;
							colind1[nnz_tot] = jCol_T_q[j] - nn2_T_q;
							++nnz_tot;
						}
					}

				}
			}
			//cout << "myid is " << myid << endl;
			//cout << "nnz_dm = " << nnz_dm << endl;
			nn2_T_q = mm2_T_q;
			nnz_C(n1, n2) = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_T_q[n1] = nnz_dm;
		}
		nn1_T_q = mm1_T_q;
	}

	delete[] Me_full_T_q;
	delete[] iRow_T_q;
	delete[] jCol_T_q;

	p_pardiso_T_q = new double[nnz_tot]();
	ip_T_q = new int[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1] + 1]();
	jp_T_q = new int[nnz_tot]();

	r_dm_T_q = new int[num_domain];
	r_dm_T_q[0] = 0;
	nn3_dm_T_q = new int[num_domain];
	nn3_dm_T_q[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_T_q[i] = r_dm_T_q[i - 1] + Vertice_Boundary_T_q[i - 1];
		nn3_dm_T_q[i] = nn3_dm_T_q[i - 1];
		for (int j = 0; j < num_domain; j++) {
			nn3_dm_T_q[i] += nnz_C(i - 1, j);
		}
	}


	cout << "num_nzero_Pmatrix_T_q[" << myid << "] = " << num_nzero_Pmatrix_T_q[myid] << endl;




	return 0;
}

int Matrix_Generator_T_q_E(Element* element, Element* oppo_element, int myid) {
	ofstream ofs;
	cout.precision(16);
	long long temp;
	int ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int  nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;
	fbr_T_q.resize(num_unKnown_Thermal);
	fbr_T_q.setZero();
	fbr_T_q1 = new double[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]];

	for (int j = 0; j < num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1]; j++) {
		fbr_T_q1[j] = 0;
	}
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii_glb = element[el].node[ii].ver; node_ii_loc = element[el].node[ii].unknown_T;
			ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];
			for (jj = 0; jj < 4; jj++) {
				node_jj_glb = element[el].node[jj].ver; node_jj_loc = element[el].node[jj].unknown_T;
				aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
					ma_full_T_q[temp] += (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
				if ((node_ii_loc != 0) && (node_jj_loc == 0)) {
					fbr_T_q1[(node_ii_loc - 1)] -= Td * (element[el].Kh[0] * bi * bj + element[el].Kh[1] * ci * cj + element[el].Kh[2] * di * dj) / (36.0 * element[el].Ve);
				}
			}
		}
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if ((op1 == 0) && (ofn1 == -7)) {
				for (ii = 0; ii < 3; ii++) {
					node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
					for (jj = 0; jj < 3; jj++) {
						node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
						temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
						if (node_ii_loc != 0 && node_jj_loc != 0) {
							ma_full_T_q[temp] += Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
						if (node_ii_loc != 0 && node_jj_loc == 0) {
							fbr_T_q1[(node_ii_loc - 1)] -= Td * Hcvt * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
						}
					}
					fbr_T_q1[(node_ii_loc - 1)] += Hcvt * Ta * 1.0 / 3 * element[el].face[nn].Area;

				}
			}
		}

		//for (nn = 0; nn < 4; nn++) {
		//	op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
		//	if ((op1 == 0) && (ofn1 == -7 )) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = element[el].node[node_ii].unknown_T;
		//			for (jj = 0; jj < 3; jj++) {
		//				node_jj = face_node[nn][jj] - 1; node_jj_loc = element[el].node[node_jj].unknown_T;
		//				temp = (long long)(node_ii_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj_loc - 1;
		//				if (node_ii_loc != 0 && node_jj_loc != 0) {
		//					ma_full_T_q[temp] += 0 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//				if (node_ii_loc != 0 && node_jj_loc == 0) {
		//					fbr_T_q1[(node_ii_loc - 1)] -= Td * 0 * 1.0 / 12 * element[el].face[nn].Area * (1 + Flabel_T_q(node_ii, node_jj));
		//				}
		//			}
		//			fbr_T_q1[(node_ii_loc - 1)] += 0 * Ta * 1.0 / 3 * element[el].face[nn].Area;
		//		}
		//	}
		//}

	}
	int matsize_T_q_temp = ma_full_T_q.size();


	cout << "matsize_T_q_temp = " << matsize_T_q_temp << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (nn = 0; nn < 4; nn++) {
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];
			if (op1 != 0) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					if (count == 10) {
						cout << "yess" << endl;
					}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}
					//k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2] + oppo_element[count_oppo_element-1].Kh[0] + oppo_element[count_oppo_element-1].Kh[1] + oppo_element[count_oppo_element-1].Kh[2])/6;
					k = (element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2]) / 3;
					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_T;
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_q;

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;     //FTT
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] -= 4.0 * factor;
								ma_full_T_q[temp] -= 4.0 * factor;
							}


							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) += Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] += Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;     //FTq


							double temp_Kh = element[el].Kh[0] + element[el].Kh[1] + element[el].Kh[2];

							temp_Kh = temp_Kh / 3;
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 0.5 * 1 * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj1_loc - 1;  //FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								ma_full_T_q[temp] += k * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[myid] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}
					for (ii = 0; ii < 3; ii++) {
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_T;
						node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_q;
						for (jj = 0; jj < 3; jj++) {
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_T;
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_q;

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;   //FTT

							factor = element[el].face[nn].Area / 12 * (1 + Flabel_T_q(node_jj_glb, node_ii_glb));
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								//ma_full_T_q[temp] += 4.0 * factor;
								ma_full_T_q[temp] += 4.0 * factor;

							}
							if (node_ii1_loc != 0 && node_jj1_loc == 0) {
								//fbr_T_q(node_ii1_loc - 1) -= Td * 4.0 * factor;
								fbr_T_q1[(node_ii1_loc - 1)] -= Td * 4.0 * factor;
							}
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//FTq
							double temp_Kh = oppo_element[count_oppo_element - 1].Kh[0] + oppo_element[count_oppo_element - 1].Kh[1] + oppo_element[count_oppo_element - 1].Kh[2];
							temp_Kh = temp_Kh / 3;

							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								ma_full_T_q[temp] += -0.5 * 1 * factor;

							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj1_loc - 1;//FqT
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								//ma_full[temp] += -1.0 * k * factor;

								ma_full_T_q[temp] += -1.0 * k * factor;
							}
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_Thermal + Accumulated_unknowns_T_q[oppo_element[count_oppo_element - 1].domain - 1] + node_jj2_loc - 1;//Fqq
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								//mat_full1[temp] += 0;
								//mat_full2[temp] += temp_Kh * factor;
								ma_full_T_q[temp] += 1 * factor;
							}
						}
					}
				}
			}
		}
	}
	int nnz_T_q = 0;
	int matsize_T_q = ma_full_T_q.size();

	int* iRow_T_q = new int[matsize_T_q];
	int* jCol_T_q = new int[matsize_T_q];

	double* Me_full_T_q = new double[matsize_T_q];

	for (auto full : ma_full_T_q) {
		Me_full_T_q[nnz_T_q] = full.second;
		iRow_T_q[nnz_T_q] = full.first / num_unKnown_Thermal + 1;
		jCol_T_q[nnz_T_q++] = full.first % num_unKnown_Thermal + 1;
	}

	int cc = 0; nnz_tot = 0;
	cout << "matsize_T_q = " << matsize_T_q << endl;
	//cout <<"nnz_T_q = "<< nnz_T_q << endl;
	int mm1_T_q = 0; int nn1_T_q = 0; int mm2_T_q = 0; int nn2_T_q = 0;

	m_T_q = new double[matsize_T_q];
	mrow_T_q = new int[matsize_T_q];
	mcol_T_q = new int[matsize_T_q];
	acoo_T_q = new double[matsize_T_q];
	rowind1 = new int[matsize_T_q];
	colind1 = new int[matsize_T_q];

	nnz_C.resize(num_domain, num_domain);
	nnz_C.setZero();

	num_nzero_Pmatrix_T_q = new int[num_domain];
	for (int n1 = myid; n1 < myid + 1; n1++) {
		mm1_T_q = nn1_T_q + num_unknown_subdomain_T_q[n1][0] + num_unknown_subdomain_T_q[n1][1];
		nn2_T_q = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2_T_q = nn2_T_q + num_unknown_subdomain_T_q[n2][0] + num_unknown_subdomain_T_q[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz_T_q; j++) {
				if (iRow_T_q[j] > nn1_T_q && iRow_T_q[j] <= mm1_T_q && jCol_T_q[j] > nn2_T_q && jCol_T_q[j] <= mm2_T_q) {
#						pragma omp critical
					{
						++nnz_dm;
						m_T_q[cc] = Me_full_T_q[j];
						mrow_T_q[cc] = iRow_T_q[j] - nn1_T_q;
						mcol_T_q[cc] = jCol_T_q[j] - nn2_T_q;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo_T_q[nnz_tot] = Me_full_T_q[j];
							rowind1[nnz_tot] = iRow_T_q[j] - nn1_T_q;
							colind1[nnz_tot] = jCol_T_q[j] - nn2_T_q;
							++nnz_tot;
						}
					}

				}
			}
			//cout << "myid is " << myid << endl;
			//cout << "nnz_dm = " << nnz_dm << endl;
			nn2_T_q = mm2_T_q;
			nnz_C(n1, n2) = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_T_q[n1] = nnz_dm;
		}
		nn1_T_q = mm1_T_q;
	}

	delete[] Me_full_T_q;
	delete[] iRow_T_q;
	delete[] jCol_T_q;

	p_pardiso_T_q = new double[nnz_tot]();
	ip_T_q = new int[num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1] + 1]();
	jp_T_q = new int[nnz_tot]();

	r_dm_T_q = new int[num_domain];
	r_dm_T_q[0] = 0;
	nn3_dm_T_q = new int[num_domain];
	nn3_dm_T_q[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_T_q[i] = r_dm_T_q[i - 1] + Vertice_Boundary_T_q[i - 1];
		nn3_dm_T_q[i] = nn3_dm_T_q[i - 1];
		for (int j = 0; j < num_domain; j++) {
			nn3_dm_T_q[i] += nnz_C(i - 1, j);
		}
	}


	//cout << "num_nzero_Pmatrix_T_q["<<myid<<"] = " << num_nzero_Pmatrix_T_q[myid] << endl;




	return 0;
}



double Flabel_T_q(int vnode1, int vnode2) {
	double Flable0;
	if (vnode1 == vnode2) {
		Flable0 = 1.0;
	}
	else {
		Flable0 = 0.0;
	}
	return Flable0;
}
