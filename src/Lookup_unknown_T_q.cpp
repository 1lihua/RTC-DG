
#include "Lookup_unknown_T_q.h"
#include<fstream>
using namespace Eigen;

using namespace std;

Eigen::MatrixXi num_unKnown_subdomain;

Eigen::VectorXi Node_Boundary_Thermal;


int Read_BounaryInfo2(vector<int>& DIE1_info, vector<int>& DIE2_info, vector<int>& DIE3_info,
	vector<int>& DIE4_info,
	vector<int>& DIE5_info,
	vector<int>& DIE6_info,
	vector<int>& DIE7_info,
	vector<int>& DIE8_info,
	vector<int>& DIE9_info,
	vector<int>& DIE10_info,
	vector<int>& DIE11_info,
	vector<int>& DIE12_info,
	vector<int>& DIE13_info,
	vector<int>& DIE14_info,
	vector<int>& DIE15_info,
	vector<int>& DIE16_info);
void s2i2(const string& info, const string& sta,
	vector<int>& DIE1_info,
	vector<int>& DIE2_info,
	vector<int>& DIE3_info,
	vector<int>& DIE4_info,
	vector<int>& DIE5_info,
	vector<int>& DIE6_info,
	vector<int>& DIE7_info,
	vector<int>& DIE8_info,
	vector<int>& DIE9_info,
	vector<int>& DIE10_info,
	vector<int>& DIE11_info,
	vector<int>& DIE12_info,
	vector<int>& DIE13_info,
	vector<int>& DIE14_info,
	vector<int>& DIE15_info,
	vector<int>& DIE16_info);





struct Sum1 {
	size_t idx;
	size_t sum;
	//按math从大到小排序
	inline bool operator < (const Sum1& x) const {
		return sum < x.sum;
	}

};

struct EdgePro1 {
	size_t idx;
	size_t  pro;
	//按math从大到小排序
	inline bool operator < (const EdgePro1& x) const {
		return pro < x.pro;
	}

};

#define Material_1 2   
#define Material_2 100   


int LookUp_glb_Edge_Index(Element* el, vector<vector<int>>& Edge_GBNO) {

	int num_edge_T_q = 0;
	vector<Sum1> opp_Sum((size_t)num_element_subdomain * 6);
#pragma omp parallel for
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (int ii = 0; ii < 6; ++ii) {
			int f_n1 = edge_node_local[ii][0] - 1, f_n2 = edge_node_local[ii][1] - 1;
			opp_Sum[(size_t)i * 6 + ii].idx = (size_t)i * 6 + ii;
			opp_Sum[(size_t)i * 6 + ii].sum = (size_t)el[i].node[f_n1].ver + (size_t)el[i].node[f_n2].ver;
		}
	}
	sort(opp_Sum.begin(), opp_Sum.end());
	for (auto beg = opp_Sum.begin(); beg != opp_Sum.end(); ) {
		size_t sum1 = (*beg).sum;
		vector<EdgePro1> v2;
		while (sum1 == (*beg).sum) {
			int ii = ((*beg).idx) % 6; int ith = ((*beg).idx) / 6; EdgePro1 temp1;
			int f_n1 = edge_node_local[ii][0] - 1, f_n2 = edge_node_local[ii][1] - 1;
			long long pro_1 = (size_t)el[ith].node[f_n1].ver * el[ith].node[f_n2].ver;
			temp1.idx = (*beg).idx; temp1.pro = pro_1;
			v2.push_back(temp1);
			beg++;
		}
		int v2Size = v2.size();
		//#pragma omp parallel for
		for (int i = 0; i < v2Size; ++i) {
			int ith = v2[i].idx / 6; int nn = v2[i].idx % 6;
			if ((Edge_GBNO[ith][nn] == 0)) { Edge_GBNO[ith][nn] = ++num_edge_T_q; }
#pragma omp parallel for
			for (int j = i + 1; j < v2Size; ++j) {
				if (v2[i].pro == v2[j].pro) {
					int op1 = v2[j].idx / 6; int ofn1 = v2[j].idx % 6;
					Edge_GBNO[op1][ofn1] = Edge_GBNO[ith][nn];
				}
			}
		}
	}
	//num_edge_T_q++;

	for (int i = 0; i < num_element_subdomain; ++i) {
		for (int ii = 0; ii < 6; ++ii) {
			Edge_GBNO[i][ii] -= 1;
		}
	}

	return num_edge_T_q;

}


//int unknownIndex_Thermal_2(Element* element,int myid) {
//
//	cout << "test0" << endl;
//
//	vector<int> vecInit(6, 0);
//	//for (int i = 0; i < 6; ++i) {
//	//	vecInit[i] = 1;
//	//}
//	vector<vector<int>>Edge_GBNO(num_element_subdomain, vecInit);
//
//	int num_edge_T_q=LookUp_glb_Edge_Index(element, Edge_GBNO);
//
//
//	cout << "num_edge_T_q = "<< num_edge_T_q << endl;
//	vector<int>Edge_BC(num_edge_T_q, 0);
//
//	vector<int> Edge_Loc(num_edge_T_q,0);
//	vector<int> Edge_Loc_q(num_edge_T_q, 0);
//
//	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
//		for (int nn = 0; nn < 4; nn++) {
//			int op1 = element[el0].face[nn].opp[1]; int op0 = element[el0].face[nn].opp[0];
//			if (op1 == -100) {  //
//				for (int ii = 0; ii < 3; ii++) {
//					int edgeii = Edge_GBNO[el0][face_edge[nn][ii]-1];
//					Edge_BC[edgeii] = -1;// Dirichlet boundary conditions  of each edge
//				}
//			}
//		}
//	}
//	vector<int> boundary_element;
//	for (int i = 0; i < num_element_subdomain; i++) {
//		for (int nn = 0; nn < 4; nn++) {
//			int op0 = element[i].face[nn].opp[0]; int ofn = element[i].face[nn].opp[1];
//			if (op0 != 0) {
//				if (element[i].face[nn].whether_boundary) {
//					boundary_element.push_back(i);
//					for (int ii = 0; ii < 3; ++ii) {
//						int edgeii = Edge_GBNO[i][face_edge[nn][ii]-1];
//						if (Edge_BC[edgeii] != -1)
//						{
//							Edge_Loc[edgeii] = -1;//not Dirichlet boundary conditions and on the boundary
//						}
//					}
//				}
//			}
//		}
//	}
//	int edge_tot = 0;
//	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
//		for (int ii = 0; ii < 6; ii++) {
//			int edgeii = Edge_GBNO[el0][ii];
//			if ((Edge_Loc[edgeii] == 0) && (Edge_BC[edgeii] != -1)) {//not Dirichlet boundary conditions and Not on the boundary
//				edge_tot += 1;
//				Edge_Loc[edgeii] = edge_tot;
//			}
//		}
//	}
//
//
//	cout << "test2" << endl;
//
//	Eigen::MatrixXi Node_RTC_Thermal;
//	Node_RTC_Thermal.resize(num_node + 1, num_domain);
//
//	int* Node_RTC_T_q = new int[num_node];
//
//	for (int i = 0; i < num_node; i++) {
//		Node_RTC_T_q[i] = 0;
//	}
//
//
//	int el, nn, op, ofn, ii, jj, node1, node2, node3, ndm, node_tot, node_ii, ncnt, num_node_Dirichlet, nUnknown_tot, nodeii;
//	int ofn1, op1, op2, op3, op4, opd1, opd2, opd3, opd4;
//	VectorXi Node_Dirichlet, num_node_SubDomain;
//	MatrixXi Node_Loc;
//	num_unKnown_eb = 0;
//	int num_unKnown_Thermal = 0;
//	Node_Dirichlet.resize(num_node);
//	Node_Dirichlet.setZero();
//	int kdirichlet = 0;
//
//	for (el = 0; el < num_element_subdomain; el++) {
//		for (nn = 0; nn < 4; nn++) {
//			op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
//			if ((op == 0) && (ofn == -100)) {
//				node1 = element[el].node[face_node[nn][0] - 1].ver - 1;
//				node2 = element[el].node[face_node[nn][1] - 1].ver - 1;
//				node3 = element[el].node[face_node[nn][2] - 1].ver - 1;
//				kdirichlet++;
//				Node_Dirichlet(node1) = -1;
//				Node_Dirichlet(node2) = -1;
//				Node_Dirichlet(node3) = -1;
//			}
//		}
//	}
//
//	//cout << "kdirichlet = " << kdirichlet << endl;
//	//cout << "kconvection = " << kconvection << endl;
//	num_node_SubDomain.resize(num_domain);
//	Node_Loc.resize(num_node, num_domain);
//	Node_Loc.setZero();
//
//	for (ndm = myid; ndm < myid+1; ndm++) {
//		node_tot = edge_tot;
//		int cnt_boundary = 0;
//		for (el = 0; el < num_element_subdomain; el++) {
//			if (element[el].domain - 1 == ndm) {
//				for (nn = 0; nn < 4; nn++) {
//					op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
//					if (op != 0) {
//						if(element[el].face[nn].whether_boundary){
//							for (ii = 0; ii < 3; ii++) {
//								node_ii = element[el].node[face_node[nn][ii] - 1].ver - 1;
//								if (Node_Dirichlet(node_ii) != -1) {
//									Node_Loc(node_ii, ndm) = -1;
//								}
//							}
//						}
//					}
//				}
//			}
//			else {
//				cout << "wrong domain" << endl;
//			}
//
//
//		}
//		for (int el = 0; el < num_element_subdomain; el++) {
//			if (element[el].domain - 1 == ndm) {
//				for (ii = 0; ii < 4; ii++) {
//					node_ii = element[el].node[ii].ver - 1;
//					if (Node_Dirichlet(node_ii) != -1) {
//						if (Node_Loc(node_ii, ndm) == 0) {
//							Node_Loc(node_ii, ndm) = ++node_tot;//vertice not on the boundary
//						}
//					}
//				}
//			}
//			else {
//				cout << "wrong domain" << endl;
//			}
//
//
//		}
//		for (int el = 0; el < num_element_subdomain; el++) {
//			if (element[el].domain - 1 == ndm) {
//				for (ii = 0; ii < 4; ii++) {
//					node_ii = element[el].node[ii].ver - 1;
//					if (Node_Dirichlet(node_ii) != -1) {
//						if (Node_Loc(node_ii, ndm) == -1) {
//							Node_Loc(node_ii, ndm) = ++node_tot;
//							cnt_boundary++;//vertice on the boundary
//						}
//					}
//				}
//			}
//			else {
//				cout << "wrong domain" << endl;
//			}
//		}
//		
//		for (int e = 0; e < boundary_element.size(); e++) {
//			int el0 = boundary_element[e];
//			for (int ii = 0; ii < 6; ii++) {
//				int edgeii = Edge_GBNO[el0][ii];
//				if ((Edge_Loc[edgeii] == -1) && (Edge_BC[edgeii] != -1)) {//not Dirichlet and on the boundary
//					node_tot += 1;
//					Edge_Loc[edgeii] = node_tot;
//					++cnt_boundary;
//
//				}
//			}
//
//		}
//
//
//
//		Vertice_Boundary_T_q[ndm] = cnt_boundary;
//		num_unknown_subdomain_T_q[ndm][0] = node_tot;
//		num_unKnown_eb += cnt_boundary;
//		cout << "nNode_IN=" << ndm << '\t' << cnt_boundary << '\t' << node_tot << endl;
//	}
//	cout << "test3" << endl;
//	for (el = 0; el < num_element_subdomain; el++) {
//		ndm = element[el].domain;
//		for (ii = 0; ii < 4; ii++) {
//			node_ii = element[el].node[ii].ver;
//			element[el].node[ii].unknown_T = Node_Loc(node_ii - 1, ndm - 1);
//		}
//		for (ii = 0; ii < 4; ii++) {
//			node_ii = element[el].node[ii].ver;
//			element[el].Eedge_T[ii] = Node_Loc(node_ii - 1, ndm - 1);
//		}
//
//
//	}
//
//	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
//		for (int ii = 0; ii < 6; ii++) {
//			int edgeii = Edge_GBNO[el0][ii];
//			if (Edge_BC[edgeii] != -1) {
//				element[el0].Eedge_T[ii+4] = Edge_Loc[edgeii];
//				
//			}
//		}
//	}
//
//
//	Node_RTC_Thermal.setZero();
//
//	cout << "test4" << endl;
//
//	for (ndm = myid; ndm < myid+1; ndm++) {
//		ncnt = 0;
//		for (el = 0; el < num_element_subdomain; el++) {
//			if (element[el].domain - 1 == ndm) {
//				for (nn = 0; nn < 4; nn++) {
//					op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
//					if (op != 0) {
//						if (element[el].face[nn].whether_boundary) {
//							for (ii = 0; ii < 3; ii++) {
//								nodeii = element[el].node[face_node[nn][ii] - 1].ver - 1;
//								if (Node_RTC_T_q[(nodeii)] == 0) {
//									ncnt++;
//									Node_RTC_T_q[(nodeii)] = ncnt;
//								}
//								if (Node_Dirichlet(nodeii) == -1) {      //error
//									cout << "ff" << endl;
//									cout << "op = " << op << endl;
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//
//
//		for (el = 0; el < num_element_subdomain; el++) {
//			for (nn = 0; nn < 4; nn++) {
//				if (element[el].face[nn].whether_boundary) {
//					for (ii = 0; ii < 3; ii++) {
//						int edge_num = face_edge[nn][ii];
//						int edgeii = Edge_GBNO[el][edge_num-1];
//						if (Edge_Loc_q[edgeii] == 0) {
//							Edge_Loc_q[edgeii] = ++ncnt;
//						}
//					}
//					
//
//					
//				}
//
//			}
//		}
//
//
//
//
//
//		num_unknown_subdomain_T_q[ndm][1] = ncnt;
//		Vertice_Boundary_T_q[ndm] += ncnt;
//		num_unKnown_eb += ncnt;
//		cout << "nNode_Bc=" << ndm << '\t' << ncnt << endl;
//	}
//
//
//	for (el = 0; el < num_element_subdomain; el++) {
//		ndm = element[el].domain;
//		for (ii = 0; ii < 4; ii++) {
//			node_ii = element[el].node[ii].ver;
//			element[el].node[ii].unknown_q = Node_RTC_T_q[(node_ii - 1)];
//		}
//		for (ii = 0; ii < 4; ii++) {
//			node_ii = element[el].node[ii].ver;
//			element[el].Eedge_q[ii] = Node_RTC_T_q[(node_ii - 1)];
//		}
//
//	}
//	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
//		for (int ii = 0; ii < 6; ii++) {
//			int edgeii = Edge_GBNO[el0][ii];
//			element[el0].Eedge_q[ii + 4] = Edge_Loc_q[edgeii];
//
//			
//		}
//	}
//
//
//	nUnknown_tot = 0;
//	int unknown_q = 0;
//	//int unknown_T = 0;
//
//	ndm = myid;
//	nUnknown_tot += num_unknown_subdomain_T_q[ndm][0];
//	for (el = 0; el < num_element_subdomain; el++) {
//		for (ii = 0; ii < 10; ii++) {
//			if (element[el].Eedge_q[ii]!= 0) {
//				element[el].Eedge_q[ii] += nUnknown_tot;
//			}
//
//		}
//	}
//	nUnknown_tot += num_unknown_subdomain_T_q[ndm][1];
//	unknown_q += num_unknown_subdomain_T_q[ndm][1];
//	num_unKnown_Thermal += num_unknown_subdomain_T_q[ndm][0] + num_unknown_subdomain_T_q[ndm][1];
//
//	//for (ndm = myid; ndm < myid+1; ndm++) {
//	//	for (el = 0; el < num_element_subdomain; el++) {
//	//		if (element[el].domain - 1 == ndm) {
//	//			for (ii = 0; ii < 4; ii++) {
//	//				if (element[el].node[ii].unknown_T != 0) {
//	//					element[el].node[ii].unknown_T += nUnknown_tot;
//	//				}
//	//			}
//	//		}
//	//	}
//	//	nUnknown_tot += num_unknown_subdomain_T_q[ndm][0];
//	//	for (el = 0; el < num_element_subdomain; el++) {
//	//		for (ii = 0; ii < 4; ii++) {
//	//			if (element[el].node[ii].unknown_q != 0) {
//	//				element[el].node[ii].unknown_q += nUnknown_tot;
//	//			}
//	//		
//	//		}
//	//	}
//	//	nUnknown_tot += num_unknown_subdomain_T_q[ndm][1];
//	//	unknown_q += num_unknown_subdomain_T_q[ndm][1];
//	//	num_unKnown_Thermal += num_unknown_subdomain_T_q[ndm][0] + num_unknown_subdomain_T_q[ndm][1];
//
//	//}
//
//	//cout << "num_unKnown=" << num_unKnown_Thermal << '\t' << nUnknown_tot << endl;
//	cout << "unknown_T = " << nUnknown_tot - unknown_q << endl;
//	cout << "unknown_q = " << unknown_q << endl;
//	cout << "unknown_T + unknown_q = " << nUnknown_tot  << endl;
//	return 0;
//}


int unknownIndex_Thermal_2(Element* element, int myid) {

	vector<int> vecInit(6, 0);
	vector<vector<int>>Edge_GBNO(num_element_subdomain, vecInit);
	int num_edge_T_q = LookUp_glb_Edge_Index(element, Edge_GBNO);
	vector<int>Edge_BC(num_edge_T_q, 0);
	vector<int> Edge_Loc(num_edge_T_q, 0);
	vector<int> Edge_Loc_q(num_edge_T_q, 0);
	cout << "num_edge_T_q = " << num_edge_T_q << endl;
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		for (int nn = 0; nn < 4; nn++) {
			int op1 = element[el0].face[nn].opp[1]; int op0 = element[el0].face[nn].opp[0];
			if (op1 == -100) {  //
				for (int ii = 0; ii < 3; ii++) {
					int edgeii = Edge_GBNO[el0][face_edge[nn][ii] - 1];
					Edge_BC[edgeii] = -1;// Dirichlet boundary conditions  of each edge
				}
			}
		}
	}
	vector<int> boundary_element;
	for (int i = 0; i < num_element_subdomain; i++) {
		for (int nn = 0; nn < 4; nn++) {
			int op0 = element[i].face[nn].opp[0]; int ofn = element[i].face[nn].opp[1];
			if (op0 != 0) {
				if (element[i].face[nn].whether_boundary) {
					boundary_element.push_back(i);
					for (int ii = 0; ii < 3; ++ii) {
						int edgeii = Edge_GBNO[i][face_edge[nn][ii] - 1];
						if (Edge_BC[edgeii] != -1)
						{
							Edge_Loc[edgeii] = -1;//not Dirichlet boundary conditions and on the boundary
						}
					}
				}
			}
		}
	}


	int* Node_RTC_T_q = new int[num_node];

	for (int i = 0; i < num_node; i++) {
		Node_RTC_T_q[i] = 0;
	}


	int el, nn, op, ofn, ii, jj, node1, node2, node3, ndm, node_tot, node_ii, ncnt, num_node_Dirichlet, nUnknown_tot, nodeii;
	int ofn1, op1, op2, op3, op4, opd1, opd2, opd3, opd4;
	VectorXi Node_Dirichlet, num_node_SubDomain;
	VectorXi Node_Loc;
	num_unKnown_eb = 0;
	//int num_unKnown_Thermal = 0;
	Node_Dirichlet.resize(num_node);
	Node_Dirichlet.setZero();
	int kdirichlet = 0;

	for (el = 0; el < num_element_subdomain; el++) {
		for (nn = 0; nn < 4; nn++) {
			op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
			if ((op == 0) && (ofn == -100)) {
				node1 = element[el].node[face_node[nn][0] - 1].ver - 1;
				node2 = element[el].node[face_node[nn][1] - 1].ver - 1;
				node3 = element[el].node[face_node[nn][2] - 1].ver - 1;
				kdirichlet++;
				Node_Dirichlet(node1) = -1;
				Node_Dirichlet(node2) = -1;
				Node_Dirichlet(node3) = -1;
			}
		}
	}

	//cout << "kdirichlet = " << kdirichlet << endl;
	//cout << "kconvection = " << kconvection << endl;
	num_node_SubDomain.resize(num_domain);
	Node_Loc.resize(num_node);
	Node_Loc.setZero();


	for (el = 0; el < num_element_subdomain; el++) {
		for (nn = 0; nn < 4; nn++) {
			op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
			if (op != 0) {
				if (element[el].face[nn].whether_boundary) {
					for (ii = 0; ii < 3; ii++) {
						node_ii = element[el].node[face_node[nn][ii] - 1].ver - 1;
						if (Node_Dirichlet(node_ii) != -1) {
							Node_Loc(node_ii) = -1;
						}
					}
				}
			}
		}
	}


	node_tot = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		for (int ii = 0; ii < 6; ii++) {
			int edgeii = Edge_GBNO[el][ii];
			if ((Edge_Loc[edgeii] == 0) && (Edge_BC[edgeii] != -1)) {//not Dirichlet boundary conditions and Not on the boundary
				node_tot += 1;
				Edge_Loc[edgeii] = node_tot;
			}
		}
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver - 1;
			if (Node_Dirichlet(node_ii) != -1) {
				if (Node_Loc(node_ii) == 0) {
					Node_Loc(node_ii) = ++node_tot;//vertice not on the boundary
				}
			}
		}

	}
	int cnt_boundary = 0;
	for (int e = 0; e < boundary_element.size(); e++) {
		int el0 = boundary_element[e];
		for (int ii = 0; ii < 6; ii++) {
			int edgeii = Edge_GBNO[el0][ii];
			if ((Edge_Loc[edgeii] == -1) && (Edge_BC[edgeii] != -1)) {//not Dirichlet and on the boundary
				node_tot += 1;
				Edge_Loc[edgeii] = node_tot;
				++cnt_boundary;

			}
		}
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el0].node[ii].ver - 1;
			if (Node_Dirichlet(node_ii) != -1) {
				if (Node_Loc(node_ii) == -1) {
					Node_Loc(node_ii) = ++node_tot;
					cnt_boundary++;//vertice on the boundary
				}
			}
		}
	}
	ndm = myid;
	Vertice_Boundary_T_q[ndm] = cnt_boundary;
	num_unknown_subdomain_T_q[ndm][0] = node_tot;
	num_unKnown_eb += cnt_boundary;
	cout << "nNode_IN=" << ndm << '\t' << cnt_boundary << '\t' << node_tot << endl;




	for (el = 0; el < num_element_subdomain; el++) {
		//ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			element[el].Eedge_T[ii] = Node_Loc(node_ii - 1);
		}
	}

	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		for (int ii = 0; ii < 6; ii++) {
			int edgeii = Edge_GBNO[el0][ii];
			if (Edge_BC[edgeii] != -1) {
				element[el0].Eedge_T[ii + 4] = Edge_Loc[edgeii];

			}
		}
	}



	for (ndm = myid; ndm < myid + 1; ndm++) {
		ncnt = 0;
		for (el = 0; el < num_element_subdomain; el++) {
			for (nn = 0; nn < 4; nn++) {
				if (element[el].face[nn].whether_boundary) {
					for (ii = 0; ii < 3; ii++) {
						int edge_num = face_edge[nn][ii];
						int edgeii = Edge_GBNO[el][edge_num - 1];
						if (Edge_Loc_q[edgeii] == 0) {
							Edge_Loc_q[edgeii] = ++ncnt;
						}
					}
					for (ii = 0; ii < 3; ii++) {
						nodeii = element[el].node[face_node[nn][ii] - 1].ver - 1;
						if (Node_RTC_T_q[(nodeii)] == 0) {
							ncnt++;
							Node_RTC_T_q[(nodeii)] = ncnt;
						}
						if (Node_Dirichlet(nodeii) == -1) {      //error
							cout << "ff" << endl;
							cout << "op = " << op << endl;
						}
					}
				}
			}
		}
		num_unknown_subdomain_T_q[ndm][1] = ncnt;
		Vertice_Boundary_T_q[ndm] += ncnt;
		num_unKnown_eb += ncnt;
		cout << "nNode_Bc=" << ndm << '\t' << ncnt << endl;
	}


	for (el = 0; el < num_element_subdomain; el++) {
		//ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			element[el].Eedge_q[ii] = Node_RTC_T_q[(node_ii - 1)];
		}
	}
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		for (int ii = 0; ii < 6; ii++) {
			int edgeii = Edge_GBNO[el0][ii];
			element[el0].Eedge_q[ii + 4] = Edge_Loc_q[edgeii];
		}
	}


	nUnknown_tot = 0;
	int unknown_q = 0;
	//int unknown_T = 0;

	ndm = myid;
	nUnknown_tot += num_unknown_subdomain_T_q[ndm][0];
	for (el = 0; el < num_element_subdomain; el++) {
		for (ii = 0; ii < 10; ii++) {
			if (element[el].Eedge_q[ii] != 0) {
				element[el].Eedge_q[ii] += nUnknown_tot;
			}

		}
	}
	nUnknown_tot += num_unknown_subdomain_T_q[ndm][1];
	unknown_q += num_unknown_subdomain_T_q[ndm][1];
	//num_unKnown_Thermal += num_unknown_subdomain_T_q[ndm][0] + num_unknown_subdomain_T_q[ndm][1];


	//cout << "num_unKnown=" << num_unKnown_Thermal << '\t' << nUnknown_tot << endl;
	cout << "unknown_T = " << nUnknown_tot - unknown_q << endl;
	cout << "unknown_q = " << unknown_q << endl;
	cout << "unknown_T + unknown_q = " << nUnknown_tot << endl;
	return 0;
}



int unknownIndex_Thermal(Element* element, int myid) {

	Eigen::MatrixXi Node_RTC_Thermal;
	Node_RTC_Thermal.resize(num_node + 1, num_domain);

	int* Node_RTC_T_q = new int[num_node];

	for (int i = 0; i < num_node; i++) {
		Node_RTC_T_q[i] = 0;
	}

	int el, nn, op, ofn, ii, jj, node1, node2, node3, ndm, node_tot, node_ii, ncnt, num_node_Dirichlet, nUnknown_tot, nodeii;
	int ofn1, op1, op2, op3, op4, opd1, opd2, opd3, opd4;
	VectorXi Node_Dirichlet, num_node_SubDomain;
	MatrixXi Node_Loc;
	num_unKnown_eb = 0;
	int num_unKnown_Thermal = 0;
	Node_Dirichlet.resize(num_node);
	Node_Dirichlet.setZero();
	int kdirichlet = 0;

	for (el = 0; el < num_element_subdomain; el++) {
		for (nn = 0; nn < 4; nn++) {
			op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
			if ((op == 0) && (ofn == -100)) {
				node1 = element[el].node[face_node[nn][0] - 1].ver - 1;
				node2 = element[el].node[face_node[nn][1] - 1].ver - 1;
				node3 = element[el].node[face_node[nn][2] - 1].ver - 1;
				kdirichlet++;
				Node_Dirichlet(node1) = -1;
				Node_Dirichlet(node2) = -1;
				Node_Dirichlet(node3) = -1;
			}
		}
	}

	cout << "kdirichlet = " << kdirichlet << endl;
	//cout << "kconvection = " << kconvection << endl;
	num_node_SubDomain.resize(num_domain);
	Node_Loc.resize(num_node, num_domain);
	Node_Loc.setZero();

	for (ndm = myid; ndm < myid + 1; ndm++) {
		node_tot = 0;
		int cnt_boundary = 0;
		for (el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (nn = 0; nn < 4; nn++) {
					op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
					if (op != 0) {
						if (element[el].face[nn].whether_boundary) {
							for (ii = 0; ii < 3; ii++) {
								node_ii = element[el].node[face_node[nn][ii] - 1].ver - 1;
								if (Node_Dirichlet(node_ii) != -1) {
									Node_Loc(node_ii, ndm) = -1;
								}
							}
						}
					}
				}
			}
			else {
				cout << "wrong domain" << endl;
			}


		}
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (ii = 0; ii < 4; ii++) {
					node_ii = element[el].node[ii].ver - 1;
					if (Node_Dirichlet(node_ii) != -1) {
						if (Node_Loc(node_ii, ndm) == 0) {
							Node_Loc(node_ii, ndm) = ++node_tot;//vertice not on the boundary
						}
					}
				}
			}
			else {
				cout << "wrong domain" << endl;
			}


		}
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (ii = 0; ii < 4; ii++) {
					node_ii = element[el].node[ii].ver - 1;
					if (Node_Dirichlet(node_ii) != -1) {
						if (Node_Loc(node_ii, ndm) == -1) {
							Node_Loc(node_ii, ndm) = ++node_tot;
							cnt_boundary++;//vertice on the boundary
						}
					}
				}
			}
			else {
				cout << "wrong domain" << endl;
			}
		}

		Vertice_Boundary_T_q[ndm] = cnt_boundary;
		num_unknown_subdomain_T_q[ndm][0] = node_tot;
		num_unKnown_eb += cnt_boundary;
		cout << "nNode_IN=" << ndm << '\t' << cnt_boundary << '\t' << node_tot << endl;
	}

	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			element[el].node[ii].unknown_T = Node_Loc(node_ii - 1, ndm - 1);
		}
	}

	Node_RTC_Thermal.setZero();
	for (ndm = myid; ndm < myid + 1; ndm++) {
		ncnt = 0;
		for (el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (nn = 0; nn < 4; nn++) {
					op = element[el].face[nn].opp[0]; ofn = element[el].face[nn].opp[1];
					if (op != 0) {
						if (element[el].face[nn].whether_boundary) {
							for (ii = 0; ii < 3; ii++) {
								nodeii = element[el].node[face_node[nn][ii] - 1].ver - 1;

								//if (element[el].node[face_node[nn][ii] - 1].unknown_q == 0) {    //wrong 
								if (Node_RTC_T_q[(nodeii)] == 0) {
									ncnt++;
									Node_RTC_T_q[(nodeii)] = ncnt;
									//element[el].node[face_node[nn][ii] - 1].unknown_q = ncnt;
								}
								//}
								if (Node_Dirichlet(nodeii) == -1) {      //error
									cout << "ff" << endl;
									cout << "op = " << op << endl;
								}
							}
						}
					}
				}
			}
		}


		num_unknown_subdomain_T_q[ndm][1] = ncnt;
		Vertice_Boundary_T_q[ndm] += ncnt;
		num_unKnown_eb += ncnt;
		cout << "nNode_Bc=" << ndm << '\t' << ncnt << endl;
	}


	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;
		for (ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			element[el].node[ii].unknown_q = Node_RTC_T_q[(node_ii - 1)];
		}
	}

	nUnknown_tot = 0;
	int unknown_q = 0;
	int unknown_T = 0;
	for (ndm = myid; ndm < myid + 1; ndm++) {
		for (el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (ii = 0; ii < 4; ii++) {
					if (element[el].node[ii].unknown_T != 0) {
						element[el].node[ii].unknown_T += nUnknown_tot;
					}
				}
			}
		}
		nUnknown_tot += num_unknown_subdomain_T_q[ndm][0];
		for (el = 0; el < num_element_subdomain; el++) {
			for (ii = 0; ii < 4; ii++) {
				if (element[el].node[ii].unknown_q != 0) {
					element[el].node[ii].unknown_q += nUnknown_tot;
				}

			}
		}
		nUnknown_tot += num_unknown_subdomain_T_q[ndm][1];
		unknown_q += num_unknown_subdomain_T_q[ndm][1];
		num_unKnown_Thermal += num_unknown_subdomain_T_q[ndm][0] + num_unknown_subdomain_T_q[ndm][1];

	}
	cout << "unknown_T = " << nUnknown_tot - unknown_q << endl;
	cout << "unknown_q = " << unknown_q << endl;
	cout << "unknown_T + unknown_q = " << nUnknown_tot << endl;
	return 0;
}



int unknownIndex_stru_older(Element* element, int myid) {
	int  node1, node2, node3, ndm, node_tot, node_ii, ncnt, num_node_Dirichlet, nUnknown_tot, nodeii;
	int* Node_Dirichlet = new int[num_node]();
	VectorXi num_node_SubDomain; num_node_SubDomain.resize(num_domain);
	vector<vector<vector<int>>> Node_Loc(num_node, vector<vector<int>>(num_domain, vector<int>(3, 0)));
	vector<vector<vector<int>>> Node_RTC_stru(num_node+1, vector<vector<int>>(num_domain, vector<int>(3, 0)));
	num_unKnown_stb = 0;
	num_unKnown_stru = 0;
	for (int el = 0; el < num_element_subdomain; el++) {
		for (int nn = 0; nn < 4; nn++) {
			//int ofn = opp_stru[el][nn + 4];

			int ofn = element[el].face[nn].opp_stru[1];
			if (ofn == -70) {
				node1 = element[el].node[face_node[nn][0] - 1].ver - 1;
				node2 = element[el].node[face_node[nn][1] - 1].ver - 1;
				node3 = element[el].node[face_node[nn][2] - 1].ver - 1;


				Node_Dirichlet[node1] = -1;
				Node_Dirichlet[node2] = -1;
				Node_Dirichlet[node3] = -1;
			}
		}
	}

	//cout << "test111" << endl;
	for (ndm = myid; ndm < myid+1; ndm++) {
		node_tot = 0;
		int cnt_boundary = 0;
		for (int el = 0; el < num_element_subdomain; el++) {
			int material_temp = element[el].Material;
			if (material_temp==1|| material_temp==100) {         //revise material
				
				for (int nn = 0; nn < 4; nn++) {
					//int op = opp_stru[el][nn]; int ofn = opp_stru[el][nn + 4];
					int op = element[el].face[nn].opp[0]; int ofn = element[el].face[nn].opp[1];
					if (op != 0) {
						//if (domain[op - 1] - 1 != ndm) {
						if (element[el].face[nn].whether_boundary) {
							for (int ii = 0; ii < 3; ii++) {
								//node_ii = ver(el, face_node[nn][ii] - 1) - 1;
								node_ii = element[el].node[face_node[nn][ii] - 1].ver - 1;
								if (Node_Dirichlet[node_ii] != -1) {
									Node_Loc[node_ii][ndm][0] = -1;
									Node_Loc[node_ii][ndm][1] = -1;
									Node_Loc[node_ii][ndm][2] = -1;
								}
							}
						}
					}
				}
			}
		}
		//cout << "test222" << endl;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int ii = 0; ii < 4; ii++) {
					//node_ii = ver(el, ii) - 1;
					node_ii = element[el].node[ii].ver - 1;
					int material_temp = element[el].Material;
					if (material_temp == 1 || material_temp == 100) {
						if (Node_Dirichlet[node_ii] != -1) {
							if (Node_Loc[node_ii][ndm][0] == 0 && Node_Loc[node_ii][ndm][1] == 0 && Node_Loc[node_ii][ndm][2] == 0) {
								Node_Loc[node_ii][ndm][0] = ++node_tot;
								Node_Loc[node_ii][ndm][1] = ++node_tot;
								Node_Loc[node_ii][ndm][2] = ++node_tot;
							}
						}

					}
					//if (Node_Dirichlet[node_ii] != -1) {
					//	if (Node_Loc[node_ii][ndm][0] == 0 && Node_Loc[node_ii][ndm][1] == 0 && Node_Loc[node_ii][ndm][2] == 0) {
					//		Node_Loc[node_ii][ndm][0] = ++node_tot;
					//		Node_Loc[node_ii][ndm][1] = ++node_tot;
					//		Node_Loc[node_ii][ndm][2] = ++node_tot;
					//	}
					//}
				}
			}
		}
		//cout << "test333" << endl;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int ii = 0; ii < 4; ii++) {
					//node_ii = ver(el, ii) - 1;
					node_ii = element[el].node[ii].ver - 1;
					int material_temp = element[el].Material;
					if (material_temp == 1 || material_temp == 100) {
						if (Node_Dirichlet[node_ii] != -1) {
							if (Node_Loc[node_ii][ndm][0] == -1 && Node_Loc[node_ii][ndm][1] == -1 && Node_Loc[node_ii][ndm][2] == -1) {
								Node_Loc[node_ii][ndm][0] = ++node_tot;
								Node_Loc[node_ii][ndm][1] = ++node_tot;
								Node_Loc[node_ii][ndm][2] = ++node_tot;
								cnt_boundary++; cnt_boundary++; cnt_boundary++;
							}
						}
					}

					//if (Node_Dirichlet[node_ii] != -1) {
					//	if (Node_Loc[node_ii][ndm][0] == -1 && Node_Loc[node_ii][ndm][1] == -1 && Node_Loc[node_ii][ndm][2] == -1) {
					//		Node_Loc[node_ii][ndm][0] = ++node_tot;
					//		Node_Loc[node_ii][ndm][1] = ++node_tot;
					//		Node_Loc[node_ii][ndm][2] = ++node_tot;
					//		cnt_boundary++; cnt_boundary++; cnt_boundary++;
					//	}
					//}
				}
			}
		}
		Node_Boundary_stru[(ndm)] = cnt_boundary;
		num_node_SubDomain(ndm) = node_tot;
		num_unKnown_subdomain_stru[ndm][0] = node_tot;
		num_unKnown_stb += cnt_boundary;
		//cout << "nNode_IN_stru = " << ndm << '\t' << cnt_boundary << '\t' << node_tot << endl;
	}
	for (int el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (int ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			int material_temp = element[el].Material;
			if (material_temp == 1 || material_temp == 100) {
				element[el].node[ii].unknown_uvw[0] = Node_Loc[node_ii - 1][ndm - 1][0];
				element[el].node[ii].unknown_uvw[1] = Node_Loc[node_ii - 1][ndm - 1][1];
				element[el].node[ii].unknown_uvw[2] = Node_Loc[node_ii - 1][ndm - 1][2];
			}

			//element[el].node[ii].unknown_uvw[0] = Node_Loc[node_ii - 1][ndm - 1][0];
			//element[el].node[ii].unknown_uvw[1] = Node_Loc[node_ii - 1][ndm - 1][1];
			//element[el].node[ii].unknown_uvw[2] = Node_Loc[node_ii - 1][ndm - 1][2];

			//ver(el, ii + 8) = Node_Loc[node_ii - 1][ndm - 1][0];
			//ver(el, ii + 12) = Node_Loc[node_ii - 1][ndm - 1][1];
			//ver(el, ii + 16) = Node_Loc[node_ii - 1][ndm - 1][2];
		}
	}
	//cout << "test444" << endl;
	ndm = myid;
	for (ndm = myid; ndm < myid+1; ndm++) {
		ncnt = 0;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int nn = 0; nn < 4; nn++) {
					//int op = opp_stru[el][nn]; int ofn = opp_stru[el][nn + 4];
					int op = element[el].face[nn].opp[0]; int ofn = element[el].face[nn].opp[1];

					if (op != 0) {
						if (element[el].face[nn].whether_boundary) {
							for (int ii = 0; ii < 3; ii++) {
								//nodeii = ver(el, face_node[nn][ii] - 1) - 1;
								nodeii = element[el].node[face_node[nn][ii] - 1].ver - 1;
								int material_temp = element[el].Material;
								if (material_temp == 1 || material_temp == 100) {
									if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
										Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
										Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
										Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
									}
								}


								//if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
								//	Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
								//	Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
								//	Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
								//}




								//if (Node_Dirichlet[nodeii] != -1) {
								//	if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
								//		Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
								//		Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
								//		Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
								//	}
								//}
							}
						}
					}
				}
			}
		}
		//Node_RTC_stru[num_node][ndm][0] = ncnt;
		num_unKnown_subdomain_stru[ndm][1] = ncnt;
		//num_unknown_subdomain_uvw
		Node_Boundary_stru[(ndm)] += ncnt;
		num_unKnown_stb += ncnt;
		cout << "nNode_Bc_stru = " << ndm << '\t' << ncnt << endl;
	}
	//cout << "test6666" << endl;
	nUnknown_tot = 0;
	int unknown_q = 0;
	int unknown_T = 0;
	for (ndm = myid; ndm < myid+1; ndm++) {
		//num_unKnown_subdomain_stru[ndm][0] = num_node_SubDomain(ndm) + Node_RTC_stru[num_node][ndm][0];
		nUnknown_tot += num_node_SubDomain(ndm);
		for (int el = 0; el < num_node; el++) {
			if (Node_RTC_stru[el][ndm][0] != 0 && Node_RTC_stru[el][ndm][1] != 0 && Node_RTC_stru[el][ndm][2] != 0) {
				Node_RTC_stru[el][ndm][0] += nUnknown_tot;
				Node_RTC_stru[el][ndm][1] += nUnknown_tot;
				Node_RTC_stru[el][ndm][2] += nUnknown_tot;
			}
		}
		nUnknown_tot += num_unKnown_subdomain_stru[ndm][1];
		unknown_q += num_unKnown_subdomain_stru[ndm][1];
		num_unKnown_stru += num_node_SubDomain(ndm) + Node_RTC_stru[num_node][ndm][0];
	}

	for (int el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (int ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			int material_temp = element[el].Material;
			if (material_temp == 1 || material_temp == 100) {
				element[el].node[ii].unknown_abc[0] = Node_RTC_stru[node_ii - 1][ndm - 1][0];
				element[el].node[ii].unknown_abc[1] = Node_RTC_stru[node_ii - 1][ndm - 1][1];
				element[el].node[ii].unknown_abc[2] = Node_RTC_stru[node_ii - 1][ndm - 1][2];
			}

			//element[el].node[ii].unknown_abc[0] = Node_RTC_stru[node_ii - 1][ndm - 1][0];
			//element[el].node[ii].unknown_abc[1] = Node_RTC_stru[node_ii - 1][ndm - 1][1];
			//element[el].node[ii].unknown_abc[2] = Node_RTC_stru[node_ii - 1][ndm - 1][2];

		}
	}

	cout << "unknown_uvw=" << nUnknown_tot- unknown_q << endl;
	cout << "unknown_abc=" << unknown_q << endl;




	//cout << nUnknown_tot * 3 << endl;
	return 0;
}



int unknownIndex_stru(Element* element, int myid) {
	int  node1, node2, node3, ndm, node_tot, node_ii, ncnt, num_node_Dirichlet, nUnknown_tot, nodeii;
	int* Node_Dirichlet = new int[num_node]();
	VectorXi num_node_SubDomain; num_node_SubDomain.resize(num_domain);
	vector<vector<vector<int>>> Node_Loc(num_node, vector<vector<int>>(num_domain, vector<int>(3, 0)));
	vector<vector<vector<int>>> Node_RTC_stru(num_node + 1, vector<vector<int>>(num_domain, vector<int>(3, 0)));
	num_unKnown_stb = 0;
	num_unKnown_stru = 0;

	enum Boundary { DIE1, DIE2, DIE3, DIE4, DIE5, DIE6, DIE7, DIE8, DIE9, DIE10, DIE11, DIE12, DIE13, DIE14, DIE15, DIE16 };
	vector<int>  DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info, DIE6_info, DIE7_info, DIE8_info,
		DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info, bd_flag;


	Read_BounaryInfo2(DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info, DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info,
		DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);










	for (int el = 0; el < num_element_subdomain; el++) {
		for (int nn = 0; nn < 4; nn++) {
			//int ofn = opp_stru[el][nn + 4];

			int ofn = element[el].face[nn].opp_stru[1];
			if (ofn == -70) {
				node1 = element[el].node[face_node[nn][0] - 1].ver - 1;
				node2 = element[el].node[face_node[nn][1] - 1].ver - 1;
				node3 = element[el].node[face_node[nn][2] - 1].ver - 1;


				Node_Dirichlet[node1] = -1;
				Node_Dirichlet[node2] = -1;
				Node_Dirichlet[node3] = -1;
			}
		}
	}

	//cout << "test111" << endl;
	for (ndm = myid; ndm < myid + 1; ndm++) {
		node_tot = 0;
		int cnt_boundary = 0;
		for (int el = 0; el < num_element_subdomain; el++) {

			int material_temp = element[el].Material;
			bool flag = 0;
			for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
				flag = 1;
			}
			//if (material_temp == Material_1 || material_temp == Material_2) {         //revise material
			if (flag) {
				for (int nn = 0; nn < 4; nn++) {
					//int op = opp_stru[el][nn]; int ofn = opp_stru[el][nn + 4];
					int op = element[el].face[nn].opp[0]; int ofn = element[el].face[nn].opp[1];
					if (op != 0) {
						//if (domain[op - 1] - 1 != ndm) {
						if (element[el].face[nn].whether_boundary) {
							for (int ii = 0; ii < 3; ii++) {
								//node_ii = ver(el, face_node[nn][ii] - 1) - 1;
								node_ii = element[el].node[face_node[nn][ii] - 1].ver - 1;
								if (Node_Dirichlet[node_ii] != -1) {
									Node_Loc[node_ii][ndm][0] = -1;
									Node_Loc[node_ii][ndm][1] = -1;
									Node_Loc[node_ii][ndm][2] = -1;
								}
							}
						}
					}
				}
			}
		}
		//cout << "test222" << endl;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int ii = 0; ii < 4; ii++) {
					//node_ii = ver(el, ii) - 1;
					node_ii = element[el].node[ii].ver - 1;
					int material_temp = element[el].Material;
					bool flag = 0;
					for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
						flag = 1;
					}
					//if (material_temp == Material_1 || material_temp == Material_2) {
					if (flag) {
						if (Node_Dirichlet[node_ii] != -1) {
							if (Node_Loc[node_ii][ndm][0] == 0 && Node_Loc[node_ii][ndm][1] == 0 && Node_Loc[node_ii][ndm][2] == 0) {
								Node_Loc[node_ii][ndm][0] = ++node_tot;
								Node_Loc[node_ii][ndm][1] = ++node_tot;
								Node_Loc[node_ii][ndm][2] = ++node_tot;
							}
						}

					}
					//if (Node_Dirichlet[node_ii] != -1) {
					//	if (Node_Loc[node_ii][ndm][0] == 0 && Node_Loc[node_ii][ndm][1] == 0 && Node_Loc[node_ii][ndm][2] == 0) {
					//		Node_Loc[node_ii][ndm][0] = ++node_tot;
					//		Node_Loc[node_ii][ndm][1] = ++node_tot;
					//		Node_Loc[node_ii][ndm][2] = ++node_tot;
					//	}
					//}
				}
			}
		}
		//cout << "test333" << endl;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int ii = 0; ii < 4; ii++) {
					//node_ii = ver(el, ii) - 1;
					node_ii = element[el].node[ii].ver - 1;
					int material_temp = element[el].Material;
					bool flag = 0;
					for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
						flag = 1;
					}
					//if (material_temp == Material_1 || material_temp == Material_2) {
					if (flag) {
						if (Node_Dirichlet[node_ii] != -1) {
							if (Node_Loc[node_ii][ndm][0] == -1 && Node_Loc[node_ii][ndm][1] == -1 && Node_Loc[node_ii][ndm][2] == -1) {
								Node_Loc[node_ii][ndm][0] = ++node_tot;
								Node_Loc[node_ii][ndm][1] = ++node_tot;
								Node_Loc[node_ii][ndm][2] = ++node_tot;
								cnt_boundary++; cnt_boundary++; cnt_boundary++;
							}
						}
					}

					//if (Node_Dirichlet[node_ii] != -1) {
					//	if (Node_Loc[node_ii][ndm][0] == -1 && Node_Loc[node_ii][ndm][1] == -1 && Node_Loc[node_ii][ndm][2] == -1) {
					//		Node_Loc[node_ii][ndm][0] = ++node_tot;
					//		Node_Loc[node_ii][ndm][1] = ++node_tot;
					//		Node_Loc[node_ii][ndm][2] = ++node_tot;
					//		cnt_boundary++; cnt_boundary++; cnt_boundary++;
					//	}
					//}
				}
			}
		}
		Node_Boundary_stru[(ndm)] = cnt_boundary;
		num_node_SubDomain(ndm) = node_tot;
		num_unKnown_subdomain_stru[ndm][0] = node_tot;
		num_unKnown_stb += cnt_boundary;
		//cout << "nNode_IN_stru = " << ndm << '\t' << cnt_boundary << '\t' << node_tot << endl;
	}
	for (int el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (int ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			int material_temp = element[el].Material;
			bool flag = 0;
			for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
				flag = 1;
			}
			//if (material_temp == Material_1 || material_temp == Material_2) {
			if (flag) {
				element[el].node[ii].unknown_uvw[0] = Node_Loc[node_ii - 1][ndm - 1][0];
				element[el].node[ii].unknown_uvw[1] = Node_Loc[node_ii - 1][ndm - 1][1];
				element[el].node[ii].unknown_uvw[2] = Node_Loc[node_ii - 1][ndm - 1][2];
			}

			//element[el].node[ii].unknown_uvw[0] = Node_Loc[node_ii - 1][ndm - 1][0];
			//element[el].node[ii].unknown_uvw[1] = Node_Loc[node_ii - 1][ndm - 1][1];
			//element[el].node[ii].unknown_uvw[2] = Node_Loc[node_ii - 1][ndm - 1][2];

			//ver(el, ii + 8) = Node_Loc[node_ii - 1][ndm - 1][0];
			//ver(el, ii + 12) = Node_Loc[node_ii - 1][ndm - 1][1];
			//ver(el, ii + 16) = Node_Loc[node_ii - 1][ndm - 1][2];
		}
	}
	//cout << "test444" << endl;
	ndm = myid;
	for (ndm = myid; ndm < myid + 1; ndm++) {
		ncnt = 0;
		for (int el = 0; el < num_element_subdomain; el++) {
			if (element[el].domain - 1 == ndm) {
				for (int nn = 0; nn < 4; nn++) {
					//int op = opp_stru[el][nn]; int ofn = opp_stru[el][nn + 4];
					int op = element[el].face[nn].opp[0]; int ofn = element[el].face[nn].opp[1];

					if (op != 0) {
						if (element[el].face[nn].whether_boundary) {
							for (int ii = 0; ii < 3; ii++) {
								//nodeii = ver(el, face_node[nn][ii] - 1) - 1;
								nodeii = element[el].node[face_node[nn][ii] - 1].ver - 1;
								int material_temp = element[el].Material;
								bool flag = 0;
								for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
									flag = 1;
								}
								//if (material_temp == Material_1 || material_temp == Material_2) {
								if (flag) {
									if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
										Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
										Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
										Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
									}
								}


								//if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
								//	Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
								//	Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
								//	Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
								//}




								//if (Node_Dirichlet[nodeii] != -1) {
								//	if (Node_RTC_stru[nodeii][ndm][0] == 0 && Node_RTC_stru[nodeii][ndm][1] == 0 && Node_RTC_stru[nodeii][ndm][2] == 0) {
								//		Node_RTC_stru[nodeii][ndm][0] = ++ncnt;
								//		Node_RTC_stru[nodeii][ndm][1] = ++ncnt;
								//		Node_RTC_stru[nodeii][ndm][2] = ++ncnt;
								//	}
								//}
							}
						}
					}
				}
			}
		}
		//Node_RTC_stru[num_node][ndm][0] = ncnt;
		num_unKnown_subdomain_stru[ndm][1] = ncnt;
		//num_unknown_subdomain_uvw
		Node_Boundary_stru[(ndm)] += ncnt;
		num_unKnown_stb += ncnt;
		cout << "nNode_Bc_stru = " << ndm << '\t' << ncnt << endl;
	}
	//cout << "test6666" << endl;
	nUnknown_tot = 0;
	int unknown_q = 0;
	int unknown_T = 0;
	for (ndm = myid; ndm < myid + 1; ndm++) {
		//num_unKnown_subdomain_stru[ndm][0] = num_node_SubDomain(ndm) + Node_RTC_stru[num_node][ndm][0];
		nUnknown_tot += num_node_SubDomain(ndm);
		for (int el = 0; el < num_node; el++) {
			if (Node_RTC_stru[el][ndm][0] != 0 && Node_RTC_stru[el][ndm][1] != 0 && Node_RTC_stru[el][ndm][2] != 0) {
				Node_RTC_stru[el][ndm][0] += nUnknown_tot;
				Node_RTC_stru[el][ndm][1] += nUnknown_tot;
				Node_RTC_stru[el][ndm][2] += nUnknown_tot;
			}
		}
		nUnknown_tot += num_unKnown_subdomain_stru[ndm][1];
		unknown_q += num_unKnown_subdomain_stru[ndm][1];
		num_unKnown_stru += num_node_SubDomain(ndm) + Node_RTC_stru[num_node][ndm][0];
	}

	for (int el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (int ii = 0; ii < 4; ii++) {
			node_ii = element[el].node[ii].ver;
			int material_temp = element[el].Material;

			bool flag = 0;
			for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material_temp == *beg1) {
				flag = 1;
			}
			//if (material_temp == Material_1 || material_temp == Material_2) {
			if (flag) {
				element[el].node[ii].unknown_abc[0] = Node_RTC_stru[node_ii - 1][ndm - 1][0];
				element[el].node[ii].unknown_abc[1] = Node_RTC_stru[node_ii - 1][ndm - 1][1];
				element[el].node[ii].unknown_abc[2] = Node_RTC_stru[node_ii - 1][ndm - 1][2];
			}

			//element[el].node[ii].unknown_abc[0] = Node_RTC_stru[node_ii - 1][ndm - 1][0];
			//element[el].node[ii].unknown_abc[1] = Node_RTC_stru[node_ii - 1][ndm - 1][1];
			//element[el].node[ii].unknown_abc[2] = Node_RTC_stru[node_ii - 1][ndm - 1][2];

		}
	}

	cout << "unknown_uvw=" << nUnknown_tot - unknown_q << endl;
	cout << "unknown_abc=" << unknown_q << endl;




	//cout << nUnknown_tot * 3 << endl;
	return 0;
}






int Read_BounaryInfo2(vector<int>& DIE1_info, vector<int>& DIE2_info, vector<int>& DIE3_info,
	vector<int>& DIE4_info,
	vector<int>& DIE5_info,
	vector<int>& DIE6_info,
	vector<int>& DIE7_info,
	vector<int>& DIE8_info,
	vector<int>& DIE9_info,
	vector<int>& DIE10_info,
	vector<int>& DIE11_info,
	vector<int>& DIE12_info,
	vector<int>& DIE13_info,
	vector<int>& DIE14_info,
	vector<int>& DIE15_info,
	vector<int>& DIE16_info)
{

	ifstream infile;
	string x;
	infile.open("material_T_F.txt");
	string sta;

	while (getline(infile, x) && (!x.empty()))
	{
		if (sta.empty()) {
			sta = x;
		}
		else {
			string pmc;       //Variable name does not contain physical meaning
			for (auto beg = x.cbegin(); beg != x.cend(); ++beg) {
				if (*beg != ' ') {
					if (*beg != ',') {
						pmc = pmc + *beg;
					}
					else {
						s2i2(pmc, sta, DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info,
							DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);
						pmc.clear();
					}
				}
			}
			s2i2(pmc, sta, DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info,
				DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);
			sta.clear();
		}
	}
	infile.close();
	return 0;
}


bool str_compare2(string str1, string str2) {
	int len = str2.length();
	for (int i = 0; i < len; i++) {
		if (str1[i] != str2[i]) {
			return false;
		}
	}
	return true;
}

void s2i2(const string& info, const string& sta,
	vector<int>& DIE1_info,
	vector<int>& DIE2_info,
	vector<int>& DIE3_info,
	vector<int>& DIE4_info,
	vector<int>& DIE5_info,
	vector<int>& DIE6_info,
	vector<int>& DIE7_info,
	vector<int>& DIE8_info,
	vector<int>& DIE9_info,
	vector<int>& DIE10_info,
	vector<int>& DIE11_info,
	vector<int>& DIE12_info,
	vector<int>& DIE13_info,
	vector<int>& DIE14_info,
	vector<int>& DIE15_info,
	vector<int>& DIE16_info)
{
	int i = 0; string temp;
	std::string::size_type sz;   // alias of size_t
	int i_dec = std::stoi(info, &sz);
	int idec_end = 0;

	while (sz + i + 1 < info.length()) {
		temp += info[sz + i + 1];
		++i;
	}
	if (i)  idec_end = std::stoi(temp);
	else   idec_end = i_dec;
	for (int j = i_dec; j <= idec_end; ++j) {
		if (str_compare2(sta, "DOMAINF")) DIE1_info.push_back(j);
		else if (str_compare2(sta, "DOMAINT"))  DIE2_info.push_back(j);
	}
	temp.clear();
	sz = 0; idec_end = 0; i_dec = 0;

}