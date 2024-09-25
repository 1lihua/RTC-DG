#include"Element.h"
#include"Global_Data.h"
#include<vector>
#include<fstream>
#include"time.h"



using namespace std;
typedef struct Sum {
	size_t idx;
	size_t sum;
	//按math从大到小排序
	inline bool operator < (const Sum& x) const {
		return sum < x.sum;
	}

};
typedef struct EdgePro {
	size_t idx;
	size_t  pro;
	//按math从大到小排序
	inline bool operator < (const EdgePro& x) const {
		return pro < x.pro;
	}

};

//#define Material_E1 1   
constexpr auto Material_E1 = 1;
constexpr auto Material_E2 = 2   ;

int LookUp_glbIndex(Element* el, vector<vector<int>>& Edge_GBNO) {

	num_edge = -1;
	vector<Sum> opp_Sum((size_t)num_element_subdomain * 6);
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
		vector<EdgePro> v2;
		while (sum1 == (*beg).sum) {
			int ii = ((*beg).idx) % 6; int ith = ((*beg).idx) / 6; EdgePro temp1;
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
			if ((Edge_GBNO[ith][nn] == 0) && (Edge_GBNO[ith][nn + 6] == 0)) { Edge_GBNO[ith][nn] = ++num_edge;  Edge_GBNO[ith][nn + 6] = (++num_edge); }
			//if (Edge_GBNO[ith][nn + 6] == 0) { Edge_GBNO[ith][nn] = ++num_edge;  Edge_GBNO[ith][nn + 6] = (++num_edge); }
#pragma omp parallel for
			for (int j = i + 1; j < v2Size; ++j) {
				if (v2[i].pro == v2[j].pro) {
					int op1 = v2[j].idx / 6; int ofn1 = v2[j].idx % 6;
					Edge_GBNO[op1][ofn1] = Edge_GBNO[ith][nn];
					Edge_GBNO[op1][ofn1 + 6] = Edge_GBNO[ith][nn + 6];
				}
			}
		}
	}

	num_edge++;

	cout << "num_edge _E_H = " << num_edge << endl;
	int node0 = 0, node1 = 0;

#pragma omp for
	for (int ith = 0; ith < num_element_subdomain; ++ith) {
		for (int ii = 0; ii < 6; ++ii) {
			node0 = edge_node_local[ii][0] - 1; node1 = edge_node_local[ii][1] - 1;
			if (el[ith].node[node0].ver < el[ith].node[node1].ver) { Edge_GBNO[ith][ii + 12] = -1; }//vector basis function Nlk l<k
		}

	}
	return 0;//added by song

}

int LookUp_glbIndexE(Element* el, vector<vector<int>>& Edge_GBNO) {
	int op0, op1, edgeii, edge_tot, ndm;
	int face_edge[4][6] = { { 0, 1, 3, 6, 7, 9 }, { 1, 2, 5, 7, 8, 11 }, { 0, 2, 4, 6, 8, 10 }, { 3, 4, 5, 9, 10, 11 } }; //Four faces represented by six edges
	vector<int>Edge_BC(num_edge, 0);
	vector<vector<int>> Edge_Loc(num_edge, vector<int>(num_domain, 0));
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		for (int nn = 0; nn < 4; nn++) {
			op1 = el[el0].face[nn].opp[1]; op0 = el[el0].face[nn].opp[0];
			if (op1 == -5) {  //PEC
				for (int ii = 0; ii < 6; ii++) {
					edgeii = Edge_GBNO[el0][face_edge[nn][ii]];
					Edge_BC[edgeii] = -1;//PEC Boundary condition of each edge
				}
			}
		}
	}
	int cnt_boundary = 0;
	for (ndm = 0; ndm < num_domain; ndm++) {
		edge_tot = 0;//Number of sides on each subdomain(Not on the boundary)
		cnt_boundary = 0;//The number of edges at the sub domain interface
		vector<int> boundary_element;

		for (int i = 0; i < num_element_subdomain; i++) {
			if (el[i].domain - 1 == ndm) {
				for (int nn = 0; nn < 4; nn++) {
					op0 = el[i].face[nn].opp[0]; int ofn = el[i].face[nn].opp[1];
					if (op0 != 0) {
						if (el[i].face[nn].whether_boundary) {
							boundary_element.push_back(i);
							int material_temp = el[i].Material;
							if (material_temp == Material_E1 || material_temp == Material_E2||1) {
								for (int ii = 0; ii < 6; ++ii) {
									edgeii = Edge_GBNO[i][face_edge[nn][ii]];
									if (Edge_BC[edgeii] != -1)
									{
										Edge_Loc[edgeii][ndm] = -1;//not PEC and on the boundary
									}
								}
							}
							//for (int ii = 0; ii < 6; ++ii) {
							//	edgeii = Edge_GBNO[i][face_edge[nn][ii]];
							//	if (Edge_BC[edgeii] != -1)
							//	{
							//		Edge_Loc[edgeii][ndm] = -1;//not PEC and on the boundary
							//	}
							//}
						}
					}
				}
			}
		}
		for (int el0 = 0; el0 < num_element_subdomain; el0++) {
			if (el[el0].domain == ndm + 1) {
				for (int ii = 0; ii < 12; ii++) {
					edgeii = Edge_GBNO[el0][ii];
					int material_temp = el[el0].Material;
					if ((Edge_Loc[edgeii][ndm] == 0) && (Edge_BC[edgeii] != -1)&&(material_temp == Material_E1 || material_temp == Material_E2 || 1)) {//not pec and Not on the boundary
						edge_tot += 1;
						Edge_Loc[edgeii][ndm] = edge_tot;
					}
					//if ((Edge_Loc[edgeii][ndm] == 0) && (Edge_BC[edgeii] != -1)) {//not pec and Not on the boundary
					//	edge_tot += 1;
					//	Edge_Loc[edgeii][ndm] = edge_tot;
					//}
				}
			}
		}
		for (int e = 0; e < boundary_element.size(); e++) {
			int el0 = boundary_element[e];
			for (int ii = 0; ii < 12; ii++) {
				edgeii = Edge_GBNO[el0][ii];
				int material_temp = el[el0].Material;
				if ((Edge_Loc[edgeii][ndm] == -1) && (Edge_BC[edgeii] != -1) && (material_temp == Material_E1 || material_temp == Material_E2 || 1)) {//not pec and on the boundary
					edge_tot += 1;
					Edge_Loc[edgeii][ndm] = edge_tot;
					++cnt_boundary;
				}
				//if ((Edge_Loc[edgeii][ndm] == -1) && (Edge_BC[edgeii] != -1)) {//not pec and on the boundary
				//	edge_tot += 1;
				//	Edge_Loc[edgeii][ndm] = edge_tot;
				//	++cnt_boundary;
				//}
			}

		}
		Edge_Boundary[ndm] = cnt_boundary; //The total number of non PEC edges in each subdomain boundary
		//cout << "The total number of non PEC edges in each subdomain boundary is " << cnt_boundary << endl;
		num_unknown_subdomain[ndm][0] = edge_tot;//The total number of edges in each subdomain
		num_unKnown_b += cnt_boundary;

		int count_temp = 0;
		for (int i = 0; i < num_element_subdomain; i++) {
			if (el[i].domain - 1 == ndm) {
				for (int nn = 0; nn < 4; nn++) {
					op0 = el[i].face[nn].opp[0];
					if (op0 != 0) {
						if (el[i].face[nn].whether_boundary) {
							for (int ii = 0; ii < 6; ++ii) {
								edgeii = Edge_GBNO[i][face_edge[nn][ii]];
								if (Edge_Loc[edgeii][ndm] == -1)  //not PMC
								{
									count_temp++;
								}
							}

						}
					}
				}
			}
		}

		//cout << "The total number of non PEC edges in each subdomain boundary is(count_temp) " << count_temp<< endl;


	}
	for (int i = 0; i < num_element_subdomain; i++) {
		for (int j = 0; j < 24; j++) {
			el[i].Eedge_GBNO[j] = 0;
		}
	}
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		ndm = el[el0].domain;
		for (int ii = 0; ii < 12; ii++) {
			edgeii = Edge_GBNO[el0][ii];
			int material_temp = el[el0].Material;
			if (Edge_BC[edgeii] != -1&&(material_temp == Material_E1 || material_temp == Material_E2 || 1)) {
				el[el0].Eedge_GBNO[ii] = Edge_Loc[edgeii][ndm - 1];
				el[el0].Eedge_GBNO[ii + 12] = Edge_GBNO[el0][ii + 12];
			}

			//if (Edge_BC[edgeii] != -1) {
			//	el[el0].Eedge_GBNO[ii] = Edge_Loc[edgeii][ndm - 1];
			//	el[el0].Eedge_GBNO[ii + 12] = Edge_GBNO[el0][ii + 12];
			//}
		}
	}




	return 0;
}

int LookUp_glbIndexH(Element* el, vector<vector<int>>& Edge_GBNO) {
	int op0, op1, edgeii, edge_tot, ndm;
	int face_edge[4][6] = { { 0, 1, 3, 6, 7, 9 }, { 1, 2, 5, 7, 8, 11 }, { 0, 2, 4, 6, 8, 10 }, { 3, 4, 5, 9, 10, 11 } };
	Eigen::VectorXi Edge_BC(num_edge);
	Eigen::MatrixXi Edge_Loc(num_edge, num_domain);
	Edge_BC.setZero();
	Edge_Loc.setZero();
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		//if (el[el0].Material == 2) {
		//	continue;
		//}
		double zt = (el[el0].node[0].zb[2] + el[el0].node[1].zb[2] + el[el0].node[2].zb[2] + el[el0].node[3].zb[2]) / 4.0;
		//if (el[el0].Material == wave_material) {
		if (zt > wave_zzz) {
			continue;
		}
		for (int nn = 0; nn < 4; nn++) {
			op1 = el[el0].face[nn].opp[1]; op0 = el[el0].face[nn].opp[0];
			if ((-215 <= op1 && op1 <= -200)) {
				for (int ii = 0; ii < 6; ii++) {
					edgeii = Edge_GBNO[el0][face_edge[nn][ii]];
					Edge_BC[edgeii] = -1;//PMC
				}
			}
		}
	}
	mesh_at_RTC = (int*)malloc(sizeof(int) * num_element_subdomain);
	for (int i = 0; i < num_element_subdomain; i++) {
		mesh_at_RTC[i] = 0;
	}
	for (int ndm = 0; ndm < num_domain; ndm++) {
		edge_tot = 0;
		for (int i = 0; i < num_element_subdomain; i++) {
			if (el[i].domain - 1 == ndm) {
				for (int nn = 0; nn < 4; nn++) {
					op0 = el[i].face[nn].opp[0];
					if (op0 != 0) {
						if (el[i].face[nn].whether_boundary) {
							mesh_at_RTC[i] = ndm + 1;
							for (int ii = 0; ii < 6; ++ii) {
								edgeii = Edge_GBNO[i][face_edge[nn][ii]];
								if (Edge_BC[edgeii] != -1)  //not PMC
								{
									int material_temp = el[i].Material;

									if (Edge_Loc(edgeii, ndm) == 0 && (material_temp == Material_E1 || material_temp == Material_E2 || 1)) {
										edge_tot = edge_tot + 1;//edge_tot++
										Edge_Loc(edgeii, ndm) = edge_tot;
									}
									//if (Edge_Loc(edgeii, ndm) == 0) {
									//	edge_tot = edge_tot + 1;//edge_tot++
									//	Edge_Loc(edgeii, ndm) = edge_tot;
									//}
								}
							}
						}
					}
				}
			}
		}
		Edge_Boundary[ndm] += edge_tot;
		//cout << "The total number of non PMC edges in each subdomain boundary is " << edge_tot << endl;
		num_unknown_subdomain[ndm][1] = edge_tot;
		num_unKnown_b += edge_tot;
	}
	for (int i = 0; i < num_element_subdomain; i++) {
		for (int j = 0; j < 24; j++) {
			el[i].Hedge_GBNO[j] = 0;
		}
	}
//#	pragma omp parallel for num_threads(112) 
#	pragma omp parallel for
	for (int el0 = 0; el0 < num_element_subdomain; el0++) {
		for (int nn = 0; nn < 4; nn++) {
			if (el[el0].face[nn].opp[0] != 0) {
				if (el[el0].face[nn].whether_boundary) {
					for (int ii = 0; ii < 6; ii++) {
						int material_temp = el[el0].Material;
						if (Edge_BC[Edge_GBNO[el0][face_edge[nn][ii]]] != -1&& (material_temp == Material_E1 || material_temp == Material_E2 || 1)) {
							el[el0].Hedge_GBNO[face_edge[nn][ii]] = Edge_Loc(Edge_GBNO[el0][face_edge[nn][ii]], el[el0].domain - 1);
							el[el0].Hedge_GBNO[face_edge[nn][ii] + 12] = Edge_GBNO[el0][face_edge[nn][ii] + 12];
						}
						//if (Edge_BC[Edge_GBNO[el0][face_edge[nn][ii]]] != -1) {
						//	el[el0].Hedge_GBNO[face_edge[nn][ii]] = Edge_Loc(Edge_GBNO[el0][face_edge[nn][ii]], el[el0].domain - 1);
						//	el[el0].Hedge_GBNO[face_edge[nn][ii] + 12] = Edge_GBNO[el0][face_edge[nn][ii] + 12];
						//}
					}
				}
			}
		}
	}
	return 0;//added by song
}



int EdgeIndex_Finalized(Element* el) {//Arrange the numbers of all subdomain edges in order
	int nUnknown_tot = 0;
	num_unknown = 0;
	int nUnknown_tot_E = 0;
	int nUnknown_tot_H = 0;
	int nUnknown_boundary = 0;
	for (int ndm = 0; ndm < num_domain; ndm++) {
		for (int el0 = 0; el0 < num_element_subdomain; el0++) {
			if (el[el0].domain == (1 + ndm)) {
				for (int ii = 0; ii < 12; ii++) {
					if (el[el0].Eedge_GBNO[ii] != 0)
						el[el0].Eedge_GBNO[ii] += nUnknown_tot;
				}
			}
		}
		nUnknown_tot += num_unknown_subdomain[ndm][0];
		nUnknown_tot_E += num_unknown_subdomain[ndm][0];
		for (int el0 = 0; el0 < num_element_subdomain; el0++) {
			if (mesh_at_RTC[el0] == (ndm + 1)) {
				for (int ii = 0; ii < 12; ii++) {
					if (el[el0].Hedge_GBNO[ii] != 0)
						el[el0].Hedge_GBNO[ii] += nUnknown_tot;
				}
			}
		}
		nUnknown_tot += num_unknown_subdomain[ndm][1];
		nUnknown_tot_H += num_unknown_subdomain[ndm][1];
		num_unknown = num_unknown + num_unknown_subdomain[ndm][0] + num_unknown_subdomain[ndm][1];
		nUnknown_boundary += Edge_Boundary[ndm];
	}
	cout << "Total number of unknowns is:     " << nUnknown_tot << "   " << num_unknown << endl;
	cout << "Total number of unknowns_E is:     " << nUnknown_tot_E << "   " << nUnknown_tot_E << endl;
	cout << "Total number of unknowns_H is:     " << nUnknown_tot_H << "   " << nUnknown_tot_H << endl;
	cout << "Total number of unknowns_Boundary is:     " << nUnknown_boundary << "   " << nUnknown_tot_H << endl;
	return 0;

}

int LookUp_EdgeIndex(Element* el_subdomain) {

	vector<int> vecInit(24, 0);
	for (int i = 12; i < 24; ++i) {
		vecInit[i] = 1;
	}
	vector<vector<int>>Edge_GBNO(num_element_subdomain, vecInit);

	Edge_Boundary = (int*)malloc(sizeof(int) * num_domain);
	num_unknown_subdomain = (int**)malloc(sizeof(int*) * num_domain);
	for (int i = 0; i < num_domain; i++) { num_unknown_subdomain[i] = (int*)malloc(sizeof(int) * 2); }

	for (int i = 0; i < num_domain; i++) {
		for (int j = 0; j < 2; j++) {
			num_unknown_subdomain[i][j] = 0;
		}
		Edge_Boundary[i] = 0;
	}

	LookUp_glbIndex(el_subdomain, Edge_GBNO);
	cout << "Total number of physical edges is =   " << num_edge << endl;
	LookUp_glbIndexE(el_subdomain, Edge_GBNO);
	LookUp_glbIndexH(el_subdomain, Edge_GBNO);
	EdgeIndex_Finalized(el_subdomain);



	Edge_GBNO.clear();
	return 0;
}
