
#include "Element.h"
#include "Global_Data.h"
#include"ReadFile.h"
#include <algorithm>
#include <numeric>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include"time.h"
#include"Global_Data.h"
#include<omp.h>
using namespace std;
int Read_BounaryInfo(vector<int>& pec_info, vector<int>& pmc_info, vector<int>& abc_info, vector<int>& lump1_info, vector<int>& lump2_info, vector<int>& lump3_info,
	vector<int>& lump4_info,
	vector<int>& lump5_info,
	vector<int>& lump6_info,
	vector<int>& lump7_info,
	vector<int>& lump8_info,
	vector<int>& lump9_info,
	vector<int>& lump10_info,
	vector<int>& lump11_info,
	vector<int>& lump12_info,
	vector<int>& lump13_info,
	vector<int>& lump14_info,
	vector<int>& lump15_info,
	vector<int>& lump16_info,
	vector<int>& wave1_info, vector<int>& wave2_info, vector<int>& wave3_info,
	vector<int>& wave4_info,
	vector<int>& wave5_info,
	vector<int>& wave6_info,
	vector<int>& wave7_info,
	vector<int>& wave8_info,
	vector<int>& wave9_info,
	vector<int>& wave10_info,
	vector<int>& wave11_info,
	vector<int>& wave12_info,
	vector<int>& wave13_info,
	vector<int>& wave14_info,
	vector<int>& wave15_info,
	vector<int>& wave16_info,
	vector<int>& neumann1_info,
	vector<int>& neumann2_info,
	vector<int>& conv_info);
void s2i(const string& info, const string& sta, vector<int>& pec_info, vector<int>& pmc_info, vector<int>& abc_info,
	vector<int>& lump1_info,
	vector<int>& lump2_info,
	vector<int>& lump3_info,
	vector<int>& lump4_info,
	vector<int>& lump5_info,
	vector<int>& lump6_info,
	vector<int>& lump7_info,
	vector<int>& lump8_info,
	vector<int>& lump9_info,
	vector<int>& lump10_info,
	vector<int>& lump11_info,
	vector<int>& lump12_info,
	vector<int>& lump13_info,
	vector<int>& lump14_info,
	vector<int>& lump15_info,
	vector<int>& lump16_info,
	vector<int>& wave1_info,
	vector<int>& wave2_info,
	vector<int>& wave3_info,
	vector<int>& wave4_info,
	vector<int>& wave5_info,
	vector<int>& wave6_info,
	vector<int>& wave7_info,
	vector<int>& wave8_info,
	vector<int>& wave9_info,
	vector<int>& wave10_info,
	vector<int>& wave11_info,
	vector<int>& wave12_info,
	vector<int>& wave13_info,
	vector<int>& wave14_info,
	vector<int>& wave15_info,
	vector<int>& wave16_info,
	vector<int>& neumann1_info,
	vector<int>& neumann2_info,
	vector<int>& conv_info);


double minval(double a, double b, double c) {
	a = (a < b ? a : b);
	return (a < c ? a : c);
}

double maxval(double a, double b, double c) {
	a = (a > b ? a : b);
	return (a > c ? a : c);
}

double sumval(double a, double b, double c) {
	return a + b + c;
}

int Flg_opp(const Element& eli, const Element& elj) {

	//wether eli and elj are neighbor ?
	for (size_t i = 0; i != 4; ++i) {
		for (size_t j = 0; j != 4; ++j) {
			if (eli.node[i].ver == elj.node[j].ver)   return 1;//Only one point has the same coordinates
		}
	}
	return  0;
}


typedef struct Sum {
	size_t idx;
	size_t sum;
	//按math从大到小排序
	inline bool operator < (const Sum& x) const {
		return this->sum < x.sum;  //revise   return sum < x.sum;
	}

};
typedef struct FacePro {
	size_t idx;
	long long  pro;
	//按math从大到小排序
	inline bool operator < (const FacePro& x) const {
		return pro < x.pro;
	}

};
int Opp_Sct(Element*& el,int num_element_subdomain,int myid) {
	vector<Sum> opp_Sum((size_t)num_element_subdomain * 4);
#   pragma omp parallel for
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (int nn = 0; nn < 4; ++nn) {
			int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
			opp_Sum[(size_t)i * 4 + nn].idx = i * 4 + nn;
			opp_Sum[(size_t)i * 4 + nn].sum = el[i].node[f_n1].ver + el[i].node[f_n2].ver + el[i].node[f_n3].ver;
		}
	}
	sort(opp_Sum.begin(), opp_Sum.end());//The sum of the nodes on all faces is arranged from small to large
	for (auto beg = opp_Sum.begin(); beg != opp_Sum.end(); ) {
		size_t sum1 = (*beg).sum;
		

		vector<FacePro> v2;
		
		while (sum1 == (*beg).sum) {     
			int nn = ((*beg).idx) % 4; int ith = ((*beg).idx) / 4; FacePro temp1;
			int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
			long long pro_1 = (long long)el[ith].node[f_n1].ver * el[ith].node[f_n2].ver * el[ith].node[f_n3].ver;
			temp1.idx = (*beg).idx; temp1.pro = pro_1;
			v2.push_back(temp1);
			beg++;
		}
		int v2Size = v2.size();
//#       pragma omp parallel for num_threads((num_domain-1)*process_threads)
#   pragma omp parallel for
		for (int i = 0; i < v2Size; ++i) {
			for (int j = i + 1; j < v2Size; ++j) {
				if (v2[i].pro == v2[j].pro) {
					int ith = v2[i].idx / 4; int nn = v2[i].idx % 4;
					int op1 = v2[j].idx / 4; int ofn1 = v2[j].idx % 4;

					int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
					int f_n4 = face_node[ofn1][0] - 1, f_n5 = face_node[ofn1][1] - 1, f_n6 = face_node[ofn1][2] - 1;
					if ((el[ith].node[f_n1].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver) && \
						(el[ith].node[f_n2].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver) && \
						(el[ith].node[f_n3].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver)) {
						el[ith].face[nn].opp[0] = op1 + 1; el[ith].face[nn].opp[1] = ofn1 + 1;//opp[1] Stores the local number of an element
						el[op1].face[ofn1].opp[0] = ith + 1; el[op1].face[ofn1].opp[1] = nn + 1;//opp[0] Stores the global number of a common coplanar tetrahedron
					}
				}
			}
		}
	}

	int kabc1 = 0, kpec1 = 0, kabc2 = 0, kpmc1 = 0, kpmc2 = 0,kwave1=0, kwave2 = 0; int elAdj = 0;
	enum Boundary { PEC, PMC, ABC, LUMP1, LUMP2, LUMP3, LUMP4, LUMP5, LUMP6, LUMP7, LUMP8, LUMP9, LUMP10, LUMP11, LUMP12, LUMP13, LUMP14, LUMP15, LUMP16,
		CONV,WAVE1,WAVE2 , WAVE3, WAVE4, WAVE5, WAVE6, WAVE7, WAVE8, WAVE9, WAVE10, WAVE11, WAVE12, WAVE13, WAVE14, WAVE15, WAVE16,NEUMANN1, NEUMANN2};
	vector<int> pec_info, pmc_info, abc_info, lump1_info, lump2_info, lump3_info, lump4_info, lump5_info, lump6_info, lump7_info, lump8_info, 
		lump9_info, lump10_info, lump11_info, lump12_info, lump13_info, lump14_info, lump15_info, lump16_info,
		wave1_info, wave2_info, wave3_info, wave4_info, wave5_info, wave6_info, wave7_info, wave8_info,
		wave9_info, wave10_info, wave11_info, wave12_info, wave13_info, wave14_info, wave15_info, wave16_info,
		conv_info, neumann1_info, neumann2_info, bd_flag, bd_flag1, bd_flag2;
	vector <Boundary> bd_class(num_element_boundary);
	vector <Boundary> bd_class1(num_element_boundary);
	vector <Boundary> bd_class2(num_element_boundary);
	Read_BounaryInfo(pec_info, pmc_info, abc_info, lump1_info, lump2_info, lump3_info,lump4_info,lump5_info,lump6_info,lump7_info,lump8_info,lump9_info,lump10_info,
		lump11_info,lump12_info,lump13_info,lump14_info,lump15_info,lump16_info,
		wave1_info, wave2_info, wave3_info, wave4_info, wave5_info, wave6_info, wave7_info, wave8_info, wave9_info, wave10_info,
		wave11_info, wave12_info, wave13_info, wave14_info, wave15_info, wave16_info, neumann1_info, neumann2_info,
		conv_info);
	//for (auto beg1 = wave1_info.cbegin(); beg1 != wave1_info.cend(); ++beg1)cout << "wave1 = " << *beg1 << endl;
	//for (auto beg1 = wave2_info.cbegin(); beg1 != wave2_info.cend(); ++beg1)cout << "wave2 = " << *beg1 << endl;
	//for (auto beg1 = wave11_info.cbegin(); beg1 != wave11_info.cend(); ++beg1)cout << "wave11 = " << *beg1 << endl;

	for (int i = 0; i < num_element_boundary; ++i) {
		int s = boundaryInfo[i];
		for (auto beg1 = abc_info.cbegin(); beg1 != abc_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = ABC; ++kabc1; }//Number of faces belonging to ABC boundary
		for (auto beg1 = pmc_info.cbegin(); beg1 != pmc_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = PMC; ++kpmc2; }//Number of faces belonging to PMC boundary
		for (auto beg1 = pec_info.cbegin(); beg1 != pec_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = PEC; ++kpec1; }//Number of faces belonging to PEC boundary
		for (auto beg1 = lump1_info.cbegin(); beg1 != lump1_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = LUMP1; ++kpmc1; }//Number of faces belonging to LUMP boundary
		for (auto beg1 = wave1_info.cbegin(); beg1 != wave1_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE1; kwave1++;}
		for (auto beg1 = wave2_info.cbegin(); beg1 != wave2_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE2; }
		for (auto beg1 = wave3_info.cbegin(); beg1 != wave3_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE3; }//Number of faces belonging to CONV boundary
		for (auto beg1 = wave4_info.cbegin(); beg1 != wave4_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE4; }
		for (auto beg1 = wave5_info.cbegin(); beg1 != wave5_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE5; }
		for (auto beg1 = wave6_info.cbegin(); beg1 != wave6_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE6; }
		for (auto beg1 = wave7_info.cbegin(); beg1 != wave7_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE7; }
		for (auto beg1 = wave8_info.cbegin(); beg1 != wave8_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE8; }
		for (auto beg1 = wave9_info.cbegin(); beg1 != wave9_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE9; }
		for (auto beg1 = wave10_info.cbegin(); beg1 != wave10_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE10; }
		for (auto beg1 = wave11_info.cbegin(); beg1 != wave11_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE11; }
		for (auto beg1 = wave12_info.cbegin(); beg1 != wave12_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE12; }
		for (auto beg1 = wave13_info.cbegin(); beg1 != wave13_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE13; }
		for (auto beg1 = wave14_info.cbegin(); beg1 != wave14_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE14; }
		for (auto beg1 = wave15_info.cbegin(); beg1 != wave15_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE15; }
		for (auto beg1 = wave16_info.cbegin(); beg1 != wave16_info.cend(); ++beg1)   if (s == *beg1) { bd_flag.push_back(i); bd_class[i] = WAVE16; kwave2++;}
	}

	for (int i = 0; i < num_element_boundary; ++i) {
		int s = boundaryInfo[i];
		for (auto beg1 = neumann1_info.cbegin(); beg1 != neumann1_info.cend(); ++beg1)   if (s == *beg1) { bd_flag1.push_back(i); bd_class1[i] = NEUMANN1; }
		for (auto beg1 = neumann2_info.cbegin(); beg1 != neumann2_info.cend(); ++beg1)   if (s == *beg1) { bd_flag1.push_back(i); bd_class1[i] = NEUMANN2; }
	}

	for (int i = 0; i < num_element_boundary; ++i) {
		int s = boundaryInfo[i];
		for (auto beg1 = conv_info.cbegin(); beg1 != conv_info.cend(); ++beg1)   if (s == *beg1) { bd_flag2.push_back(i); bd_class2[i] = CONV; ++kpmc1; }
	}


	for (auto beg1 = bd_flag.cbegin(); beg1 != bd_flag.cend(); ++beg1) {
		int numb = (*beg1), verb_1 = Bondver[numb][0] - 1, verb_2 = Bondver[numb][1] - 1, verb_3 = Bondver[numb][2] - 1;
		double zxb = zb_boundary(verb_1, 0) + zb_boundary(verb_2, 0) + zb_boundary(verb_3, 0);
		double zyb = zb_boundary(verb_1, 1) + zb_boundary(verb_2, 1) + zb_boundary(verb_3, 1);
		double zzb = zb_boundary(verb_1, 2) + zb_boundary(verb_2, 2) + zb_boundary(verb_3, 2);
		Eigen::Vector3d vb(zxb, zyb, zzb);
//#      pragma omp parallel for num_threads((num_domain-1)*process_threads)
#   pragma omp parallel for
		for (int i = 0; i < num_element_subdomain; i++) {
			double yt = el[i].node[0].zb[1] + el[i].node[1].zb[1] + el[i].node[2].zb[1] + el[i].node[3].zb[1];
			for (int nn = 0; nn < 4; nn++) {
				int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
				double zx = el[i].node[f_n1].zb[0] + el[i].node[f_n2].zb[0] + el[i].node[f_n3].zb[0];
				double zy = el[i].node[f_n1].zb[1] + el[i].node[f_n2].zb[1] + el[i].node[f_n3].zb[1];
				double zz = el[i].node[f_n1].zb[2] + el[i].node[f_n2].zb[2] + el[i].node[f_n3].zb[2];
				Eigen::Vector3d v1(zx, zy, zz);
				if (abs(v1[0] - vb[0]) < 1e-6 && abs(v1[1] - vb[1]) < 1e-6&&abs(v1[2] - vb[2]) < 1e-6 && abs(v1.norm() - vb.norm()) < 1e-6) {
					switch (bd_class[numb]) {
					case ABC: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -7; break;
					case PEC: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -5; break;
					case PMC: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -6; break;
					case LUMP1:el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -100; break;
					case LUMP2:el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -101; break;
					case LUMP3: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -102; break;
					case LUMP4: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -103; break;
					case LUMP5: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -104; break;
					case LUMP6: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -105; break;
					case LUMP7: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -106; break;
					case LUMP8: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -107; break;
					case LUMP9: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -108; break;
					case LUMP10: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -109; break;
					case LUMP11: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -110; break;
					case LUMP12: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -111; break;
					case LUMP13: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -112; break;
					case LUMP14: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -113; break;
					case LUMP15: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -114; break;
					case LUMP16: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -115; break;
					case WAVE1:el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -200; break;
					case WAVE2:el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -201; break;
					case WAVE3: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -202; break;
					case WAVE4: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -203; break;
					case WAVE5: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -204; break;
					case WAVE6: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -205; break;
					case WAVE7: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -206; break;
					case WAVE8: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -207; break;
					case WAVE9: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -208; break;
					case WAVE10: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -209; break;
					case WAVE11: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -210; break;
					case WAVE12: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -211; break;
					case WAVE13: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -212; break;
					case WAVE14: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -213; break;
					case WAVE15: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -214; break;
					case WAVE16: el[i].face[nn].opp[0] = 0; el[i].face[nn].opp[1] = -215; break;
					case CONV: el[i].face[nn].opp_stru[0] = 0; el[i].face[nn].opp_stru[1] = -70; break;
					}
				}
			}
		}
	}

	for (auto beg1 = bd_flag1.cbegin(); beg1 != bd_flag1.cend(); ++beg1) {
		int numb = (*beg1), verb_1 = Bondver[numb][0] - 1, verb_2 = Bondver[numb][1] - 1, verb_3 = Bondver[numb][2] - 1;
		double zxb = zb_boundary(verb_1, 0) + zb_boundary(verb_2, 0) + zb_boundary(verb_3, 0);
		double zyb = zb_boundary(verb_1, 1) + zb_boundary(verb_2, 1) + zb_boundary(verb_3, 1);
		double zzb = zb_boundary(verb_1, 2) + zb_boundary(verb_2, 2) + zb_boundary(verb_3, 2);
		Eigen::Vector3d vb(zxb, zyb, zzb);
		//#      pragma omp parallel for num_threads((num_domain-1)*process_threads)
#   pragma omp parallel for
		for (int i = 0; i < num_element_subdomain; i++) {
			double yt = el[i].node[0].zb[1] + el[i].node[1].zb[1] + el[i].node[2].zb[1] + el[i].node[3].zb[1];
			for (int nn = 0; nn < 4; nn++) {
				int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
				double zx = el[i].node[f_n1].zb[0] + el[i].node[f_n2].zb[0] + el[i].node[f_n3].zb[0];
				double zy = el[i].node[f_n1].zb[1] + el[i].node[f_n2].zb[1] + el[i].node[f_n3].zb[1];
				double zz = el[i].node[f_n1].zb[2] + el[i].node[f_n2].zb[2] + el[i].node[f_n3].zb[2];
				Eigen::Vector3d v1(zx, zy, zz);
				if (abs(v1[0] - vb[0]) < 1e-6 && abs(v1[1] - vb[1]) < 1e-6 && abs(v1[2] - vb[2]) < 1e-6 && abs(v1.norm() - vb.norm()) < 1e-6) {
					switch (bd_class1[numb]) {

					case NEUMANN1: el[i].face[nn].opp_T_q[0] = 0; el[i].face[nn].opp_T_q[1] = -7; break;
					case NEUMANN2:el[i].face[nn].opp_T_q[0] = 0; el[i].face[nn].opp_T_q[1] = -8; break;
					}
				}
			}
		}
	}

	for (auto beg1 = bd_flag2.cbegin(); beg1 != bd_flag2.cend(); ++beg1) {
		int numb = (*beg1), verb_1 = Bondver[numb][0] - 1, verb_2 = Bondver[numb][1] - 1, verb_3 = Bondver[numb][2] - 1;
		double zxb = zb_boundary(verb_1, 0) + zb_boundary(verb_2, 0) + zb_boundary(verb_3, 0);
		double zyb = zb_boundary(verb_1, 1) + zb_boundary(verb_2, 1) + zb_boundary(verb_3, 1);
		double zzb = zb_boundary(verb_1, 2) + zb_boundary(verb_2, 2) + zb_boundary(verb_3, 2);
		Eigen::Vector3d vb(zxb, zyb, zzb);
		//#      pragma omp parallel for num_threads((num_domain-1)*process_threads)
#   pragma omp parallel for
		for (int i = 0; i < num_element_subdomain; i++) {
			double yt = el[i].node[0].zb[1] + el[i].node[1].zb[1] + el[i].node[2].zb[1] + el[i].node[3].zb[1];
			for (int nn = 0; nn < 4; nn++) {
				int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
				double zx = el[i].node[f_n1].zb[0] + el[i].node[f_n2].zb[0] + el[i].node[f_n3].zb[0];
				double zy = el[i].node[f_n1].zb[1] + el[i].node[f_n2].zb[1] + el[i].node[f_n3].zb[1];
				double zz = el[i].node[f_n1].zb[2] + el[i].node[f_n2].zb[2] + el[i].node[f_n3].zb[2];
				Eigen::Vector3d v1(zx, zy, zz);
				if (abs(v1[0] - vb[0]) < 1e-6 && abs(v1[1] - vb[1]) < 1e-6 && abs(v1[2] - vb[2]) < 1e-6 && abs(v1.norm() - vb.norm()) < 1e-6) {
					switch (bd_class2[numb]) {
					case CONV: el[i].face[nn].opp_stru[0] = 0; el[i].face[nn].opp_stru[1] = -70; break;
					}
				}
			}
		}
	}




	cout << "kwave1 = " << kwave1 << endl;
	cout << "kwave2 = " << kwave2 << endl;
	kpmc2_GB = kpmc2;
	kpmc1_GB = kpmc1;
	zb.resize(0, 0);
	delete[] Bondver;
	zb_boundary.setZero();
	int op0;

	for (int i = 0; i < num_element_subdomain; i++) {
		for (int nn = 0; nn < 4; nn++) {
			op0 = el[i].face[nn].opp[0]; int ofn = el[i].face[nn].opp[1];
			if (op0 == 0&& ofn==0) {
				el[i].face[nn].whether_boundary = 1;
				num_opp_element[myid]++;
				
			}
		}
	}

	Global_num_element_boundary = new int [num_opp_element[myid]];
	Local_num_face_boundary = new int[num_opp_element[myid]];
	int count = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		for (int nn = 0; nn < 4; nn++) {
			op0 = el[i].face[nn].opp[0]; int ofn = el[i].face[nn].opp[1];

			int op10=el[i].face[nn].opp_T_q[0]; int op11 = el[i].face[nn].opp_T_q[1];

			if (op0 == 0 && ofn == 0 && op10 == 0 && op11 == 0) { //gbhfbgnbg
				Global_num_element_boundary[count] = el[i].Global_num;
				Local_num_face_boundary[count] = nn+1;
				count++;
			}
		}
	}



	return 0;
}

//The following two functions are used to determine boundary information  

int Read_BounaryInfo(vector<int>&pec_info, vector<int> &pmc_info, vector<int> &abc_info,vector<int> & lump1_info, vector<int>& lump2_info, vector<int>& lump3_info,
	vector<int>& lump4_info,
	vector<int>& lump5_info,
	vector<int>& lump6_info,
	vector<int>& lump7_info,
	vector<int>& lump8_info,
	vector<int>& lump9_info,
	vector<int>& lump10_info,
	vector<int>& lump11_info,
	vector<int>& lump12_info,
	vector<int>& lump13_info,
	vector<int>& lump14_info,
	vector<int>& lump15_info,
	vector<int>& lump16_info,
	vector<int>& wave1_info,
	vector<int>& wave2_info,
	vector<int>& wave3_info,
	vector<int>& wave4_info,
	vector<int>& wave5_info,
	vector<int>& wave6_info,
	vector<int>& wave7_info,
	vector<int>& wave8_info,
	vector<int>& wave9_info,
	vector<int>& wave10_info,
	vector<int>& wave11_info,
	vector<int>& wave12_info,
	vector<int>& wave13_info,
	vector<int>& wave14_info,
	vector<int>& wave15_info,
	vector<int>& wave16_info,

	vector<int>& neumann1_info,
	vector<int>& neumann2_info,
	vector<int>& conv_info) 
{
	
		ifstream infile;
		string x;
		infile.open(BoundaryInfo);
		string sta;

		while (getline(infile, x)&&(!x.empty()))
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
							s2i(pmc, sta, pec_info, pmc_info, abc_info, lump1_info, lump2_info, lump3_info, lump4_info, lump5_info,
								lump6_info, lump7_info, lump8_info, lump9_info, lump10_info, lump11_info, lump12_info, lump13_info, lump14_info, lump15_info, lump16_info,
								wave1_info, wave2_info, wave3_info, wave4_info, wave5_info,
								wave6_info, wave7_info, wave8_info, wave9_info, wave10_info, wave11_info, wave12_info, wave13_info, wave14_info, wave15_info, wave16_info,
								neumann1_info, neumann2_info,
								conv_info);
							pmc.clear();
						}
					}
				}
				s2i(pmc, sta, pec_info, pmc_info, abc_info, lump1_info, lump2_info, lump3_info, lump4_info, lump5_info,
					lump6_info, lump7_info, lump8_info, lump9_info, lump10_info, lump11_info, lump12_info, lump13_info, lump14_info, lump15_info, lump16_info,
					wave1_info, wave2_info, wave3_info, wave4_info, wave5_info,
					wave6_info, wave7_info, wave8_info, wave9_info, wave10_info, wave11_info, wave12_info, wave13_info, wave14_info, wave15_info, wave16_info,
					neumann1_info, neumann2_info,
					conv_info);
				sta.clear();
			}
		}
		infile.close();
		return 0;
	}


bool str_compare(string str1, string str2) {
	int len = str2.length();
	for (int i = 0; i < len; i++) {
		if (str1[i] != str2[i]) {
			return false;
		}
	}
	return true;
}

void s2i(const string& info, const string& sta, vector<int>&pec_info, vector<int> &pmc_info, vector<int> &abc_info, 
	vector<int>& lump1_info, 
	vector<int>& lump2_info, 
	vector<int>& lump3_info,
	vector<int>& lump4_info,
	vector<int>& lump5_info,
	vector<int>& lump6_info,
	vector<int>& lump7_info,
	vector<int>& lump8_info,
	vector<int>& lump9_info,
	vector<int>& lump10_info,
	vector<int>& lump11_info,
	vector<int>& lump12_info,
	vector<int>& lump13_info,
	vector<int>& lump14_info,
	vector<int>& lump15_info,
	vector<int>& lump16_info,
	vector<int>& wave1_info, 
	vector<int>& wave2_info, 
	vector<int>& wave3_info,
	vector<int>& wave4_info,
	vector<int>& wave5_info,
	vector<int>& wave6_info,
	vector<int>& wave7_info,
	vector<int>& wave8_info,
	vector<int>& wave9_info,
	vector<int>& wave10_info,
	vector<int>& wave11_info,
	vector<int>& wave12_info,
	vector<int>& wave13_info,
	vector<int>& wave14_info,
	vector<int>& wave15_info,
	vector<int>& wave16_info,
	vector<int>& neumann1_info,
	vector<int>& neumann2_info,
	vector<int>& conv_info) 
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
			if (str_compare(sta, "PMC"))  pmc_info.push_back(j);
			else if (str_compare(sta, "PEC"))  pec_info.push_back(j);
			else if (str_compare(sta, "ABC"))  abc_info.push_back(j);
			else if (str_compare(sta, "LUMP10")) lump10_info.push_back(j);
			else if (str_compare(sta, "LUMP11")) lump11_info.push_back(j);
			else if (str_compare(sta, "LUMP12")) lump12_info.push_back(j);
			else if (str_compare(sta, "LUMP13")) lump13_info.push_back(j);
			else if (str_compare(sta, "LUMP14")) lump14_info.push_back(j);
			else if (str_compare(sta, "LUMP15")) lump15_info.push_back(j);
			else if (str_compare(sta, "LUMP16")) lump16_info.push_back(j);
			else if (str_compare(sta, "LUMP1"))  lump1_info.push_back(j);
			else if (str_compare(sta, "LUMP2"))  lump2_info.push_back(j);
			else if (str_compare(sta, "LUMP3"))  lump3_info.push_back(j);
			else if (str_compare(sta, "LUMP4"))  lump4_info.push_back(j);
			else if (str_compare(sta, "LUMP5"))  lump5_info.push_back(j);
			else if (str_compare(sta, "LUMP6"))  lump6_info.push_back(j);
			else if (str_compare(sta, "LUMP7"))  lump7_info.push_back(j);
			else if (str_compare(sta, "LUMP8"))  lump8_info.push_back(j);
			else if (str_compare(sta, "LUMP9"))  lump9_info.push_back(j);
			else if (str_compare(sta, "WAVE10")) wave10_info.push_back(j);
			else if (str_compare(sta, "WAVE11")) wave11_info.push_back(j);
			else if (str_compare(sta, "WAVE12")) wave12_info.push_back(j);
			else if (str_compare(sta, "WAVE13")) wave13_info.push_back(j);
			else if (str_compare(sta, "WAVE14")) wave14_info.push_back(j);
			else if (str_compare(sta, "WAVE15")) {
				wave15_info.push_back(j); //cout << "ok" << endl;
			}
			else if (str_compare(sta, "WAVE16")) wave16_info.push_back(j);
			else if (str_compare(sta, "WAVE1"))  wave1_info.push_back(j);
			else if (str_compare(sta, "WAVE2"))  wave2_info.push_back(j);
			else if (str_compare(sta, "WAVE3"))  wave3_info.push_back(j);
			else if (str_compare(sta, "WAVE4"))  wave4_info.push_back(j);
			else if (str_compare(sta, "WAVE5"))  wave5_info.push_back(j);
			else if (str_compare(sta, "WAVE6"))  wave6_info.push_back(j);
			else if (str_compare(sta, "WAVE7"))  wave7_info.push_back(j);
			else if (str_compare(sta, "WAVE8"))  wave8_info.push_back(j);
			else if (str_compare(sta, "WAVE9"))  wave9_info.push_back(j);
			else if (str_compare(sta, "NEUMANN1"))  neumann1_info.push_back(j);
			else if (str_compare(sta, "NEUMANN2"))  neumann2_info.push_back(j);
			else if (str_compare(sta, "CONV"))  conv_info.push_back(j);
		}
		temp.clear();
		sz = 0; idec_end = 0; i_dec = 0;

	}



int Opp_Sct1(Element*& el, int num_element_boundary_total, int myid, int* Global_num_element_boundary_total, int* Local_num_face_boundary_total, int* address_offset_face) {
	vector<Sum> opp_Sum((size_t)num_element_boundary_total);
#   pragma omp parallel for
	for (int i = 0; i < num_element_boundary_total; ++i) {
		int num_el = Global_num_element_boundary_total[i] - 1;
		int nn = Local_num_face_boundary_total[i] - 1;
		int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
		opp_Sum[(size_t)i].idx = i ;
		opp_Sum[(size_t)i].sum = el[num_el].node[f_n1].ver + el[num_el].node[f_n2].ver + el[num_el].node[f_n3].ver;
	}

	sort(opp_Sum.begin(), opp_Sum.end());//The sum of the nodes on all faces is arranged from small to large
	for (auto beg = opp_Sum.begin(); beg != opp_Sum.end(); ) {
		size_t sum1 = (*beg).sum;


		vector<FacePro> v2;

		while (sum1 == (*beg).sum) {
			int idx_temp = (*beg).idx;
			int nn = Local_num_face_boundary_total[idx_temp] - 1; int ith = Global_num_element_boundary_total[idx_temp] - 1;
			FacePro temp1;
			int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
			long long pro_1 = (long long)el[ith].node[f_n1].ver * el[ith].node[f_n2].ver * el[ith].node[f_n3].ver;
			temp1.idx = (*beg).idx; temp1.pro = pro_1;
			v2.push_back(temp1);
			beg++;
		}
		int v2Size = v2.size();
#   pragma omp parallel for
		for (int i = 0; i < v2Size; ++i) {
			for (int j = i + 1; j < v2Size; ++j) {
				if (v2[i].pro == v2[j].pro) {
					int num_eli = v2[i].idx; int num_elj = v2[j].idx;
					int ith = Global_num_element_boundary_total[num_eli] - 1; int nn = Local_num_face_boundary_total[num_eli] - 1;
					int op1 = Global_num_element_boundary_total[num_elj] - 1; int ofn1 = Local_num_face_boundary_total[num_elj] - 1;

					int f_n1 = face_node[nn][0] - 1, f_n2 = face_node[nn][1] - 1, f_n3 = face_node[nn][2] - 1;
					int f_n4 = face_node[ofn1][0] - 1, f_n5 = face_node[ofn1][1] - 1, f_n6 = face_node[ofn1][2] - 1;
					if ((el[ith].node[f_n1].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver) && \
						(el[ith].node[f_n2].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver) && \
						(el[ith].node[f_n3].ver == el[op1].node[f_n4].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n5].ver) || (el[ith].node[f_n1].ver == el[op1].node[f_n6].ver)) {
						el[ith].face[nn].opp[0] = op1 + 1; el[ith].face[nn].opp[1] = ofn1 + 1;//opp[1] Stores the local number of an element
						el[op1].face[ofn1].opp[0] = ith + 1; el[op1].face[ofn1].opp[1] = nn + 1;//opp[0] Stores the global number of a common coplanar tetrahedron
					}

					//el[ith].face[nn].opp[0] = op1 + 1; el[ith].face[nn].opp[1] = ofn1 + 1;//opp[1] Stores the local number of an element
					//el[op1].face[ofn1].opp[0] = ith + 1; el[op1].face[ofn1].opp[1] = nn + 1;//opp[0] Stores the global number of a common coplanar tetrahedron 
				}
			}
		}
	}



	return 0;
}








