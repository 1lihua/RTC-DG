

#include"Epsil_Mu.h"

//int Epsil_Mu(Element* el,int num_element_domain) {
//	int material;
//	//int cnt1, cnt2, cnt3;
//	//cnt1 = 0; cnt2 = 0; cnt3 = 0;
//	double epsil_r = 1.0;
//	for (int ith = 0; ith < num_element_domain; ++ith) {
//		material = el[ith].Material;
//		if (material == 1) {
//			el[ith].epsil = epsilon0;
//			el[ith].eta = eta0;
//			el[ith].sigma = 0;
//		}
//		else if (material == 3) {
//			epsil_r = 1.0;
//			el[ith].epsil = epsilon0 * epsil_r;
//			el[ith].eta = eta0 / sqrt(epsil_r);
//			el[ith].sigma = 0.0;
//			//el[ith].sigma = 0;
//		}
//		else if (material == 2) {
//			epsil_r = 4.4;
//			el[ith].epsil = epsilon0 * epsil_r;
//			el[ith].eta = eta0 / sqrt(epsil_r);
//			el[ith].sigma = 0.00;
//		}
//		else {
//			el[ith].epsil = epsilon0;
//			el[ith].eta = eta0;
//			el[ith].sigma = 0;
//		}
//
//
//
//	}
//	return 0;
//}

int Read_BounaryInfo1(vector<int>& DIE1_info, vector<int>& DIE2_info, vector<int>& DIE3_info,
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
void s2i1(const string& info, const string& sta,
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

int Epsil_Mu(Element* el, int num_element_domain) {



	enum Boundary { DIE1, DIE2, DIE3, DIE4, DIE5, DIE6, DIE7, DIE8, DIE9, DIE10, DIE11, DIE12, DIE13, DIE14, DIE15, DIE16};
	vector<int>  DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info, DIE6_info, DIE7_info, DIE8_info,
		DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info, bd_flag;


	Read_BounaryInfo1(DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info, DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info,
		DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);
	//for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)
	//	cout << *beg1;
	//cout<< endl;
	//for (auto beg1 = DIE2_info.cbegin(); beg1 != DIE2_info.cend(); ++beg1)
	//	cout << *beg1 << endl;
	//cout << "lkmi" << endl;

	int material;
	//int cnt1, cnt2, cnt3;
	//cnt1 = 0; cnt2 = 0; cnt3 = 0;
	double epsil_r = 1.0;
	for (int ith = 0; ith < num_element_domain; ++ith) {
		el[ith].epsil = epsilon0;
		el[ith].eta = eta0;
		el[ith].sigma = 0;
	}

	for (int ith = 0; ith < num_element_domain; ++ith) {
		material = el[ith].Material;



		for (auto beg1 = DIE1_info.cbegin(); beg1 != DIE1_info.cend(); ++beg1)   if (material == *beg1) {
			el[ith].epsil = epsilon0* DIEeps[0];
			el[ith].eta = eta0 /sqrt( DIEeps[0]);
			el[ith].sigma = DIEsigma[0];
			el[ith].sigma_temp = DIEsigma[0];
			el[ith].Kh[0] = DIEk1[0];
			el[ith].Kh[1] = DIEk2[0];
			el[ith].Kh[2] = DIEk3[0];
			el[ith].Rho = DIERho[0];
			el[ith].Cp = DIESHC[0];
			el[ith].Q = DIEQ[0];
			el[ith].E = DIEE[0];
			el[ith].nu = DIEnu[0];
			el[ith].alpha = DIEalpha[0];
		}
		for (auto beg1 = DIE2_info.cbegin(); beg1 != DIE2_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE2];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE2]);
			el[ith].sigma = DIEsigma[DIE2];
			el[ith].sigma_temp = DIEsigma[DIE2];
			el[ith].Kh[0] = DIEk1[DIE2];
			el[ith].Kh[1] = DIEk2[DIE2];
			el[ith].Kh[2] = DIEk3[DIE2];
			el[ith].Rho = DIERho[DIE2];
			el[ith].Cp = DIESHC[DIE2];
			el[ith].Q = DIEQ[DIE2];
			el[ith].E = DIEE[DIE2];
			el[ith].nu = DIEnu[DIE2];
			el[ith].alpha = DIEalpha[DIE2];
		}
		for (auto beg1 = DIE3_info.cbegin(); beg1 != DIE3_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE3];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE3]);
			el[ith].sigma = DIEsigma[DIE3];
			el[ith].sigma_temp = DIEsigma[DIE3];
			el[ith].Kh[0] = DIEk1[DIE3];
			el[ith].Kh[1] = DIEk2[DIE3];
			el[ith].Kh[2] = DIEk3[DIE3];
			el[ith].Rho = DIERho[DIE3];
			el[ith].Cp = DIESHC[DIE3];
			el[ith].Q = DIEQ[DIE3];
			el[ith].E = DIEE[DIE3];
			el[ith].nu = DIEnu[DIE3];
			el[ith].alpha = DIEalpha[DIE3];
		}
		for (auto beg1 = DIE4_info.cbegin(); beg1 != DIE4_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE4];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE4]);
			el[ith].sigma = DIEsigma[DIE4];
			el[ith].sigma_temp = DIEsigma[DIE4];
			el[ith].Kh[0] = DIEk1[DIE4];
			el[ith].Kh[1] = DIEk2[DIE4];
			el[ith].Kh[2] = DIEk3[DIE4];
			el[ith].Rho = DIERho[DIE4];
			el[ith].Q = DIEQ[DIE4];
			el[ith].Cp = DIESHC[DIE4];
			el[ith].E = DIEE[DIE4];
			el[ith].nu = DIEnu[DIE4];
			el[ith].alpha = DIEalpha[DIE4];
		}
		for (auto beg1 = DIE5_info.cbegin(); beg1 != DIE5_info.cend(); ++beg1)   if (material == *beg1)
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE5];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE5]);
			el[ith].sigma = DIEsigma[DIE5];
			el[ith].sigma_temp = DIEsigma[DIE5];
			el[ith].Kh[0] = DIEk1[DIE5];
			el[ith].Kh[1] = DIEk2[DIE5];
			el[ith].Kh[2] = DIEk3[DIE5];
			el[ith].Rho = DIERho[DIE5];
			el[ith].Cp = DIESHC[DIE5];
			el[ith].Q = DIEQ[DIE5];
			el[ith].E = DIEE[DIE5];
			el[ith].nu = DIEnu[DIE5];
			el[ith].alpha = DIEalpha[DIE5];
		}
		for (auto beg1 = DIE6_info.cbegin(); beg1 != DIE6_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE6];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE6]);
			el[ith].sigma = DIEsigma[DIE6];
			el[ith].sigma_temp = DIEsigma[DIE6];
			el[ith].Kh[0] = DIEk1[DIE6];
			el[ith].Kh[1] = DIEk2[DIE6];
			el[ith].Kh[2] = DIEk3[DIE6];
			el[ith].Rho = DIERho[DIE6];
			el[ith].Cp = DIESHC[DIE6];
			el[ith].Q = DIEQ[DIE6];
			el[ith].E = DIEE[DIE6];
			el[ith].nu = DIEnu[DIE6];
			el[ith].alpha = DIEalpha[DIE6];
		}
		for (auto beg1 = DIE7_info.cbegin(); beg1 != DIE7_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE7];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE7]);
			el[ith].sigma = DIEsigma[DIE7];
			el[ith].sigma_temp = DIEsigma[DIE6];
			el[ith].Kh[0] = DIEk1[DIE7];
			el[ith].Kh[1] = DIEk2[DIE7];
			el[ith].Kh[2] = DIEk3[DIE7];
			el[ith].Rho = DIERho[DIE7];
			el[ith].Cp = DIESHC[DIE7];
			el[ith].Q = DIEQ[DIE7];
			el[ith].E = DIEE[DIE7];
			el[ith].nu = DIEnu[DIE7];
			el[ith].alpha = DIEalpha[DIE7];
		}
		for (auto beg1 = DIE8_info.cbegin(); beg1 != DIE8_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE8];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE8]);
			el[ith].sigma = DIEsigma[DIE8];
			el[ith].sigma_temp = DIEsigma[DIE8];
			el[ith].Kh[0] = DIEk1[DIE8];
			el[ith].Kh[1] = DIEk2[DIE8];
			el[ith].Kh[2] = DIEk3[DIE8];
			el[ith].Rho = DIERho[DIE8];
			el[ith].Cp = DIESHC[DIE8];
			el[ith].Q = DIEQ[DIE8];
			el[ith].E = DIEE[DIE8];
			el[ith].nu = DIEnu[DIE8];
			el[ith].alpha = DIEalpha[DIE8];
		}
		for (auto beg1 = DIE9_info.cbegin(); beg1 != DIE9_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE9];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE9]);
			el[ith].sigma = DIEsigma[DIE9];
			el[ith].sigma_temp = DIEsigma[DIE9];
			el[ith].Kh[0] = DIEk1[DIE9];
			el[ith].Kh[1] = DIEk2[DIE9];
			el[ith].Kh[2] = DIEk3[DIE9];
			el[ith].Rho = DIERho[DIE9];
			el[ith].Cp = DIESHC[DIE9];
			el[ith].Q = DIEQ[DIE9];
			el[ith].E = DIEE[DIE9];
			el[ith].nu = DIEnu[DIE9];
			el[ith].alpha = DIEalpha[DIE9];
		}
		for (auto beg1 = DIE10_info.cbegin(); beg1 != DIE10_info.cend(); ++beg1)   if (material == *beg1) 
		{
			el[ith].epsil = epsilon0 * DIEeps[DIE10];
			el[ith].eta = eta0 / sqrt(DIEeps[DIE10]);
			el[ith].sigma = DIEsigma[DIE10];
			el[ith].sigma_temp = DIEsigma[DIE10];
			el[ith].Kh[0] = DIEk1[DIE10];
			el[ith].Kh[1] = DIEk2[DIE10];
			el[ith].Kh[2] = DIEk3[DIE10];
			el[ith].Rho = DIERho[DIE10];
			el[ith].Cp = DIESHC[DIE10];
			el[ith].Q = DIEQ[DIE10];
			el[ith].E = DIEE[DIE10];
			el[ith].nu = DIEnu[DIE10];
			el[ith].alpha = DIEalpha[DIE10];
		}
		for (auto beg1 = DIE11_info.cbegin(); beg1 != DIE11_info.cend(); ++beg1)   if (material == *beg1);
		for (auto beg1 = DIE12_info.cbegin(); beg1 != DIE12_info.cend(); ++beg1)   if (material == *beg1);
		for (auto beg1 = DIE13_info.cbegin(); beg1 != DIE13_info.cend(); ++beg1)   if (material == *beg1);
		for (auto beg1 = DIE14_info.cbegin(); beg1 != DIE14_info.cend(); ++beg1)   if (material == *beg1);
		for (auto beg1 = DIE15_info.cbegin(); beg1 != DIE15_info.cend(); ++beg1)   if (material == *beg1);
		for (auto beg1 = DIE16_info.cbegin(); beg1 != DIE16_info.cend(); ++beg1)   if (material == *beg1);

		//if (material == 1) {
		//	el[ith].epsil = epsilon0;
		//	el[ith].eta = eta0;
		//	el[ith].sigma = 0;
		//}
		//else if (material == 3) {
		//	epsil_r = 1.0;
		//	el[ith].epsil = epsilon0 * epsil_r;
		//	el[ith].eta = eta0 / sqrt(epsil_r);
		//	el[ith].sigma = 0.0;
		//	//el[ith].sigma = 0;
		//}
		//else if (material == 2) {
		//	epsil_r = 4.4;
		//	el[ith].epsil = epsilon0 * epsil_r;
		//	el[ith].eta = eta0 / sqrt(epsil_r);
		//	el[ith].sigma = 0.00;
		//}
		//else {
		//	el[ith].epsil = epsilon0;
		//	el[ith].eta = eta0;
		//	el[ith].sigma = 0;
		//}



	}
	return 0;
}



int Read_BounaryInfo1( vector<int>& DIE1_info, vector<int>& DIE2_info, vector<int>& DIE3_info,
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
	infile.open("material_info.txt");
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
						s2i1(pmc, sta,  DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info,
							DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);
						pmc.clear();
					}
				}
			}
			s2i1(pmc, sta,  DIE1_info, DIE2_info, DIE3_info, DIE4_info, DIE5_info,
				DIE6_info, DIE7_info, DIE8_info, DIE9_info, DIE10_info, DIE11_info, DIE12_info, DIE13_info, DIE14_info, DIE15_info, DIE16_info);
			sta.clear();
		}
	}
	infile.close();
	return 0;
}


bool str_compare1(string str1, string str2) {
	int len = str2.length();
	for (int i = 0; i < len; i++) {
		if (str1[i] != str2[i]) {
			return false;
		}
	}
	return true;
}

void s2i1(const string& info, const string& sta,
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
		if (str_compare1(sta, "DIE10")) DIE10_info.push_back(j);
		else if (str_compare1(sta, "DIE11")) DIE11_info.push_back(j);
		else if (str_compare1(sta, "DIE12")) DIE12_info.push_back(j);
		else if (str_compare1(sta, "DIE13")) DIE13_info.push_back(j);
		else if (str_compare1(sta, "DIE14")) DIE14_info.push_back(j);
		else if (str_compare1(sta, "DIE15")) DIE15_info.push_back(j);
		else if (str_compare1(sta, "DIE16")) DIE16_info.push_back(j);
		else if (str_compare1(sta, "DIE1"))  DIE1_info.push_back(j);
		else if (str_compare1(sta, "DIE2"))  DIE2_info.push_back(j);
		else if (str_compare1(sta, "DIE3"))  DIE3_info.push_back(j);
		else if (str_compare1(sta, "DIE4"))  DIE4_info.push_back(j);
		else if (str_compare1(sta, "DIE5"))  DIE5_info.push_back(j);
		else if (str_compare1(sta, "DIE6"))  DIE6_info.push_back(j);
		else if (str_compare1(sta, "DIE7"))  DIE7_info.push_back(j);
		else if (str_compare1(sta, "DIE8"))  DIE8_info.push_back(j);
		else if (str_compare1(sta, "DIE9"))  DIE9_info.push_back(j);
	}
	temp.clear();
	sz = 0; idec_end = 0; i_dec = 0;

}
