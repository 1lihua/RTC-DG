#include "Material_Parm_T_q.h"

using namespace std;

int Material_Parm_T_q(Element* el_subdomain, int num_element_subdomain1) {
	int el, material, cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0, cnt5 = 0;

	for (el = 0; el < num_element_subdomain1; el++) {
		//material = element[el];
		material = el_subdomain[el].Material;

		switch (material)
		{
		case 1:
			el_subdomain[el].Kh[0] = 0.024;
			el_subdomain[el].Kh[1] = 0.024;
			el_subdomain[el].Kh[2] = 0.024;
			el_subdomain[el].Rho = 1.29;
			el_subdomain[el].Cp = 29.3;
			el_subdomain[el].E = 3.2e9;
			el_subdomain[el].nu = 0.35;
			el_subdomain[el].alpha = 7.0e-5;
			break;
		case 2:
			el_subdomain[el].Kh[0] = 0.3;
			el_subdomain[el].Kh[1] = 0.3;
			el_subdomain[el].Kh[2] = 0.3;
			el_subdomain[el].Rho = 2200;
			el_subdomain[el].Cp = 1050;
			el_subdomain[el].E = 3.2e9;
			el_subdomain[el].nu = 0.35;
			el_subdomain[el].alpha = 7.0e-5;
			break;
		case 3:
			el_subdomain[el].Kh[0] = 238;
			el_subdomain[el].Kh[1] = 238;
			el_subdomain[el].Kh[2] = 238;
			el_subdomain[el].Rho = 2700;
			el_subdomain[el].Cp = 900;
			el_subdomain[el].E = 110e9;
			el_subdomain[el].nu = 0.35;
			el_subdomain[el].alpha = 17e-6;
			break;

		default:
			cnt5++;
			el_subdomain[el].Kh[0] = 400;
			el_subdomain[el].Kh[1] = 400;
			el_subdomain[el].Kh[2] = 400;
			el_subdomain[el].Rho = 8700;
			el_subdomain[el].Cp = 385;
			el_subdomain[el].E = 110e9;
			el_subdomain[el].nu = 0.35;
			el_subdomain[el].alpha = 17e-6;
			break;

		}

	}
	//cout << "cnt5 = " << cnt5 << endl;

	return 0;
}