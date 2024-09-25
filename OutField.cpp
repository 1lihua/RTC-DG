#include"Element.h"
#include"Global_Data.h"
#include<fstream>
#include<complex>
#include<vector>
#include<mpi.h>
#include<omp.h>

using namespace std;
double sign_judge(int edgeNo1, int edgeNo2);

int Out_field(Element* el,const int& nfreq) {
	cout.precision(16);

	long double zbx, zby, zbz, L1, L2;
	int  nn, edge_glb, node10, node11, node_glb, ii, kk , node12;
	int node1, node2, node3, edge_loc;
	complex<double> Ex0(0, 0), Ey0(0, 0), Ez0(0,0);
	complex<double> N0[3];
	complex<double>* Ex = new complex<double>[num_node];
	complex<double>* Ey = new complex<double>[num_node];
	complex<double>* Ez = new complex<double>[num_node];
	int edge_node_local[12][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 3, 1 }, { 2, 3 }, { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 3, 1 }, { 2, 3 } };
	for (int i = 0; i < num_node; i++) {
		Ex[i].real(0.0); Ex[i].imag(0.0);
		Ey[i].real(0.0); Ey[i].imag(0.0);
		Ez[i].real(0.0); Ez[i].imag(0.0);
	}
	double yt = 0.0;  int op = 0, ofn = 0;
	for (int ith = 0; ith != num_element; ++ith) {	
			yt = (el[ith].node[0].zb[zbY] + el[ith].node[1].zb[zbY] + el[ith].node[2].zb[zbY] + el[ith].node[3].zb[zbY]) / 4.0;
			for (int nn = 0; nn != 4; ++nn) {
				op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
				if (ofn == -6 ) {
					node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
					zbx = (el[ith].node[node1].zb[zbX] + el[ith].node[node2].zb[zbX] + el[ith].node[node3].zb[zbX]) / 3.0;
					zby = (el[ith].node[node1].zb[zbY] + el[ith].node[node2].zb[zbY] + el[ith].node[node3].zb[zbY]) / 3.0;
					zbz = (el[ith].node[node1].zb[zbZ] + el[ith].node[node2].zb[zbZ] + el[ith].node[node3].zb[zbZ]) / 3.0;
					for (int ii = 0; ii != 6; ++ii) {
						edge_loc = face_edge[nn][ii] - 1;
						edge_glb = el[ith].Eedge_GBNO[edge_loc];
						node11 = edge_node_local[edge_loc][0]; node12 = edge_node_local[edge_loc][1];

						if (edge_glb != 0) {
							//!Output field
							L1 = el[ith].a[node11] + el[ith].b[node11] * zbx + el[ith].c[node11] * zby + el[ith].d[node11] * zbz;

							L2 = el[ith].a[node12] + el[ith].b[node12] * zbx + el[ith].c[node12] * zby + el[ith].d[node12] * zbz;
							if (edge_loc <= 5) {
								//N0[0] = el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) - L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2)
								N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								//N1 = Eedge_GBNO(el, edge_loc + 12) * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) + L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2);
								N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							//Ez1 = Ez1 + Xsol1(edge_glb) * dcmplx(N1(3))
							Ez0 += Xsol[edge_glb - 1] * N0[2];
						}
					}
				}
			}
		
	}
	double aa; double  bb; double cc;
	
	aa = Hsubstrate; bb = Wmicrostrip * factor;  cc = kpmc1_GB;
	V11.push_back(((Ez0 * -Hsubstrate)) / cc);
	cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	S11 = (V11[nfreq] - factor) / factor;
	
	Ez0 = 0;
	for (int ith = 0; ith != num_element; ++ith) {
			yt = (el[ith].node[0].zb[zbY] + el[ith].node[1].zb[zbY] + el[ith].node[2].zb[zbY] + el[ith].node[3].zb[zbY]) / 4.0;
			for (int nn = 0; nn != 4; ++nn) {
				op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
				if ((ofn == -66)) {
					node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
					zbx = (el[ith].node[node1].zb[zbX] + el[ith].node[node2].zb[zbX] + el[ith].node[node3].zb[zbX]) / 3.0;
					zby = (el[ith].node[node1].zb[zbY] + el[ith].node[node2].zb[zbY] + el[ith].node[node3].zb[zbY]) / 3.0;
					zbz = (el[ith].node[node1].zb[zbZ] + el[ith].node[node2].zb[zbZ] + el[ith].node[node3].zb[zbZ]) / 3.0;
					for (int ii = 0; ii != 6; ++ii) {
						edge_loc = face_edge[nn][ii] - 1;
						edge_glb = el[ith].Eedge_GBNO[edge_loc];
						node11 = edge_node_local[edge_loc][0]; node12 = edge_node_local[edge_loc][1];
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * zbx + el[ith].c[node11] * zby + el[ith].d[node11] * zbz;
							L2 = el[ith].a[node12] + el[ith].b[node12] * zbx + el[ith].c[node12] * zby + el[ith].d[node12] * zbz;
							if (edge_loc <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							Ez0 += Xsol[edge_glb - 1] * N0[2];
						}
					}
				}
			}
		
	}

	
	aa = Hsubstrate; bb = Wmicrostrip * factor; cc = kpmc2_GB;
	V21.push_back(((Ez0 * -Hsubstrate)) / cc);
	cout << Freq / 1e9 << "GHz" << ',' << " V21[" << nfreq << "] = " << V21[nfreq] << endl;
	S21 = V21[nfreq] / factor;
	

	free(Xsol);
	return 0;
}

int Out_field_subdomain1(Element* el, const int& nfreq) {
	cout.precision(16);

	long double zbx, zby, zbz, L1, L2;
	int  nn, edge_glb, node11, ii, node12;
	int node1, node2, node3, edge_loc;
	complex<double> Ex0(0, 0), Ey0(0, 0), Ez0(0, 0);
	complex<double> E1x0(0, 0), E1y0(0, 0), E1z0(0, 0);//
	complex<double> E2x0(0, 0), E2y0(0, 0), E2z0(0, 0);//
	complex<double> Ecx(0, 0), Ecy(0, 0), Ecz(0, 0);//
	complex<double> N0[3];
	complex<double> N1[3];
	Ecz.real(factor / Hsubstrate);//

	int edge_node_local[12][2] = { { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 3, 1 }, { 2, 3 }, { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 2 }, { 3, 1 }, { 2, 3 } };
	int singal = 0;
	int op = 0, ofn = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
			if (ofn == -6) {
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zbx = (el[ith].node[node1].zb[zbX] + el[ith].node[node2].zb[zbX] + el[ith].node[node3].zb[zbX]) / 3.0;
				zby = (el[ith].node[node1].zb[zbY] + el[ith].node[node2].zb[zbY] + el[ith].node[node3].zb[zbY]) / 3.0;
				zbz = (el[ith].node[node1].zb[zbZ] + el[ith].node[node2].zb[zbZ] + el[ith].node[node3].zb[zbZ]) / 3.0;
				for (int ii = 0; ii != 6; ++ii) {

					edge_loc = face_edge[nn][ii] - 1;
					edge_glb = el[ith].Eedge_GBNO[edge_loc];
					node11 = edge_node_local[edge_loc][0]; node12 = edge_node_local[edge_loc][1];

					if (edge_glb != 0) {
						//!Output field
						L1 = el[ith].a[node11] + el[ith].b[node11] * zbx + el[ith].c[node11] * zby + el[ith].d[node11] * zbz;
						L2 = el[ith].a[node12] + el[ith].b[node12] * zbx + el[ith].c[node12] * zby + el[ith].d[node12] * zbz;
						if (edge_loc <= 5) {
							N1[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (el[ith].d[node12] - el[ith].d[node11])/6*2* el[ith].face[nn].Area/6/ el[ith].Ve);//

							complex<double> temp_X_con;
							temp_X_con.real(Xsol[edge_glb - 1].real());
							temp_X_con.imag(Xsol[edge_glb - 1].imag());
							E1z0 += N1[2] * Ecz * temp_X_con;

							N1[2].real( (el[ith].d[node12]* el[ith].d[node12] - el[ith].d[node12]*el[ith].d[node11]+ el[ith].d[node11]* el[ith].d[node11]) / 24 * 4 * el[ith].face[nn].Area);//
							double temp1 = Xsol[edge_glb - 1].imag() * Xsol[edge_glb - 1].imag() + Xsol[edge_glb - 1].real() * Xsol[edge_glb - 1].real();
							N1[2].real(N1[2].real() * temp1);
							E1z0 -= N1[2] * Ecz;
							E2z0 += N1[2] * Ecz;

							//N0[0] = el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) - L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2)
							N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
							N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
							N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
						}
						else {
							N1[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (el[ith].d[node12] + el[ith].d[node11]) / 6 * 2 * el[ith].face[nn].Area);//
							E1z0 += N1[2] * Ecz * Xsol[edge_glb - 1];

							N1[2].real((el[ith].d[node12] * el[ith].d[node12] + el[ith].d[node12] * el[ith].d[node11] + el[ith].d[node11] * el[ith].d[node11]) / 24 * 4 * el[ith].face[nn].Area);//
							double temp1 = Xsol[edge_glb - 1].imag() * Xsol[edge_glb - 1].imag() + Xsol[edge_glb - 1].real() * Xsol[edge_glb - 1].real();
							N1[2].real(N1[2].real() * temp1);
							E1z0 -= N1[2] * Ecz;
							E2z0 += N1[2] * Ecz;


							//N1 = Eedge_GBNO(el, edge_loc + 12) * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) + L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2);
							N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
							N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
							N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
						}
						//Ez1 = Ez1 + Xsol1(edge_glb) * dcmplx(N1(3))
						singal++;
						Ez0 += Xsol[edge_glb - 1] * N0[2];
					}
				}
			}
		}

	}
	double aa; double  bb; double cc;

	aa = Hsubstrate; bb = Wmicrostrip * factor;  cc = kpmc1_GB;


	V11.push_back(((Ez0 * -Hsubstrate)) / cc);


	S11 = (V11[nfreq] - factor) / factor;
	if (singal > 0) {
		ofstream ofs11("S11.csv", ios::app);
		cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << ',' << S11_1 << endl;
		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	}
	singal = 0;
	complex<double> S11_temp = E1z0 / E2z0;
	double S11_1_temp = 20 * log10(sqrt(S11_temp.real() * S11_temp.real() + S11_temp.imag() * S11_temp.imag()));
	//cout << (double)(Freq / 1e9) << ',' << S11_1_temp << ',' << endl;


	Ez0 = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
			if ((ofn == -66)) {
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zbx = (el[ith].node[node1].zb[zbX] + el[ith].node[node2].zb[zbX] + el[ith].node[node3].zb[zbX]) / 3.0;
				zby = (el[ith].node[node1].zb[zbY] + el[ith].node[node2].zb[zbY] + el[ith].node[node3].zb[zbY]) / 3.0;
				zbz = (el[ith].node[node1].zb[zbZ] + el[ith].node[node2].zb[zbZ] + el[ith].node[node3].zb[zbZ]) / 3.0;
				for (int ii = 0; ii != 6; ++ii) {
					edge_loc = face_edge[nn][ii] - 1;
					edge_glb = el[ith].Eedge_GBNO[edge_loc];
					node11 = edge_node_local[edge_loc][0]; node12 = edge_node_local[edge_loc][1];
					if (edge_glb != 0) {
						L1 = el[ith].a[node11] + el[ith].b[node11] * zbx + el[ith].c[node11] * zby + el[ith].d[node11] * zbz;
						L2 = el[ith].a[node12] + el[ith].b[node12] * zbx + el[ith].c[node12] * zby + el[ith].d[node12] * zbz;
						if (edge_loc <= 5) {
							N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
							N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
							N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
						}
						else {
							N0[0].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
							N0[1].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
							N0[2].real(el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[edge_loc] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
						}
						singal++;
						Ez0 += Xsol[edge_glb - 1] * N0[2];
					}
				}
			}
		}

	}


	aa = Hsubstrate; bb = Wmicrostrip * factor; cc = kpmc2_GB;
	V21.push_back(((Ez0 * -Hsubstrate)) / cc);
	S21 = V21[nfreq] / factor;
	if (singal > 0) {
		ofstream ofs21("S21.csv", ios::app);
		cout << Freq / 1e9 << "GHz" << ',' << " V21[" << nfreq << "] = " << V21[nfreq] << endl;
		//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << ',' << S21_1 << endl;
		ofs21 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
	}

	free(Xsol);
	return 0;
}


int Out_field_subdomain_LP(Element* el, const int& nfreq) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };

	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();


	complex<double>V_temp; V_temp.imag(0); V_temp.real(0);
	int singal = 0;

	double yc1;
	for (int ith = 0; ith < num_element_subdomain; ith++) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[ith].face[nn].opp[0]; ofn1 = el[ith].face[nn].opp[1];
			if ((ofn1 == -6)) {
				nVec[0].real(el[ith].face[nn].N_Vector[0]); nVec[0].imag(0.0);
				nVec[1].real(el[ith].face[nn].N_Vector[1]); nVec[1].imag(0.0);
				nVec[2].real(el[ith].face[nn].N_Vector[2]); nVec[2].imag(0.0);
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (ii = 0; ii < 6; ii++) {


					edgeii_loc = face_edge[nn][ii];

					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					edgeii_E = el[ith].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[ith].Eedge_GBNO[edgeii_loc + 12 - 1];


					if (edgeii_E != 0) {
						for (nGauss = 0; nGauss < nQuads; nGauss++) {
							rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
							rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
							rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
							//value of test basis function
							Li1[0] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].b[node_ii2 - 1];
							Li1[1] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].c[node_ii2 - 1];
							Li1[2] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].d[node_ii2 - 1];

							Li2[0] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].b[node_ii1 - 1];
							Li2[1] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].c[node_ii1 - 1];
							Li2[2] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].d[node_ii1 - 1];
							//Incident E has Z component;
							//EincX.real(0.0); EincX.imag(0.0);
							//EincY.real(0.0); EincY.imag(0.0);
							//EincZ.real(-factor / Hsubstrate); EincZ.imag(0.0);//E = V0/d

							//Einc[0].real(EincX.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[0].imag(EincX.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							//Einc[1].real(EincY.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[1].imag(EincY.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							//Einc[2].real(EincZ.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[2].imag(EincZ.imag() / (Z0 * Wmicrostrip / Hsubstrate));

							//ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
							//ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
							//ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];

							//ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
							//ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
							//ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];

							////incident H has x component
							////HincX.real(-factor / thickness);
							//HincX.real(0.0);
							//HincX.imag(0.0);
							//HincY.real(0.0); HincY.imag(0.0);
							//HincZ.real(0.0); HincZ.imag(0.0);
							//Hinc[0].real(HincX.real() / el[ith].eta); Hinc[0].imag(HincX.imag() / el[ith].eta);
							//Hinc[1].real(HincY.real() / el[ith].eta); Hinc[1].imag(HincY.imag() / el[ith].eta);
							//Hinc[2].real(HincZ.real() / el[ith].eta); Hinc[2].imag(HincZ.imag() / el[ith].eta);
							//ncrosH[0] = nVec[1] * Hinc[2] - nVec[2] * Hinc[1];
							//ncrosH[1] = nVec[2] * Hinc[0] - nVec[0] * Hinc[2];
							//ncrosH[2] = nVec[0] * Hinc[1] - nVec[1] * Hinc[0];
							singal++;
							V_temp.imag(V_temp.imag() + Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]) * el[ith].face[nn].Area * weight0[nGauss]);
							V_temp.real(V_temp.real() + Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]) * el[ith].face[nn].Area * weight0[nGauss]);

							//for (int ii0 = 0; ii0 < 3; ii0++) {
							//	singal++;
							//	//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
							//	V_temp.imag(V_temp.imag() + Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]);
							//	V_temp.real(V_temp.imag() + Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]);
							//}
						}
					}
				}
			}

		}
	}

	double aa; double  bb; double cc;
	aa = Hsubstrate; bb = Wmicrostrip * factor;  cc = kpmc1_GB;
	V11.push_back((V_temp / -Wmicrostrip));
	S11 = (V11[nfreq] - factor) / factor;

	if (singal > 0) {
		ofstream ofs11("S11.csv", ios::app);
		cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S11 = " << S11_1 << endl;
		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	}



	singal = 0;
	V_temp.imag(0); V_temp.real(0);
	for (int ith = 0; ith < num_element_subdomain; ith++) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[ith].face[nn].opp[0]; ofn1 = el[ith].face[nn].opp[1];
			if ((ofn1 == -66)) {
				nVec[0].real(el[ith].face[nn].N_Vector[0]); nVec[0].imag(0.0);
				nVec[1].real(el[ith].face[nn].N_Vector[1]); nVec[1].imag(0.0);
				nVec[2].real(el[ith].face[nn].N_Vector[2]); nVec[2].imag(0.0);
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (ii = 0; ii < 6; ii++) {


					edgeii_loc = face_edge[nn][ii];

					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					edgeii_E = el[ith].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[ith].Eedge_GBNO[edgeii_loc + 12 - 1];
					if (edgeii_E != 0) {
						for (nGauss = 0; nGauss < nQuads; nGauss++) {
							rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
							rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
							rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
							//value of test basis function
							Li1[0] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].b[node_ii2 - 1];
							Li1[1] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].c[node_ii2 - 1];
							Li1[2] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].d[node_ii2 - 1];

							Li2[0] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].b[node_ii1 - 1];
							Li2[1] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].c[node_ii1 - 1];
							Li2[2] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].d[node_ii1 - 1];
							//Incident E has Z component;
							//EincX.real(0.0); EincX.imag(0.0);
							//EincY.real(0.0); EincY.imag(0.0);
							//EincZ.real(-factor / Hsubstrate); EincZ.imag(0.0);//E = V0/d

							//Einc[0].real(EincX.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[0].imag(EincX.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							//Einc[1].real(EincY.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[1].imag(EincY.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							//Einc[2].real(EincZ.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[2].imag(EincZ.imag() / (Z0 * Wmicrostrip / Hsubstrate));

							//ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
							//ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
							//ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];

							//ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
							//ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
							//ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];

							////incident H has x component
							////HincX.real(-factor / thickness);
							//HincX.real(0.0);
							//HincX.imag(0.0);
							//HincY.real(0.0); HincY.imag(0.0);
							//HincZ.real(0.0); HincZ.imag(0.0);
							//Hinc[0].real(HincX.real() / el[ith].eta); Hinc[0].imag(HincX.imag() / el[ith].eta);
							//Hinc[1].real(HincY.real() / el[ith].eta); Hinc[1].imag(HincY.imag() / el[ith].eta);
							//Hinc[2].real(HincZ.real() / el[ith].eta); Hinc[2].imag(HincZ.imag() / el[ith].eta);
							//ncrosH[0] = nVec[1] * Hinc[2] - nVec[2] * Hinc[1];
							//ncrosH[1] = nVec[2] * Hinc[0] - nVec[0] * Hinc[2];
							//ncrosH[2] = nVec[0] * Hinc[1] - nVec[1] * Hinc[0];
							singal++;
							V_temp.imag(V_temp.imag() + Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[2] + Li2[2]) * el[ith].face[nn].Area * weight0[nGauss]);
							V_temp.real(V_temp.real() + Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[2] + Li2[2]) * el[ith].face[nn].Area * weight0[nGauss]);
							//for (int ii0 = 0; ii0 < 3; ii0++) {
							//	singal++;
							//	//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
							//	V_temp.imag(V_temp.imag() + Xsol[edgeii_E - 1].imag() * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * nVec[ii0].real() * el[ith].face[nn].Area * weight0[nGauss]);
							//	V_temp.real(V_temp.real() + Xsol[edgeii_E - 1].real() * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * nVec[ii0].real() * el[ith].face[nn].Area * weight0[nGauss]);
							//}
						}
					}
				}
			}

		}
	}


	aa = Hsubstrate; bb = Wmicrostrip * factor; cc = kpmc2_GB;
	V21.push_back(V_temp / -Wmicrostrip);
	S21 = V21[nfreq] / factor;
	if (singal > 0) {
		ofstream ofs21("S21.csv", ios::app);
		cout << Freq / 1e9 << "GHz" << ',' << " V21[" << nfreq << "] = " << V21[nfreq] << endl;
		//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S21 = " << S21_1 << endl;
		ofs21 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
	}



	//S21 = 1.0;
	return 0;
}


int Out_field_subdomain(Element* el, const int& nfreq) {
	cout.precision(16);
	int nQuads = 6;
	//double weight0[6] = { { 0.22338158 },{ 0.22338158 },{ 0.22338158 },{ 0.10995174 },{ 0.10995174 },{ 0.10995174 } };
	//double alpha0[6] = { { 0.10810301 },{ 0.44594849 },{ 0.44594849 },{ 0.81684757 },{ 0.09157621 },{ 0.09157621 } };
	//double beta0[6] = { { 0.44594849 },{ 0.10810301 },{ 0.44594849 },{ 0.09157621 },{ 0.81684757 },{ 0.09157621 } };
	//double gamma0[6] = { { 0.44594849 },{ 0.44594849 },{ 0.10810301 },{ 0.09157621 },{ 0.09157621 },{ 0.81684757 } };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -6)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;

				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				//zb1[zbX] = zb(ver(ith, node1) - 1, zbX); zb2[zbX] = zb(ver(ith, node2) - 1, zbX); zb3[zbX] = zb(ver(ith, node3) - 1, zbX);
				//zb1[zbY] = zb(ver(ith, node1) - 1, zbY); zb2[zbY] = zb(ver(ith, node2) - 1, zbY); zb3[zbY] = zb(ver(ith, node3) - 1, zbY);
				//zb1[zbZ] = zb(ver(ith, node1) - 1, zbZ); zb2[zbZ] = zb(ver(ith, node2) - 1, zbZ); zb3[zbZ] = zb(ver(ith, node3) - 1, zbZ);
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10*omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10*omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	S11 = Eout / Eincr;
	//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));

	if (signal > 0) {
		ofstream ofs11("S11.csv", ios::app);
		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S11 = " << S11_1 << endl;
		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	}



	signal = 0;

	//cout << S11_1 << endl;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -66)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10*omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10*omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];

						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	S21 = Eout / Eincr;
	if (signal > 0) {
		ofstream ofs21("S21.csv", ios::app);
		//cout << Freq / 1e9 << "GHz" << ',' << " V21[" << nfreq << "] = " << V21[nfreq] << endl;
		//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S21 = " << S21_1 << endl;
		ofs21 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
	}
	//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//cout << S21_1 << endl;

	return 0;
}

int Out_field_subdomain_new1(Element* el, const int& nfreq,int myid,int ET) {
	cout.precision(16);
	int nQuads = 6;
	//double weight0[6] = { { 0.22338158 },{ 0.22338158 },{ 0.22338158 },{ 0.10995174 },{ 0.10995174 },{ 0.10995174 } };
	//double alpha0[6] = { { 0.10810301 },{ 0.44594849 },{ 0.44594849 },{ 0.81684757 },{ 0.09157621 },{ 0.09157621 } };
	//double beta0[6] = { { 0.44594849 },{ 0.10810301 },{ 0.44594849 },{ 0.09157621 },{ 0.81684757 },{ 0.09157621 } };
	//double gamma0[6] = { { 0.44594849 },{ 0.44594849 },{ 0.10810301 },{ 0.09157621 },{ 0.09157621 },{ 0.81684757 } };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -6)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;

				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				//zb1[zbX] = zb(ver(ith, node1) - 1, zbX); zb2[zbX] = zb(ver(ith, node2) - 1, zbX); zb3[zbX] = zb(ver(ith, node3) - 1, zbX);
				//zb1[zbY] = zb(ver(ith, node1) - 1, zbY); zb2[zbY] = zb(ver(ith, node2) - 1, zbY); zb3[zbY] = zb(ver(ith, node3) - 1, zbY);
				//zb1[zbZ] = zb(ver(ith, node1) - 1, zbZ); zb2[zbZ] = zb(ver(ith, node2) - 1, zbZ); zb3[zbZ] = zb(ver(ith, node3) - 1, zbZ);
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}

	complex<double>* Eout_v = new complex<double>[num_domain];
	complex<double>* Eincr_v = new complex<double>[num_domain];

	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S11 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S11_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;

		}
		else {
			ofstream ofs11("S11.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
		}


	}


	//S11 = Eout / Eincr;
	//if (signal > 0) {
	//	ofstream ofs11("S11.csv", ios::app);
	//	//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//	double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//	//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//	cout << (double)(Freq / 1e9) << "S11 = " << S11_1 << endl;
	//	ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//}



	signal = 0;

	//cout << S11_1 << endl;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -66)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];

						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}


	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S21 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S21_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

		}
		else {
			ofstream ofs11("S21.csv", ios::app);
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
		}
		//ofstream ofs11("S21.csv", ios::app);
		////cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		////double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		//cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
		//ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

	}


	return 0;
}

int Out_field_subdomain_new233(Element* el, const int& nfreq, int myid, int ET) {
	cout.precision(16);
	int nQuads = 6;
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		//if (el[ith].Material == 1) {
		//	continue;
		//}
		double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
		//if (el[i].Material == 1) {
		if (zt < 0) {
			continue;
		}

		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -200)) {
				Kwave[2] = omega * sqrt(el[ith].epsil * mur);
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				xrod1 = 6.6e-3; yrod1 = 0;//revise



				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];

				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					vector<complex<double>> Einc{ EincX, EincY, EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}

	complex<double>* Eout_v = new complex<double>[num_domain];
	complex<double>* Eincr_v = new complex<double>[num_domain];

	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S11 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S11_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout <<  "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;

		}
		else {
			ofstream ofs11("S11.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout << "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
		}
	}


	//S11 = Eout / Eincr;
	//if (signal > 0) {
	//	ofstream ofs11("S11.csv", ios::app);
	//	//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//	double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//	//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//	cout << (double)(Freq / 1e9) << "S11 = " << S11_1 << endl;
	//	ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//}

	double W00 = 37.26e-3;
	double L00 = 27.9e-3;
	
	signal = 0;

	//cout << S11_1 << endl;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
		//if (el[i].Material == 1) {
		if (zt < 0) {
			continue;
		}
		Kwave[2] = omega * sqrt(el[ith].epsil * mur);
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -201)) {
				signal++;
				xrod1 = 6.6e-3 - 2 * L00; yrod1 = 0;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					EincX.real((cos(-Kwave[2] * rs[2]))* ((1 / r / log(rb_wave / ra_wave))* disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2]))* ((1 / r / log(rb_wave / ra_wave))* disx / r));
					EincY.real((cos(-Kwave[2] * rs[2]))* ((1 / r / log(rb_wave / ra_wave))* disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2]))* ((1 / r / log(rb_wave / ra_wave))* disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					vector<complex<double>> Einc{ EincX, EincY, EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];

						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}


	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S21 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S21_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

		}
		else {
			ofstream ofs11("S21.csv", ios::app);
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
		}
		//ofstream ofs11("S21.csv", ios::app);
		////cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		////double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		//cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
		//ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

	}






	signal = 0;

	//cout << S11_1 << endl;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		if (el[ith].Material == 1) {
			continue;
		}
		Kwave[2] = omega * sqrt(el[ith].epsil * mur);
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -215)) {
				signal++;
				xrod1 = 6.6e-3 - 6 * L00; yrod1 = W00 * 6;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					vector<complex<double>> Einc{ EincX, EincY, EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];

						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}


	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S21 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S161_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S161 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

		}
		else {
			ofstream ofs11("S161.csv", ios::app);
			double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S161 = " << S21_1 << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
		}
		//ofstream ofs11("S21.csv", ios::app);
		////cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		////double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		//cout << (double)(Freq / 1e9) << "  S21 = " << S21_1 << endl;
		//ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

	}













	return 0;
}

int Out_field_subdomain_new(Element* el, const int& nfreq, int myid, int ET) {
	cout.precision(16);
	int nQuads = 6;
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;

	//for (int ith = 0; ith != num_element_subdomain; ++ith) {
	//	//if (el[ith].Material == 1) {
	//	//	continue;
	//	//}
	//	double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
	//	if (el[ith].Material != wave_material) {
	//	//if (zt < 0) {
	//		continue;
	//	}
	//	for (int nn = 0; nn != 4; ++nn) {
	//		//ofn = opp[ith][nn + 4];
	//		ofn = el[ith].face[nn].opp[1];
	//		if ((ofn == -200)) {
	//			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
	//			signal++;
	//			node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
	//			//xrod1 = 6.6e-3; yrod1 = 0;//revise
	//			xrod1 = x_wave[0]; yrod1 = y_wave[0];
	//			zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
	//			zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
	//			zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
	//			for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
	//				rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
	//				rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
	//				rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
	//				disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
	//				r = sqrt(pow(disx, 2) + pow(disy, 2));
	//				xc = rs[0]; yc = rs[1]; zc = rs[2];
	//				EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
	//				EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
	//				EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
	//				EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
	//				EincZ.real(0.0);
	//				EincZ.imag(0.0);
	//				vector<complex<double>> Einc{ EincX, EincY, EincZ };
	//				for (int ii0 = 0; ii0 < 3; ii0++) {
	//					complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
	//					Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
	//					Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
	//				}
	//				Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
	//				for (int ii = 0; ii != 12; ++ii) {
	//					//edge_glb = Eedge_GBNO[ith][ii];
	//					edge_glb = el[ith].Eedge_GBNO[ii];
	//					node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
	//					if (edge_glb != 0) {
	//						L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
	//						L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;
	//						//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
	//						//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
	//						if (ii <= 5) {
	//							N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
	//							N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
	//							N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
	//						}
	//						else {
	//							N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
	//							N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
	//							N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
	//							//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
	//							//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
	//							//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
	//						}
	//						Ex += Xsol[edge_glb - 1] * N0[0];
	//						Ey += Xsol[edge_glb - 1] * N0[1];
	//						Ez += Xsol[edge_glb - 1] * N0[2];
	//					}
	//				}
	//				Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
	//				for (int ii0 = 0; ii0 < 3; ii0++) {
	//					complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
	//						- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
	//					Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
	//					Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
	//				}
	//			}
	//		}
	//	}
	//}
	//complex<double>* Eout_v = new complex<double>[num_domain];
	//complex<double>* Eincr_v = new complex<double>[num_domain];
	//MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	//if (myid == 0) {
	//	complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
	//	for (int ii0 = 0; ii0 < num_domain; ii0++) {
	//		Eincr_temp += Eincr_v[ii0];
	//		Eout_temp += Eout_v[ii0];
	//	}
	//	S11 = Eout_temp / Eincr_temp;
	//	if (ET) {
	//		ofstream ofs11("S11_double.csv", ios::app);
	//		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//		cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
	//		cout << "  Eout_temp = " << Eout_temp << endl;
	//		cout << "  Eincr_temp = " << Eincr_temp << endl;
	//		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//	}
	//	else {
	//		ofstream ofs11("S11.csv", ios::app);
	//		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//		cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
	//		cout << "  Eout_temp = " << Eout_temp << endl;
	//		cout << "  Eincr_temp = " << Eincr_temp << endl;
	//		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//	}
	//}



	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		//if (el[ith].Material == 1) {
		//	continue;
		//}
		double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
		//if (el[i].Material == 1) {
		if (zt < 0) {
			continue;
		}

		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -200)) {
				Kwave[2] = omega * sqrt(el[ith].epsil * mur);
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				//xrod1 = 6.6e-3; yrod1 = 0;//revise
				xrod1 = x_wave[0]; yrod1 = y_wave[0];


				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];

				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					vector<complex<double>> Einc{ powerful* EincX, powerful* EincY, powerful* EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}

	complex<double>* Eout_v = new complex<double>[num_domain];
	complex<double>* Eincr_v = new complex<double>[num_domain];

	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S11 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S11_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout << "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;

		}
		else {
			ofstream ofs11("S11.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout << "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
		}
	}


	for (int i = 1; i < num_wave; i++) {
		//cout << "x_wave[" << i << "] = " << x_wave[i] << endl;

		signal = 0;

		//cout << S11_1 << endl;
		Eincr = 0; Eout = 0;
		for (int ith = 0; ith != num_element_subdomain; ++ith) {
			double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
			//if (el[ith].Material != wave_material) {
			if (zt < wave_zzz) {
				continue;
			}
			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
			for (int nn = 0; nn != 4; ++nn) {
				//ofn = opp[ith][nn + 4];
				ofn = el[ith].face[nn].opp[1];
				if (ofn == (-200-i)) {
					signal++;
					xrod1 =  x_wave[i]; yrod1 = y_wave[i];
					node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
					zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
					zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
					zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
					for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
						rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
						rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
						rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
						xc = rs[0]; yc = rs[1]; zc = rs[2];
						disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
						r = sqrt(pow(disx, 2) + pow(disy, 2));
						EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
						EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
						EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
						EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
						EincZ.real(0.0);
						EincZ.imag(0.0);
						vector<complex<double>> Einc{ powerful* EincX, powerful* EincY, powerful* EincZ };
						for (int ii0 = 0; ii0 < 3; ii0++) {
							complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
							Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
							Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
						}
						Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
						for (int ii = 0; ii != 12; ++ii) {
							//edge_glb = Eedge_GBNO[ith][ii];
							edge_glb = el[ith].Eedge_GBNO[ii];

							node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
							if (edge_glb != 0) {
								//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
								//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
								L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
								L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

								if (ii <= 5) {
									N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
									N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
									N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

									//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
									//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
									//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
								}
								else {
									N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
									N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
									N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
									//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
									//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
									//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
								}
								Ex += Xsol[edge_glb - 1] * N0[0];
								Ey += Xsol[edge_glb - 1] * N0[1];
								Ez += Xsol[edge_glb - 1] * N0[2];
							}
						}
						Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
						for (int ii0 = 0; ii0 < 3; ii0++) {
							complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
								- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
							Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
							Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
						}
					}
				}
			}
		}


		MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


		if (myid == 0) {
			complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
			for (int ii0 = 0; ii0 < num_domain; ii0++) {
				Eincr_temp += Eincr_v[ii0];
				Eout_temp += Eout_v[ii0];

			}
			S21 = Eout_temp / Eincr_temp;
			if (ET) {
				string wenjianming = "S" + to_string(i + 1) + "1_double.csv";

				ofstream ofs11(wenjianming, ios::app);
				//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
				double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				cout << (double)(Freq / 1e9) << "  S"<<i+1<<"1 = " << S21_1 << endl;
				ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

			}
			else {
				string wenjianming = "S" + to_string(i + 1) + "1.csv";
				ofstream ofs11(wenjianming, ios::app);
				double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				cout << (double)(Freq / 1e9) << "  S" << i + 1 << "1 = " << S21_1 << endl;
				ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
			}

		}

	}


	return 0;
}


int Out_field_subdomain_new_test(Element* el, const int& nfreq, int myid, int ET) {
	cout.precision(16);
	int nQuads = 6;
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;

	//for (int ith = 0; ith != num_element_subdomain; ++ith) {
	//	//if (el[ith].Material == 1) {
	//	//	continue;
	//	//}
	//	double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
	//	if (el[ith].Material != wave_material) {
	//	//if (zt < 0) {
	//		continue;
	//	}
	//	for (int nn = 0; nn != 4; ++nn) {
	//		//ofn = opp[ith][nn + 4];
	//		ofn = el[ith].face[nn].opp[1];
	//		if ((ofn == -200)) {
	//			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
	//			signal++;
	//			node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
	//			//xrod1 = 6.6e-3; yrod1 = 0;//revise
	//			xrod1 = x_wave[0]; yrod1 = y_wave[0];
	//			zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
	//			zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
	//			zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
	//			for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
	//				rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
	//				rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
	//				rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
	//				disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
	//				r = sqrt(pow(disx, 2) + pow(disy, 2));
	//				xc = rs[0]; yc = rs[1]; zc = rs[2];
	//				EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
	//				EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disx / r));
	//				EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
	//				EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / rb_wave)) * disy / r));
	//				EincZ.real(0.0);
	//				EincZ.imag(0.0);
	//				vector<complex<double>> Einc{ EincX, EincY, EincZ };
	//				for (int ii0 = 0; ii0 < 3; ii0++) {
	//					complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
	//					Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
	//					Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
	//				}
	//				Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
	//				for (int ii = 0; ii != 12; ++ii) {
	//					//edge_glb = Eedge_GBNO[ith][ii];
	//					edge_glb = el[ith].Eedge_GBNO[ii];
	//					node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
	//					if (edge_glb != 0) {
	//						L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
	//						L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;
	//						//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
	//						//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
	//						if (ii <= 5) {
	//							N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
	//							N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
	//							N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
	//						}
	//						else {
	//							N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
	//							N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
	//							N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
	//							//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
	//							//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
	//							//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
	//						}
	//						Ex += Xsol[edge_glb - 1] * N0[0];
	//						Ey += Xsol[edge_glb - 1] * N0[1];
	//						Ez += Xsol[edge_glb - 1] * N0[2];
	//					}
	//				}
	//				Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
	//				for (int ii0 = 0; ii0 < 3; ii0++) {
	//					complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
	//						- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
	//					Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
	//					Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
	//				}
	//			}
	//		}
	//	}
	//}
	//complex<double>* Eout_v = new complex<double>[num_domain];
	//complex<double>* Eincr_v = new complex<double>[num_domain];
	//MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	//if (myid == 0) {
	//	complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
	//	for (int ii0 = 0; ii0 < num_domain; ii0++) {
	//		Eincr_temp += Eincr_v[ii0];
	//		Eout_temp += Eout_v[ii0];
	//	}
	//	S11 = Eout_temp / Eincr_temp;
	//	if (ET) {
	//		ofstream ofs11("S11_double.csv", ios::app);
	//		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//		cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
	//		cout << "  Eout_temp = " << Eout_temp << endl;
	//		cout << "  Eincr_temp = " << Eincr_temp << endl;
	//		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//	}
	//	else {
	//		ofstream ofs11("S11.csv", ios::app);
	//		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
	//		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
	//		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//		cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
	//		cout << "  Eout_temp = " << Eout_temp << endl;
	//		cout << "  Eincr_temp = " << Eincr_temp << endl;
	//		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	//	}
	//}



	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		//if (el[ith].Material == 1) {
		//	continue;
		//}
		double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
		//if (el[i].Material == 1) {
		if (zt < 0) {
			continue;
		}

		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -200)) {
				Kwave[2] = omega * sqrt(el[ith].epsil * mur);
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				//xrod1 = 6.6e-3; yrod1 = 0;//revise
				xrod1 = x_wave[0]; yrod1 = y_wave[0];


				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];

				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					vector<complex<double>> Einc{ Amplitude_wave[0] * EincX, Amplitude_wave[0] * EincY, Amplitude_wave[0] * EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}

	complex<double>* Eout_v = new complex<double>[num_domain];
	complex<double>* Eincr_v = new complex<double>[num_domain];

	MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


	if (myid == 0) {
		complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
		for (int ii0 = 0; ii0 < num_domain; ii0++) {
			Eincr_temp += Eincr_v[ii0];
			Eout_temp += Eout_v[ii0];

		}
		S11 = Eout_temp / Eincr_temp;
		if (ET) {
			ofstream ofs11("S11_double.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout << "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;

		}
		else {
			ofstream ofs11("S11.csv", ios::app);
			//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
			double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
			//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
			cout << (double)(Freq / 1e9) << "  S11 = " << S11_1 << endl;
			cout << "  Eout_temp = " << Eout_temp << endl;
			cout << "  Eincr_temp = " << Eincr_temp << endl;
			ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
		}
	}


	for (int i = 1; i < num_wave; i++) {
		//cout << "x_wave[" << i << "] = " << x_wave[i] << endl;

		signal = 0;

		//cout << S11_1 << endl;
		Eincr = 0; Eout = 0;
		for (int ith = 0; ith != num_element_subdomain; ++ith) {
			double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
			//if (el[ith].Material != wave_material) {
			if (zt < wave_zzz) {
				continue;
			}
			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
			for (int nn = 0; nn != 4; ++nn) {
				//ofn = opp[ith][nn + 4];
				ofn = el[ith].face[nn].opp[1];
				if (ofn == (-200 - i)) {
					signal++;
					xrod1 = x_wave[i]; yrod1 = y_wave[i];
					node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
					zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
					zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
					zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
					for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
						rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
						rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
						rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
						xc = rs[0]; yc = rs[1]; zc = rs[2];
						disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
						r = sqrt(pow(disx, 2) + pow(disy, 2));
						EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
						EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
						EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
						EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
						EincZ.real(0.0);
						EincZ.imag(0.0);
						vector<complex<double>> Einc{ Amplitude_wave[0] * EincX, Amplitude_wave[0] * EincY, Amplitude_wave[0] * EincZ };
						for (int ii0 = 0; ii0 < 3; ii0++) {
							complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
							Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
							Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
						}

						Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
						for (int ii = 0; ii != 12; ++ii) {
							//edge_glb = Eedge_GBNO[ith][ii];
							edge_glb = el[ith].Eedge_GBNO[ii];

							node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
							if (edge_glb != 0) {
								//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
								//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
								L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
								L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

								if (ii <= 5) {
									N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
									N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
									N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

									//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
									//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
									//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
								}
								else {
									N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
									N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
									N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
									//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
									//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
									//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
								}
								Ex += Xsol[edge_glb - 1] * N0[0];
								Ey += Xsol[edge_glb - 1] * N0[1];
								Ez += Xsol[edge_glb - 1] * N0[2];
							}
						}
						Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
						for (int ii0 = 0; ii0 < 3; ii0++) {
							//complex<double> x((Etot[ii0].real() - Einc[ii0].real())* Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							//	- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
							//Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
							//Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
							complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
								- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
							Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
							Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
						}
					}
				}
			}
		}


		MPI_Gather(&Eout, 1, MPI_DOUBLE_COMPLEX, Eout_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Gather(&Eincr, 1, MPI_DOUBLE_COMPLEX, Eincr_v, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);


		if (myid == 0) {
			complex<double> Eincr_temp(0, 0); complex<double> Eout_temp(0, 0);
			for (int ii0 = 0; ii0 < num_domain; ii0++) {
				Eincr_temp += Eincr_v[ii0];
				Eout_temp += Eout_v[ii0];

			}
			S21 = Eout_temp / Eincr_temp;
			if (ET) {
				string wenjianming = "S" + to_string(i + 1) + "1_double.csv";

				ofstream ofs11(wenjianming, ios::app);
				//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
				double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				cout << (double)(Freq / 1e9) << "  S" << i + 1 << "1 = " << S21_1 << endl;
				ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;

			}
			else {
				string wenjianming = "S" + to_string(i + 1) + "1.csv";
				ofstream ofs11(wenjianming, ios::app);
				double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
				cout << (double)(Freq / 1e9) << "  S" << i + 1 << "1 = " << S21_1 << endl;
				ofs11 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
			}

		}

	}


	return 0;
}






int Out_field_subdomain_double(Element* el, const int& nfreq) {
	cout.precision(16);
	int nQuads = 6;
	//double weight0[6] = { { 0.22338158 },{ 0.22338158 },{ 0.22338158 },{ 0.10995174 },{ 0.10995174 },{ 0.10995174 } };
	//double alpha0[6] = { { 0.10810301 },{ 0.44594849 },{ 0.44594849 },{ 0.81684757 },{ 0.09157621 },{ 0.09157621 } };
	//double beta0[6] = { { 0.44594849 },{ 0.10810301 },{ 0.44594849 },{ 0.09157621 },{ 0.81684757 },{ 0.09157621 } };
	//double gamma0[6] = { { 0.44594849 },{ 0.44594849 },{ 0.10810301 },{ 0.09157621 },{ 0.09157621 },{ 0.81684757 } };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	double Kwave[3], L1, L2;
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	complex<double> N0[3];
	vector<complex<double>> Etot(3, (0, 0));
	vector<double> rs(3, 0);
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	int nn, edge_glb, node10, node11, node_glb, ii, kk, node12, node1, node2, node3;
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;// zbX = 0, zbY = 1, zbZ = 2;
	int signal = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -6)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;

				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				//zb1[zbX] = zb(ver(ith, node1) - 1, zbX); zb2[zbX] = zb(ver(ith, node2) - 1, zbX); zb3[zbX] = zb(ver(ith, node3) - 1, zbX);
				//zb1[zbY] = zb(ver(ith, node1) - 1, zbY); zb2[zbY] = zb(ver(ith, node2) - 1, zbY); zb3[zbY] = zb(ver(ith, node3) - 1, zbY);
				//zb1[zbZ] = zb(ver(ith, node1) - 1, zbZ); zb2[zbZ] = zb(ver(ith, node2) - 1, zbZ); zb3[zbZ] = zb(ver(ith, node3) - 1, zbZ);
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));

					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	S11 = Eout / Eincr;
	//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));

	if (signal > 0) {
		ofstream ofs11("S11_double.csv", ios::app);
		//cout << " V11[" << nfreq << "] = " << V11[nfreq] << endl;
		double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S11 = " << S11_1 << endl;
		ofs11 << (double)(Freq / 1e9) << ',' << S11_1 << endl;
	}



	signal = 0;

	//cout << S11_1 << endl;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element_subdomain; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			//ofn = opp[ith][nn + 4];
			ofn = el[ith].face[nn].opp[1];
			if ((ofn == -66)) {
				signal++;
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[face_node[nn][0] - 1].zb[0]; zb2[0] = el[ith].node[face_node[nn][1] - 1].zb[0]; zb3[0] = el[ith].node[face_node[nn][2] - 1].zb[0];
				zb1[1] = el[ith].node[face_node[nn][0] - 1].zb[1]; zb2[1] = el[ith].node[face_node[nn][1] - 1].zb[1]; zb3[1] = el[ith].node[face_node[nn][2] - 1].zb[1];
				zb1[2] = el[ith].node[face_node[nn][0] - 1].zb[2]; zb2[2] = el[ith].node[face_node[nn][1] - 1].zb[2]; zb3[2] = el[ith].node[face_node[nn][2] - 1].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					EincX.real(0.0); EincX.imag(0.0);
					EincY.real(0.0); EincY.imag(0.0);
					//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[ith] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
					EincZ.real(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					EincZ.imag(A10 * omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}
					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						//edge_glb = Eedge_GBNO[ith][ii];
						edge_glb = el[ith].Eedge_GBNO[ii];

						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;
						if (edge_glb != 0) {
							//L1 = a(ith, node11) + b(ith, node11) * xc + c(ith, node11) * yc + d(ith, node11) * zc;
							//L2 = a(ith, node12) + b(ith, node12) * xc + c(ith, node12) * yc + d(ith, node12) * zc;
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);

								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) - L2 * b(ith, node11)) * l(ith, ii) / (36.0 * pow(Ve(ith), 2))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) - L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) - L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							else {
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
								//N0[0].real(Eedge_GBNO[ith][ii + 12] * (L1 * b(ith, node12) + L2 * b(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[0].imag(0.0);
								//N0[1].real(Eedge_GBNO[ith][ii + 12] * (L1 * c(ith, node12) + L2 * c(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[1].imag(0.0);
								//N0[2].real(Eedge_GBNO[ith][ii + 12] * (L1 * d(ith, node12) + L2 * d(ith, node11)) * l(ith, ii) / (36.0 * Ve(ith) * Ve(ith))); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	S21 = Eout / Eincr;
	if (signal > 0) {
		ofstream ofs21("S21_double.csv", ios::app);
		//cout << Freq / 1e9 << "GHz" << ',' << " V21[" << nfreq << "] = " << V21[nfreq] << endl;
		//double S11_1 = 20 * log10(sqrt(S11.real() * S11.real() + S11.imag() * S11.imag()));
		double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
		cout << (double)(Freq / 1e9) << "S21 = " << S21_1 << endl;
		ofs21 << (double)(Freq / 1e9) << ',' << S21_1 << endl;
	}
	//double S21_1 = 20 * log10(sqrt(S21.real() * S21.real() + S21.imag() * S21.imag()));
	//cout << S21_1 << endl;

	return 0;
}


int Out_field_sub(Element* el, const int& nfreq) {
	cout.precision(16);
	int nQuads = 6;
	double weight0[6] = { { 0.22338158 },{ 0.22338158 },{ 0.22338158 },{ 0.10995174 },{ 0.10995174 },{ 0.10995174 } };
	double alpha0[6] = { { 0.10810301 },{ 0.44594849 },{ 0.44594849 },{ 0.81684757 },{ 0.09157621 },{ 0.09157621 } };
	double beta0[6] = { { 0.44594849 },{ 0.10810301 },{ 0.44594849 },{ 0.09157621 },{ 0.81684757 },{ 0.09157621 } };
	double gamma0[6] = { { 0.44594849 },{ 0.44594849 },{ 0.10810301 },{ 0.09157621 },{ 0.09157621 },{ 0.81684757 } };
	vector<double> rs(3, 0);
	double Kwave[3];
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	cout << "Kwave = " << Kwave[2] << endl;
	complex<double> Eincr(0, 0), Eout1(0, 0), Eout2(0, 0), Eout3(0, 0), Eincr1(0, 0), Eout(0, 0);
	complex<double> EincX(0, 0), EincY(0, 0), EincZ(0, 0), HincX(0, 0), HincY(0, 0), HincZ(0, 0);
	complex<double> Ex(0, 0), Ey(0, 0), Ez(0, 0);
	vector<complex<double>> Etot(3, (0, 0));
	int  nn, edge_glb, node10, node11, node_glb, ii, kk, node12;
	int node1, node2, node3;
	double L1, L2;
	complex<double> N0[3];
	vector<double> zb1(3, 0), zb2(3, 0), zb3(3, 0);
	double xc, yc, zc;
	double disx, disy, r;
	int op = 0, ofn = 0;
	for (int ith = 0; ith != num_element; ++ith) {
		Kwave[2] = omega * sqrt(el[ith].epsil * mur);
		for (int nn = 0; nn != 4; ++nn) {
			op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
			if (ofn == -6) {
				if (ofn == -6) { xrod1 = 6.6e-3; yrod1 = 0; }
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[node1].zb[0]; zb2[0] = el[ith].node[node2].zb[0]; zb3[0] = el[ith].node[node3].zb[0];
				zb1[1] = el[ith].node[node1].zb[1]; zb2[1] = el[ith].node[node2].zb[1]; zb3[1] = el[ith].node[node3].zb[1];
				zb1[2] = el[ith].node[node1].zb[2]; zb2[2] = el[ith].node[node2].zb[2]; zb3[2] = el[ith].node[node3].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					//Incident E has R component;
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					//cout << "EincX = " << EincX << '\t' << "EincY = " << EincY << '\t' << "EincZ = " << EincZ << endl;
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}

					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;

						if (edge_glb != 0) {
							//!Output field
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								//N0[0] = el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) - L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2)
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								//N1 = Eedge_GBNO(el, ii + 12) * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) + L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, ii) / (36.d0 * Ve(el) * *2);
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real() - Einc[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag() - Einc[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);

					}
				}
			}
		}
	}
	double W00 = 37.26e-3;
	double L00 = 27.9e-3;
	S11 = Eout / Eincr;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
			if (ofn == -66) {
				if (ofn == -66) { xrod1 = 6.6e-3 - 2 * L00; yrod1 = 0; }
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[node1].zb[0]; zb2[0] = el[ith].node[node2].zb[0]; zb3[0] = el[ith].node[node3].zb[0];
				zb1[1] = el[ith].node[node1].zb[1]; zb2[1] = el[ith].node[node2].zb[1]; zb3[1] = el[ith].node[node3].zb[1];
				zb1[2] = el[ith].node[node1].zb[2]; zb2[2] = el[ith].node[node2].zb[2]; zb3[2] = el[ith].node[node3].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					//Incident E has R component;
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					//cout << "EincX = " << EincX << '\t' << "EincY = " << EincY << '\t' << "EincZ = " << EincZ << endl;
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}

					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;

						if (edge_glb != 0) {
							//!Output field
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								//N0[0] = el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) - L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2)
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								//N1 = Eedge_GBNO(el, ii + 12) * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) + L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, ii) / (36.d0 * Ve(el) * *2);
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	S21 = Eout / Eincr;
	Eincr = 0; Eout = 0;
	for (int ith = 0; ith != num_element; ++ith) {
		for (int nn = 0; nn != 4; ++nn) {
			op = el[ith].face[nn].opp[0]; ofn = el[ith].face[nn].opp[1];
			if (ofn == -616) {
				if (ofn == -616) { xrod1 = 6.6e-3 - 6 * L00; yrod1 = W00 * 6; }
				node1 = face_node[nn][0] - 1; node2 = face_node[nn][1] - 1; node3 = face_node[nn][2] - 1;
				zb1[0] = el[ith].node[node1].zb[0]; zb2[0] = el[ith].node[node2].zb[0]; zb3[0] = el[ith].node[node3].zb[0];
				zb1[1] = el[ith].node[node1].zb[1]; zb2[1] = el[ith].node[node2].zb[1]; zb3[1] = el[ith].node[node3].zb[1];
				zb1[2] = el[ith].node[node1].zb[2]; zb2[2] = el[ith].node[node2].zb[2]; zb3[2] = el[ith].node[node3].zb[2];
				for (int nGauss = 0; nGauss < nQuads; ++nGauss) {
					rs[0] = alpha0[nGauss] * zb1[0] + beta0[nGauss] * zb2[0] + gamma0[nGauss] * zb3[0];
					rs[1] = alpha0[nGauss] * zb1[1] + beta0[nGauss] * zb2[1] + gamma0[nGauss] * zb3[1];
					rs[2] = alpha0[nGauss] * zb1[2] + beta0[nGauss] * zb2[2] + gamma0[nGauss] * zb3[2];
					xc = rs[0]; yc = rs[1]; zc = rs[2];
					disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
					r = sqrt(pow(disx, 2) + pow(disy, 2));
					//Incident E has R component;
					EincX.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincX.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
					EincY.real((cos(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincY.imag((sin(-Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
					EincZ.real(0.0);
					EincZ.imag(0.0);
					//cout << "EincX = " << EincX << '\t' << "EincY = " << EincY << '\t' << "EincZ = " << EincZ << endl;
					vector<complex<double>> Einc{ EincX ,EincY ,EincZ };
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x(Einc[ii0].real() * Einc[ii0].real() + Einc[ii0].imag() * Einc[ii0].imag(), 0.0);
						Eincr.real(Eincr.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eincr.imag((Eincr.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]));
					}

					Ex = (0, 0), Ey = (0, 0), Ez = (0, 0);
					for (int ii = 0; ii != 12; ++ii) {
						edge_glb = el[ith].Eedge_GBNO[ii];
						node11 = edge_node_local[ii][0] - 1; node12 = edge_node_local[ii][1] - 1;

						if (edge_glb != 0) {
							//!Output field
							L1 = el[ith].a[node11] + el[ith].b[node11] * xc + el[ith].c[node11] * yc + el[ith].d[node11] * zc;
							L2 = el[ith].a[node12] + el[ith].b[node12] * xc + el[ith].c[node12] * yc + el[ith].d[node12] * zc;

							if (ii <= 5) {
								//N0[0] = el[ith].Eedge_GBNO[edge_loc + 12] * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) - L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, edge_loc) / (36.d0 * Ve(el) * *2)
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] - L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * pow(el[ith].Ve, 2))); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] - L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] - L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							else {
								//N1 = Eedge_GBNO(el, ii + 12) * (L1 * (/ b(el, node12), c(el, node12), d(el, node12) / ) + L2 * (/ b(el, node11), c(el, node11), d(el, node11) / )) * l(el, ii) / (36.d0 * Ve(el) * *2);
								N0[0].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].b[node12] + L2 * el[ith].b[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[0].imag(0.0);
								N0[1].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].c[node12] + L2 * el[ith].c[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[1].imag(0.0);
								N0[2].real(el[ith].Eedge_GBNO[ii + 12] * (L1 * el[ith].d[node12] + L2 * el[ith].d[node11]) * el[ith].length[ii] / (36.0 * el[ith].Ve * el[ith].Ve)); N0[2].imag(0.0);
							}
							Ex += Xsol[edge_glb - 1] * N0[0];
							Ey += Xsol[edge_glb - 1] * N0[1];
							Ez += Xsol[edge_glb - 1] * N0[2];
						}
					}
					Etot[0] = Ex; Etot[1] = Ey; Etot[2] = Ez;
					for (int ii0 = 0; ii0 < 3; ii0++) {
						complex<double> x((Etot[ii0].real()) * Einc[ii0].real() + (Etot[ii0].imag()) * Einc[ii0].imag(), \
							- (Etot[ii0].real()) * Einc[ii0].imag() + (Etot[ii0].imag()) * Einc[ii0].real());
						Eout.real(Eout.real() + x.real() * el[ith].face[nn].Area * weight0[nGauss]);
						Eout.imag(Eout.imag() + x.imag() * el[ith].face[nn].Area * weight0[nGauss]);
					}
				}
			}
		}
	}
	//S161 = Eout / Eincr;
	return 0;
}