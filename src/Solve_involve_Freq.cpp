

#include"Solve_involve_Freq.h"
using namespace Eigen;
using namespace std;

void  Solve_involve_Freq(Element* el_subdomain, int myid) {

	for (int nfreq = 0; Freq < fre_end; ++nfreq) {
		cout << "start" << endl;
		Freq = ((double)(nfreq)) * fre_step + fre_start;
		omega = 2.0 * pi * Freq;
		port_num = 1;
		//Solve_E_H_boundary(el_subdomain,myid);
		if (method == 1) {
			cout << " method = " << 1 << endl;
			cout << "method1"<< endl;
			Solve_E_H_boundary(el_subdomain, myid);
		}
		else if (method == 2) {
			cout << " method = " << 2 << endl;
			Solve_E_H_boundary_method2(el_subdomain, myid);
		}
		else if (method == 3) {
			cout << " method = " << 3 << endl;
			Solve_E_H_boundary_method3(el_subdomain, myid);
		}
		else if (method == -1) {
			cout << " method = " << -1 << endl;
			Solve_E_H_boundary_Iteration(el_subdomain, myid);
		}
		else if (method == -2) {
			cout << " method = " << -2 << endl;
			Solve_E_H_boundary_iter_long(el_subdomain, myid);
		}
		else {
			cout << " method = " << 100 << endl;
			Solve_E_H_boundary(el_subdomain, myid);
		}
		/*Solve_E_H_boundary(el_subdomain, myid);*/
		//Out_field_subdomain(el_subdomain, nfreq);

		Out_field_subdomain_new(el_subdomain, nfreq, myid, 0);

		int of1, ofn1;
		double nx, ny, nz;
		double  P_rad = 0;
		double U_max = 0;

		time_t start_FETI, end_FETI;
		start_FETI = time(NULL);


//		if (abs(Freq - 0.3e9) < 1e-6) {
//			int nQuads = 6;
//			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
//			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
//			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
//			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
//
//
//			double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
//			int num_sample = 361 * 361;
//			VectorXd Efield_far;
//			//Vector3d Ur, zb_node1, zb_node2, zb_node3, ncrosNj1, Ms1;
//			Vector3d Ur;
//			//Vector3cd ncrosNj, Js1, Noperator, Loperator;
//			Vector3cd Noperator, Loperator;
//			Efield_far.resize(num_sample);
//			Efield_far.setZero();
//			int nfar_cnt = 0;
//			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
//			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
//			complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
//
//			for (int nphi = 0; nphi < 1; nphi++) {
//				for (int ntheta = 0; ntheta < 360; ntheta++) {
//
//
//					double d_phi = nphi + 0.5; double d_ntheta = ntheta + 0.5;
//					phi0 = d_phi * pi / 180.0; theta0 = d_ntheta * pi / 180.0;
//					Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
//					Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
//
//#pragma omp  parallel for
//					for (int el = 0; el < num_element; el++) {
//						Vector3d  zb_node1, zb_node2, zb_node3, ncrosNj1, Ms1;
//
//
//
//						Vector3cd ncrosNj, Js1;
//						for (int nn = 0; nn < 4; nn++) {
//
//							of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
//							//of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
//							nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
//
//							//nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
//							if (ofn1 == -7) {
//								zb_node1(0) = el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
//								zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
//								zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
//
//
//								//zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
//								//zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
//								//zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
//								for (int pp = 0; pp < nQuads; pp++) {
//									Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
//									Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
//									Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
//									kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
//									for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//										node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
//
//										edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt];
//										Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										//edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
//										//Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										//Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										//Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										omegaepsil.real(0); omegaepsil.imag(-omega * mur);
//										ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
//										Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
//										if (edge_glb1 != 0) {
//											for (int ii0 = 0; ii0 < 3; ii0++) {
//												cossin.real(cos(kr_cosphi));
//												cossin.imag(sin(kr_cosphi));
//
//												//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//#	   pragma omp critical
//												{
//													Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//												}
//												/*Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);*/
//
//
//
//											}
//										}
//									}
//									for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//										edge_loc1 = face_edge[nn][edge_cnt] - 1;
//										edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1];
//										node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//
//										Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
//										Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
//
//
//
//										//edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
//										//node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//										//Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
//										//Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
//										if (edge_loc1 < 6) {
//											Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
//											Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
//											Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
//
//											//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
//											//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
//											//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
//										}
//										else {
//
//											Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
//											Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
//											Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
//
//
//
//											//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
//											//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
//											//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
//										}
//										ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
//										Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
//										if (edge_glb1 != 0) {
//											for (int ii0 = 0; ii0 < 3; ii0++) {
//												cossin.real(cos(kr_cosphi));
//												cossin.imag(sin(kr_cosphi));
//												//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//#	   pragma omp critical
//												{
//													Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//												}
//												/*Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];*/
//											}
//										}
//									}
//								}
//							}
//						}
//					}
//
//
//
//					Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
//					Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
//					Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
//					Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
//					Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
//					Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
//
//
//					complex<double>* Etheta_n = new complex<double>[num_domain];
//					complex<double>* Ephi_n = new complex<double>[num_domain];
//
//					MPI_Gather(&Etheta, 1, MPI_DOUBLE_COMPLEX, Etheta_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//					MPI_Gather(&Ephi, 1, MPI_DOUBLE_COMPLEX, Ephi_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//					if (myid == 0) {
//
//						Etheta.imag(0.0); Ephi.imag(0.0);
//						Etheta.real(0.0); Ephi.real(0.0);
//
//						for (int i = 0; i < num_domain; i++) {
//							Etheta.imag(Etheta.imag() + Etheta_n[i].imag());
//							Etheta.real(Etheta.real() + Etheta_n[i].real());
//							Ephi.imag(Ephi.imag() + Ephi_n[i].imag());
//							Ephi.real(Ephi.real() + Ephi_n[i].real());
//						}
//
//						P_rad += (norm(Etheta) + norm(Ephi)) * sin(phi0) * 1 / 180 * pi * 1 / 180 * pi;
//						U_max = U_max > (norm(Etheta) + norm(Ephi)) ? U_max : (norm(Etheta) + norm(Ephi));
//
//					}
//
//
//					Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
//					nfar_cnt++;
//				}
//			}
//
//
//
//
//			double E_max = abs(Efield_far(0));
//			for (int i = 0; i < num_sample; i++) {
//				E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
//			}
//
//
//			ofstream ofs111("Direction_map_0_16_245_phi3.txt");
//			for (int i = 0; i < num_sample; i++) {
//				ofs111 << Efield_far(i) / E_max << endl;
//			}
//		}


		/*  time long*/

//		if (abs(Freq - 0.3e9) < 1e-6) {
//			int nQuads = 6;
//			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
//			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
//			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
//			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
//
//			double* Efield_far_mine = new double[180 * 360];
//			for (int i = 0; i < 180 * 360; i++) {
//				Efield_far_mine[i] = 0;
//			}
//
//			double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
//			int num_sample = 361 * 361;
//			VectorXd Efield_far;
//			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
//			Vector3cd ncrosNj, Js1, Noperator, Loperator;
//			Efield_far.resize(num_sample);
//			Efield_far.setZero();
//			int nfar_cnt = 0;
//			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
//			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
//			complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
//
//			for (int nphi = 0; nphi < 180; nphi++) {
//				for (int ntheta = 0; ntheta < 360; ntheta++) {
//
//
//					double d_phi = nphi + 0.5; double d_ntheta = ntheta + 0.5;
//					phi0 = d_phi * pi / 180.0; theta0 = d_ntheta * pi / 180.0;
//					Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
//					Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
//
////#pragma omp  parallel for
//					for (int el = 0; el < num_element_subdomain; el++) {
//						for (int nn = 0; nn < 4; nn++) {
//
//							of1=el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
//							//of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
//							nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
//
//							//nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
//							if (ofn1 == -7) {
//								zb_node1(0)=el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
//								zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
//								zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
//
//
//								//zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
//								//zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
//								//zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
//								for (int pp = 0; pp < nQuads; pp++) {
//									Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
//									Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
//									Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
//									kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
//									for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//										node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
//
//										edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt + 12];
//										Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//										//edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
//										//Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										//Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										//Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
//										omegaepsil.real(0); omegaepsil.imag(-omega * mur);
//										ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
//										Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
//										if (edge_glb1 != 0) {
//											for (int ii0 = 0; ii0 < 3; ii0++) {
//												cossin.real(cos(kr_cosphi));
//												cossin.imag(sin(kr_cosphi));
//												//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//
//												Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//
//											}
//										}
//									}
//									for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//										edge_loc1 = face_edge[nn][edge_cnt] - 1;
//										edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1 + 12];
//										node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//
//										Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
//										Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
//
//
//
//										//edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
//										//node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//										//Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
//										//Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
//										if (edge_loc1 < 6) {
//											Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
//											Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
//											Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
//
//											//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
//											//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
//											//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
//										}
//										else {
//
//											Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
//											Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
//											Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
//
//
//
//											//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
//											//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
//											//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
//										}
//										ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
//										Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
//										if (edge_glb1 != 0) {
//											for (int ii0 = 0; ii0 < 3; ii0++) {
//												cossin.real(cos(kr_cosphi));
//												cossin.imag(sin(kr_cosphi));
//												//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//												Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//											}
//										}
//									}
//								}
//							}
//						}
//					}
//
//
//
//					Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
//					Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
//					Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
//					Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
//					Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
//					Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
//
//
//					complex<double>* Etheta_n = new complex<double>[num_domain];
//					complex<double>* Ephi_n = new complex<double>[num_domain];
//
//					MPI_Gather(&Etheta, 1, MPI_DOUBLE_COMPLEX, Etheta_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//					MPI_Gather(&Ephi, 1, MPI_DOUBLE_COMPLEX, Ephi_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//					if (myid == 0) {
//
//						Etheta.imag(0.0); Ephi.imag(0.0);
//						Etheta.real(0.0); Ephi.real(0.0);
//
//						for (int i = 0; i < num_domain; i++) {
//							Etheta.imag(Etheta.imag() + Etheta_n[i].imag());
//							Etheta.real(Etheta.real() + Etheta_n[i].real());
//							Ephi.imag(Ephi.imag() + Ephi_n[i].imag());
//							Ephi.real(Ephi.real() + Ephi_n[i].real());
//						}
//						Efield_far_mine[nphi + ntheta * 180] = (norm(Etheta) + norm(Ephi));
//						P_rad += (norm(Etheta) + norm(Ephi))*sin(phi0)*1/180*pi*1/180*pi;
//						U_max = U_max > (norm(Etheta) + norm(Ephi)) ? U_max : (norm(Etheta) + norm(Ephi));
//
//					}
//
//
//					Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
//					nfar_cnt++;
//				}
//			}
//
//			ofstream ofs111("Direction_map_Efield_far_mine.txt");
//			for (int i = 0; i < 180 * 360; i++) {
//				ofs111 << Efield_far_mine[i] << endl;
//			}
//
//			ofstream ofs1225("Direction_map_P_rad_mine.txt");
//
//			ofs1225 << P_rad << endl;
//
//			double E_max = abs(Efield_far(0));
//			for (int i = 0; i < num_sample; i++) {
//				E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
//			}
//
//
//			ofstream ofs1224("Direction_map_0_16_245_phi3.txt");
//			for (int i = 0; i < num_sample; i++) {
//				ofs1224 << Efield_far(i) / E_max << endl;
//			}
//		}
//
//
//
//		end_FETI = time(NULL);
//
//		double time_FETI = (double)(end_FETI - start_FETI);
//
//
//		if (myid == 0) {
//			cout << " 4 * pi * U_max / P_rad =  " << 4 * pi * U_max / P_rad << endl;
//			cout << "4 * pi * U_max / P_rad cost time is  " << time_FETI << endl;
//
//		}
		
		/*  time long*/





		time_t start_pre_pro, end_pre_pro;
		start_pre_pro = time(NULL);

		//VectorXd Noperator_free_of_theta_phi;
		//VectorXd Loperator_free_of_theta_phi;
		int num_face_ABC = 0;


		for (int el = 0; el < num_element_subdomain; el++) {
			for (int nn = 0; nn < 4; nn++) {
				of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
				if (ofn1 == -7) {
					num_face_ABC++;
				}
			}
		}
		//cout << "myid is " << myid << "    num_face_ABC  is " << num_face_ABC << endl;
		MatrixXd Noperator_free_of_theta_phi;
		MatrixXd Loperator_free_of_theta_phi;

		//Noperator_free_of_theta_phi.resize(0);
		Noperator_free_of_theta_phi.resize(num_face_ABC, 3);
		Loperator_free_of_theta_phi.resize(num_face_ABC, 3);
		Noperator_free_of_theta_phi.setZero();
		Loperator_free_of_theta_phi.setZero();

		MatrixXd Xgs_free_of_theta_phi;
		MatrixXd Ygs_free_of_theta_phi;
		MatrixXd Zgs_free_of_theta_phi;
		Xgs_free_of_theta_phi.resize(num_face_ABC, 6);
		Ygs_free_of_theta_phi.resize(num_face_ABC, 6);
		Zgs_free_of_theta_phi.resize(num_face_ABC, 6);
		Xgs_free_of_theta_phi.setZero();
		Ygs_free_of_theta_phi.setZero();
		Zgs_free_of_theta_phi.setZero();

		MatrixXd JS;
		JS.resize(num_face_ABC, 6 * 3);
		JS.setZero();
		MatrixXd MS;
		MS.resize(num_face_ABC, 6 * 3);
		MS.setZero();


		complex<double>* X_sol_temp = new complex<double>[num_face_ABC * 6];
		complex<double>* X_sol_temp1 = new complex<double>[num_face_ABC * 6];


		int count_face_ABC = -1;
		complex<double>omegaepsil2;
		omegaepsil2.real(0); omegaepsil2.imag(-omega * mur);


		for (int el = 0; el < num_element_subdomain; el++) {
			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
			//double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
			//double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
			//double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
			//double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };

			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
			Vector3d ncrosNj, Js1;
			for (int nn = 0; nn < 4; nn++) {
				of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
				nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
				if (ofn1 == -7) {
					count_face_ABC++;
					zb_node1(0) = el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
					zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
					zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
					for (int pp = 0; pp < 6; pp++) {
						Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
						Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
						Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
						Xgs_free_of_theta_phi(count_face_ABC, pp) = Xgs;
						Ygs_free_of_theta_phi(count_face_ABC, pp) = Ygs;
						Zgs_free_of_theta_phi(count_face_ABC, pp) = Zgs;
						for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
							node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
							edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt + 12];
							Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							//omegaepsil.real(0); omegaepsil.imag(-omega * mur);
							//ncrosNj(0) = ((ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = ((nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = ((nx * Njy - ny * Njx)) / omegaepsil;
							ncrosNj(0) = ((ny * Njz - nz * Njy)) / 1; ncrosNj(1) = ((nz * Njx - nx * Njz)) / 1; ncrosNj(2) = ((nx * Njy - ny * Njx)) / 1;
							Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
							JS(count_face_ABC, edge_cnt * 3 + 0) = edge_sign1 * Js1[0]; JS(count_face_ABC, edge_cnt * 3 + 1) = edge_sign1 * Js1[1]; JS(count_face_ABC, edge_cnt * 3 + 2) = edge_sign1 * Js1[2];
							if (edge_glb1 != 0) {
								X_sol_temp[count_face_ABC * 6 + edge_cnt] = Xsol[edge_glb1 - 1] / omegaepsil2;
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//cossin.real(cos(kr_cosphi));
									//cossin.imag(sin(kr_cosphi));
									//Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
									Noperator_free_of_theta_phi(count_face_ABC, ii0) = (el_subdomain[el].face[nn].Area * 1 * 1) * 1 * 1 * 1;

								}
							}
							else {
								X_sol_temp[count_face_ABC * 6 + edge_cnt] = 0.0;
							}
						}
						for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
							edge_loc1 = face_edge[nn][edge_cnt] - 1;
							edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1 + 12];
							node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
							Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
							Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
							if (edge_loc1 < 6) {
								Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
								Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
								Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
							}
							else {
								Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
								Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
								Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
							}
							ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
							Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
							MS(count_face_ABC, edge_cnt * 3 + 0) = edge_sign1 * Ms1[0]; MS(count_face_ABC, edge_cnt * 3 + 1) = edge_sign1 * Ms1[1]; MS(count_face_ABC, edge_cnt * 3 + 2) = edge_sign1 * Ms1[2];
							if (edge_glb1 != 0) {
								X_sol_temp1[count_face_ABC * 6 + edge_cnt] = Xsol[edge_glb1 - 1];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//cossin.real(cos(kr_cosphi));
									//cossin.imag(sin(kr_cosphi));
									/*Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];*/
									Loperator_free_of_theta_phi(count_face_ABC, ii0) = (el_subdomain[el].face[nn].Area * 1 * 1 * 1) * 1 * 1;

								}
							}
							else {

								X_sol_temp1[count_face_ABC * 6 + edge_cnt] = 0.0;
							}
						}
					}
				}
			}
		}
		end_pre_pro = time(NULL);

		double time_pre_pro = (double)(end_pre_pro - start_pre_pro);

		cout << "time_pre_pro cost time is  " << time_pre_pro << endl;



		/*   Serial  */


		if (abs(Freq - 1000e9) < 1e-6) {
			int nQuads = 6;
			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
			double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
			int num_sample = 361 * 361;
			double* Efield_far_mine = new double[180 * 360];
			for (int i = 0; i < 180 * 360; i++) {
				Efield_far_mine[i] = 0;
			}
			VectorXd Efield_far;
			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
			Vector3cd ncrosNj, Js1, Noperator, Loperator;
			Efield_far.resize(num_sample);
			Efield_far.setZero();
			int nfar_cnt = 0;
			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
			complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
			int num_of_nphi = 360; int num_of_ntheta = 180;
			for (int nphi = 0; nphi < num_of_nphi; nphi++) {
				for (int ntheta = 89; ntheta < 90; ntheta++) {
					double d_phi = nphi + 0.5; double d_ntheta = ntheta + 1;
					phi0 = d_phi * pi / 180.0; theta0 = d_ntheta * pi / 180.0;
					Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
					Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
					for (int face_ABC_number = 0; face_ABC_number < num_face_ABC; face_ABC_number++) {
						for (int pp = 0; pp < 6; pp++) {
							Xgs = Xgs_free_of_theta_phi(face_ABC_number, pp);
							Ygs = Ygs_free_of_theta_phi(face_ABC_number, pp);
							Zgs = Zgs_free_of_theta_phi(face_ABC_number, pp);
							kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
								Js1[0] = JS(face_ABC_number, edge_cnt * 3 + 0);
								Js1[1] = JS(face_ABC_number, edge_cnt * 3 + 1);
								Js1[2] = JS(face_ABC_number, edge_cnt * 3 + 2);
								complex<double>X_temp1 = X_sol_temp[face_ABC_number * 6 + edge_cnt];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									cossin.real(cos(kr_cosphi));
									cossin.imag(sin(kr_cosphi));
									//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
									Noperator(ii0) += complex<double>(Noperator_free_of_theta_phi(face_ABC_number, ii0) * 1 * weight0[pp]) * cossin * X_temp1 * Js1(ii0);
								}
							}
							//count_face_ABC;
							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
								Ms1[0] = MS(face_ABC_number, edge_cnt * 3 + 0);
								Ms1[1] = MS(face_ABC_number, edge_cnt * 3 + 1);
								Ms1[2] = MS(face_ABC_number, edge_cnt * 3 + 2);
								complex<double>X_temp1 = X_sol_temp1[face_ABC_number * 6 + edge_cnt];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									cossin.real(cos(kr_cosphi));
									cossin.imag(sin(kr_cosphi));
									//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
									Loperator(ii0) += complex<double>(Loperator_free_of_theta_phi(face_ABC_number, ii0) * Ms1(ii0) * 1 * weight0[pp]) * cossin * X_temp1;
								}
							}
						}
					}
					//for (int el = 0; el < num_element_subdomain; el++) {
					//	for (int nn = 0; nn < 4; nn++) {
					//		of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
					//		//of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
					//		nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
					//		//nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
					//		if (ofn1 == -7) {
					//			zb_node1(0) = el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
					//			zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
					//			zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
					//			//zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
					//			//zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
					//			//zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
					//			for (int pp = 0; pp < nQuads; pp++) {
					//				Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
					//				Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
					//				Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
					//				kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
					//				for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
					//					node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
					//					edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt + 12];
					//					Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
					//					Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
					//					Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
					//					//edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
					//					//Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
					//					//Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
					//					//Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
					//					omegaepsil.real(0); omegaepsil.imag(-omega * mur);
					//					ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
					//					Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
					//					if (edge_glb1 != 0) {
					//						for (int ii0 = 0; ii0 < 3; ii0++) {
					//							cossin.real(cos(kr_cosphi));
					//							cossin.imag(sin(kr_cosphi));
					//							//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
					//							Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
					//						}
					//					}
					//				}
					//				for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
					//					edge_loc1 = face_edge[nn][edge_cnt] - 1;
					//					edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1 + 12];
					//					node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
					//					Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
					//					Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
					//					//edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
					//					//node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
					//					//Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
					//					//Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
					//					if (edge_loc1 < 6) {
					//						Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
					//						Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
					//						Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
					//						//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
					//						//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
					//						//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
					//					}
					//					else {
					//						Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
					//						Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
					//						Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
					//						//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
					//						//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
					//						//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
					//					}
					//					ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
					//					Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
					//					if (edge_glb1 != 0) {
					//						for (int ii0 = 0; ii0 < 3; ii0++) {
					//							cossin.real(cos(kr_cosphi));
					//							cossin.imag(sin(kr_cosphi));
					//							//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
					//							Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
					//						}
					//					}
					//				}
					//			}
					//		}
					//	}
					//}
					Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
					Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
					Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
					Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
					Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
					Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
					complex<double>* Etheta_n = new complex<double>[num_domain];
					complex<double>* Ephi_n = new complex<double>[num_domain];
					MPI_Gather(&Etheta, 1, MPI_DOUBLE_COMPLEX, Etheta_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
					MPI_Gather(&Ephi, 1, MPI_DOUBLE_COMPLEX, Ephi_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
					if (myid == 0) {
						Etheta.imag(0.0); Ephi.imag(0.0);
						Etheta.real(0.0); Ephi.real(0.0);
						for (int i = 0; i < num_domain; i++) {
							Etheta.imag(Etheta.imag() + Etheta_n[i].imag());
							Etheta.real(Etheta.real() + Etheta_n[i].real());
							Ephi.imag(Ephi.imag() + Ephi_n[i].imag());
							Ephi.real(Ephi.real() + Ephi_n[i].real());
						}
						Efield_far_mine[nphi + ntheta * num_of_nphi] = (norm(Etheta) + norm(Ephi));
						P_rad += (norm(Etheta) + norm(Ephi)) * sin(theta0) * 1 / 180 * pi ;
						//P_rad += (norm(Etheta) + norm(Ephi)) * sin(phi0) * 1 / 180 * pi * 1 / 180 * pi;
						U_max = U_max > (norm(Etheta) + norm(Ephi)) ? U_max : (norm(Etheta) + norm(Ephi));
					}
					Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
					nfar_cnt++;
				}
			}
			ofstream ofs1224("Direction_map_Efield_far_mine.txt");
			for (int i = 0; i < num_of_nphi * num_of_ntheta; i++) {
				ofs1224 << Efield_far_mine[i] << endl;
			}
			ofstream ofs1225("Direction_map_P_rad_mine.txt");
			ofs1225 << P_rad << endl;
			double E_max = abs(Efield_far(0));
			for (int i = 0; i < num_sample; i++) {
				E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
			}
			ofstream ofs111("Direction_map_0_16_245_phi3.txt");
			for (int i = 0; i < num_sample; i++) {
				ofs111 << Efield_far(i) / E_max << endl;
			}
		}


		/*  Serial  */




		/*   Parallel*/
//		if (abs(Freq - 0.3e9) < 1e-6) {
//			int nQuads = 6;
//			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
//			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
//			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
//			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
//			double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
//			int num_sample = 361 * 361;
//			double* Efield_far_mine = new double[180 * 360];
//			for (int i = 0; i < 180 * 360; i++) {
//				Efield_far_mine[i] = 0;
//			}
//			VectorXd Efield_far;
//			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
//			Vector3cd ncrosNj, Js1, Noperator, Loperator;
//			Efield_far.resize(num_sample);
//			Efield_far.setZero();
//			int nfar_cnt = 0;
//			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
//			/*double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;*/
//			double phi0, theta0, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
//
//
//			complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
//			for (int nphi = 0; nphi < 1; nphi++) {
//#pragma omp  parallel for
//				for (int ntheta = 0; ntheta < 360; ntheta++) {
//					double d_phi = nphi + 0.5; double d_ntheta = ntheta + 0.5;
//					phi0 = d_phi * pi / 180.0; theta0 = d_ntheta * pi / 180.0;
//					Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
//					Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
////#pragma omp  parallel for num_threads(10)
//					for (int face_ABC_number = 0; face_ABC_number < num_face_ABC; face_ABC_number++) {
//						double  Xgs, Ygs, Zgs, kr_cosphi;
//						for (int pp = 0; pp < 6; pp++) {
//							//Xgs = Xgs_free_of_theta_phi(face_ABC_number, pp);
//							//Ygs = Ygs_free_of_theta_phi(face_ABC_number, pp);
//							//Zgs = Zgs_free_of_theta_phi(face_ABC_number, pp);
//#	   pragma omp critical
//							{
//								Xgs = Xgs_free_of_theta_phi(face_ABC_number, pp);
//								Ygs = Ygs_free_of_theta_phi(face_ABC_number, pp);
//								Zgs = Zgs_free_of_theta_phi(face_ABC_number, pp);
//								kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
//							}
//							/*kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));*/
//							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
////#	   pragma omp critical
////								{
////									Js1[0] = JS(face_ABC_number, edge_cnt * 3 + 0);
////									Js1[1] = JS(face_ABC_number, edge_cnt * 3 + 1);
////									Js1[2] = JS(face_ABC_number, edge_cnt * 3 + 2);
////									complex<double>X_temp1 = X_sol_temp[face_ABC_number * 6 + edge_cnt];
////								}
//								Js1[0] = JS(face_ABC_number, edge_cnt * 3 + 0);
//								Js1[1] = JS(face_ABC_number, edge_cnt * 3 + 1);
//								Js1[2] = JS(face_ABC_number, edge_cnt * 3 + 2);
//								complex<double>X_temp1 = X_sol_temp[face_ABC_number * 6 + edge_cnt];
//								for (int ii0 = 0; ii0 < 3; ii0++) {
//#	   pragma omp critical
//									{
//										cossin.real(cos(kr_cosphi));
//										cossin.imag(sin(kr_cosphi));
//										Noperator(ii0) += complex<double>(Noperator_free_of_theta_phi(face_ABC_number, ii0) * 1 * weight0[pp]) * cossin * X_temp1 * Js1(ii0);
//									}
//
//									//cossin.real(cos(kr_cosphi));
//									//cossin.imag(sin(kr_cosphi));
//									//Noperator(ii0) += complex<double>(Noperator_free_of_theta_phi(face_ABC_number, ii0) * 1 * weight0[pp]) * cossin * X_temp1 * Js1(ii0);
//								}
//							}
//							//count_face_ABC;
//							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
////#	   pragma omp critical
////								{
////									Ms1[0] = MS(face_ABC_number, edge_cnt * 3 + 0);
////									Ms1[1] = MS(face_ABC_number, edge_cnt * 3 + 1);
////									Ms1[2] = MS(face_ABC_number, edge_cnt * 3 + 2);
////									complex<double>X_temp1 = X_sol_temp1[face_ABC_number * 6 + edge_cnt];
////								}
//								Ms1[0] = MS(face_ABC_number, edge_cnt * 3 + 0);
//								Ms1[1] = MS(face_ABC_number, edge_cnt * 3 + 1);
//								Ms1[2] = MS(face_ABC_number, edge_cnt * 3 + 2);
//								complex<double>X_temp1 = X_sol_temp1[face_ABC_number * 6 + edge_cnt];
//								for (int ii0 = 0; ii0 < 3; ii0++) {
//#	   pragma omp critical
//									{
//										cossin.real(cos(kr_cosphi));
//										cossin.imag(sin(kr_cosphi));
//										Loperator(ii0) += complex<double>(Loperator_free_of_theta_phi(face_ABC_number, ii0) * Ms1(ii0) * 1 * weight0[pp]) * cossin * X_temp1;
//									}
//									//cossin.real(cos(kr_cosphi));
//									//cossin.imag(sin(kr_cosphi));
//									//Loperator(ii0) += complex<double>(Loperator_free_of_theta_phi(face_ABC_number, ii0) * Ms1(ii0) * 1 * weight0[pp]) * cossin * X_temp1;
//								}
//							}
//						}
//					}
//					//for (int el = 0; el < num_element_subdomain; el++) {
//					//	for (int nn = 0; nn < 4; nn++) {
//					//		of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
//					//		//of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
//					//		nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
//					//		//nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
//					//		if (ofn1 == -7) {
//					//			zb_node1(0) = el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
//					//			zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
//					//			zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
//					//			//zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
//					//			//zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
//					//			//zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
//					//			for (int pp = 0; pp < nQuads; pp++) {
//					//				Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
//					//				Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
//					//				Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
//					//				kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
//					//				for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//					//					node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
//					//					edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt + 12];
//					//					Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//					//					Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//					//					Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
//					//					//edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
//					//					//Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
//					//					//Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
//					//					//Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
//					//					omegaepsil.real(0); omegaepsil.imag(-omega * mur);
//					//					ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
//					//					Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
//					//					if (edge_glb1 != 0) {
//					//						for (int ii0 = 0; ii0 < 3; ii0++) {
//					//							cossin.real(cos(kr_cosphi));
//					//							cossin.imag(sin(kr_cosphi));
//					//							//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//					//							Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
//					//						}
//					//					}
//					//				}
//					//				for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
//					//					edge_loc1 = face_edge[nn][edge_cnt] - 1;
//					//					edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1 + 12];
//					//					node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//					//					Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
//					//					Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
//					//					//edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
//					//					//node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
//					//					//Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
//					//					//Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
//					//					if (edge_loc1 < 6) {
//					//						Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
//					//						Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
//					//						Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
//					//						//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
//					//						//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
//					//						//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
//					//					}
//					//					else {
//					//						Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
//					//						Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
//					//						Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
//					//						//Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
//					//						//Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
//					//						//Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
//					//					}
//					//					ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
//					//					Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
//					//					if (edge_glb1 != 0) {
//					//						for (int ii0 = 0; ii0 < 3; ii0++) {
//					//							cossin.real(cos(kr_cosphi));
//					//							cossin.imag(sin(kr_cosphi));
//					//							//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//					//							Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
//					//						}
//					//					}
//					//				}
//					//			}
//					//		}
//					//	}
//					//}
//					Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
//					Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
//					Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
//					Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
//					Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
//					Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
//					complex<double>* Etheta_n = new complex<double>[num_domain];
//					complex<double>* Ephi_n = new complex<double>[num_domain];
//					MPI_Gather(&Etheta, 1, MPI_DOUBLE_COMPLEX, Etheta_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//					MPI_Gather(&Ephi, 1, MPI_DOUBLE_COMPLEX, Ephi_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//					if (myid == 0) {
//						Etheta.imag(0.0); Ephi.imag(0.0);
//						Etheta.real(0.0); Ephi.real(0.0);
//						for (int i = 0; i < num_domain; i++) {
//							Etheta.imag(Etheta.imag() + Etheta_n[i].imag());
//							Etheta.real(Etheta.real() + Etheta_n[i].real());
//							Ephi.imag(Ephi.imag() + Ephi_n[i].imag());
//							Ephi.real(Ephi.real() + Ephi_n[i].real());
//						}
//						Efield_far_mine[nphi + ntheta * 180] = (norm(Etheta) + norm(Ephi));
//						P_rad += (norm(Etheta) + norm(Ephi)) * sin(phi0) * 1 / 180 * pi * 1 / 180 * pi;
//						U_max = U_max > (norm(Etheta) + norm(Ephi)) ? U_max : (norm(Etheta) + norm(Ephi));
//					}
//					Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
//					nfar_cnt++;
//				}
//			}
//			ofstream ofs1224("Direction_map_Efield_far_mine.txt");
//			for (int i = 0; i < 180 * 360; i++) {
//				ofs1224 << Efield_far_mine[i] << endl;
//			}
//			ofstream ofs1225("Direction_map_P_rad_mine.txt");
//			ofs1225 << P_rad << endl;
//			double E_max = abs(Efield_far(0));
//			for (int i = 0; i < num_sample; i++) {
//				E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
//			}
//			ofstream ofs111("Direction_map_0_16_245_phi3.txt");
//			for (int i = 0; i < num_sample; i++) {
//				ofs111 << Efield_far(i) / E_max << endl;
//			}
//		}

		/*   Parallel*/

		end_FETI = time(NULL);

		double time_FETI = (double)(end_FETI - start_FETI);


		if (myid == 0) {
			//cout << " 4 * pi * U_max / P_rad =  " << 2 * pi * U_max / P_rad << endl;
			//cout << " 4 * pi * U_max / P_rad =  " << 4 * pi * U_max / P_rad << endl;
			//cout << "4 * pi * U_max / P_rad cost time is  " << time_FETI << endl;

		}


	//	if (abs(Freq - 2.45e9) < 1e-6) {
	//		int nQuads = 6;
	//		double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//		double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//		double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//		double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	//		double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
	//		int num_sample = 361 * 361;
	//		VectorXd Efield_far;
	//		Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
	//		Vector3cd ncrosNj, Js1, Noperator, Loperator;
	//		Efield_far.resize(num_sample);
	//		Efield_far.setZero();
	//		int nfar_cnt = 0;
	//		int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
	//		double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz;
	//		complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
	//		for (int nphi = 90; nphi < 91; nphi++) {
	//			for (int ntheta = 0; ntheta < 361; ntheta++) {
	//				phi0 = nphi * pi / 180.0; theta0 = ntheta * pi / 180.0;
	//				Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
	//				Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
	//				for (int el = 0; el < num_element; el++) {
	//					for (int nn = 0; nn < 4; nn++) {
	//						of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
	//						nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
	//						if (ofn1 == -7) {
	//							zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
	//							zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
	//							zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
	//							for (int pp = 0; pp < nQuads; pp++) {
	//								Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
	//								Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
	//								Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
	//								kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
	//									edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
	//									Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									omegaepsil.real(0); omegaepsil.imag(-omega * mur);
	//									ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
	//									Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
	//										}
	//									}
	//								}
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
	//									node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
	//									Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
	//									Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
	//									if (edge_loc1 < 6) {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
	//									}
	//									else {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
	//									}
	//									ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
	//									Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
	//										}
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//				Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
	//				Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
	//				Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
	//				Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
	//				Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
	//				Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
	//				Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
	//				nfar_cnt++;
	//			}
	//		}
	//		double E_max = abs(Efield_far(0));
	//		for (int i = 0; i < num_sample; i++) {
	//			E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
	//		}
	//		ofstream ofs111("Direction_map_90_16_245_phi3.txt");
	//		for (int i = 0; i < num_sample; i++) {
	//			ofs111 << Efield_far(i) / E_max << endl;
	//		}
	//	}


	//for (int i_num = 0; i_num < 1; i_num++) {
	//	Freq = (i_num) * 0.05e9 + 2.45e9;
	//	cout << Freq << endl;
	//	omega = 2 * pi * Freq;
	//	for (int i_1 = 0; i_1 < num_unknown; i_1++) {
	//		Xsol[i_1].real(0.0); Xsol[i_1].imag(0.0);
	//	}
	//	for (int i_2 = 0; i_2 < nnz; i_2++) {
	//		acoo1[i_2].real = acoo_1[i_2].real() + acoo_2[i_2].real() * omega + acoo_3[i_2].real() * omega * omega;
	//		acoo1[i_2].imag = acoo_1[i_2].imag() + acoo_2[i_2].imag() * omega + acoo_3[i_2].imag() * omega * omega;
	//		m[i_2] = m_1[i_2] + m_2[i_2] * omega + m_3[i_2] * omega * omega;
	//	}
	//	clock_t start_time = clock();
	//	S_para1(Vector);
	//	clock_t end_time = clock();
	//	cout << "per time use: " << (end_time - start_time) / (double)CLOCKS_PER_SEC << " s" << endl;
	//	ofs_t << Freq << " " << (end_time - tt_start0) / (double)CLOCKS_PER_SEC << endl;
	//	cout << i_num << endl;
	//	if (abs(Freq - 2.45e9) < 1e-6) {
	//		int nQuads = 6;
	//		double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//		double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//		double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//		double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	//		double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
	//		int num_sample = 361 * 361;
	//		VectorXd Efield_far;
	//		Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
	//		Vector3cd ncrosNj, Js1, Noperator, Loperator;
	//		Efield_far.resize(num_sample);
	//		Efield_far.setZero();
	//		int nfar_cnt = 0;
	//		int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
	//		double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz;
	//		complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
	//		for (int nphi = 0; nphi < 1; nphi++) {
	//			for (int ntheta = 0; ntheta < 361; ntheta++) {
	//				phi0 = nphi * pi / 180.0; theta0 = ntheta * pi / 180.0;
	//				Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
	//				Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
	//				for (int el = 0; el < num_element; el++) {
	//					for (int nn = 0; nn < 4; nn++) {
	//						of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
	//						nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
	//						if (ofn1 == -7) {
	//							zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
	//							zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
	//							zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
	//							for (int pp = 0; pp < nQuads; pp++) {
	//								Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
	//								Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
	//								Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
	//								kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
	//									edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
	//									Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									omegaepsil.real(0); omegaepsil.imag(-omega * mur);
	//									ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
	//									Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
	//										}
	//									}
	//								}
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
	//									node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
	//									Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
	//									Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
	//									if (edge_loc1 < 6) {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
	//									}
	//									else {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
	//									}
	//									ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
	//									Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
	//										}
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//				Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
	//				Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
	//				Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
	//				Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
	//				Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
	//				Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
	//				Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
	//				nfar_cnt++;
	//			}
	//		}
	//		double E_max = abs(Efield_far(0));
	//		for (int i = 0; i < num_sample; i++) {
	//			E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
	//		}
	//		ofstream ofs111("Direction_map_0_16_245_phi3.txt");
	//		for (int i = 0; i < num_sample; i++) {
	//			ofs111 << Efield_far(i) / E_max << endl;
	//		}
	//	}
	//	if (abs(Freq - 2.45e9) < 1e-6) {
	//		int nQuads = 6;
	//		double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//		double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//		double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//		double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	//		double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
	//		int num_sample = 361 * 361;
	//		VectorXd Efield_far;
	//		Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
	//		Vector3cd ncrosNj, Js1, Noperator, Loperator;
	//		Efield_far.resize(num_sample);
	//		Efield_far.setZero();
	//		int nfar_cnt = 0;
	//		int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
	//		double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz;
	//		complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
	//		for (int nphi = 90; nphi < 91; nphi++) {
	//			for (int ntheta = 0; ntheta < 361; ntheta++) {
	//				phi0 = nphi * pi / 180.0; theta0 = ntheta * pi / 180.0;
	//				Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
	//				Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
	//				for (int el = 0; el < num_element; el++) {
	//					for (int nn = 0; nn < 4; nn++) {
	//						of1 = opp[el][nn]; ofn1 = opp[el][nn + 4];
	//						nx = Vector[0][nn][el]; ny = Vector[1][nn][el]; nz = Vector[2][nn][el];
	//						if (ofn1 == -7) {
	//							zb_node1(0) = zb(ver(el, face_node[nn][0] - 1) - 1, 0); zb_node1(1) = zb(ver(el, face_node[nn][0] - 1) - 1, 1); zb_node1(2) = zb(ver(el, face_node[nn][0] - 1) - 1, 2);
	//							zb_node2(0) = zb(ver(el, face_node[nn][1] - 1) - 1, 0); zb_node2(1) = zb(ver(el, face_node[nn][1] - 1) - 1, 1); zb_node2(2) = zb(ver(el, face_node[nn][1] - 1) - 1, 2);
	//							zb_node3(0) = zb(ver(el, face_node[nn][2] - 1) - 1, 0); zb_node3(1) = zb(ver(el, face_node[nn][2] - 1) - 1, 1); zb_node3(2) = zb(ver(el, face_node[nn][2] - 1) - 1, 2);
	//							for (int pp = 0; pp < nQuads; pp++) {
	//								Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
	//								Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
	//								Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
	//								kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
	//									edge_glb1 = Eedge_GBNO[el][edge_cnt]; edge_sign1 = Eedge_GBNO[el][edge_cnt + 12];
	//									Njx = 2.0 * l(el, edge_cnt) * (c(el, node11) * d(el, node12) - d(el, node11) * c(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njy = 2.0 * l(el, edge_cnt) * (d(el, node11) * b(el, node12) - b(el, node11) * d(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									Njz = 2.0 * l(el, edge_cnt) * (b(el, node11) * c(el, node12) - c(el, node11) * b(el, node12)) / (36.0 * Ve(el) * Ve(el));
	//									omegaepsil.real(0); omegaepsil.imag(-omega * mur);
	//									ncrosNj(0) = (complex<double>(ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = (complex<double>(nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = (complex<double>(nx * Njy - ny * Njx)) / omegaepsil;
	//									Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
	//										}
	//									}
	//								}
	//								for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
	//									edge_loc1 = face_edge[nn][edge_cnt] - 1;   edge_glb1 = Eedge_GBNO[el][edge_loc1]; edge_sign1 = Eedge_GBNO[el][edge_loc1 + 12];
	//									node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
	//									Li11 = (a(el, node11) + b(el, node11) * Xgs + c(el, node11) * Ygs + d(el, node11) * Zgs);
	//									Li12 = (a(el, node12) + b(el, node12) * Xgs + c(el, node12) * Ygs + d(el, node12) * Zgs);
	//									if (edge_loc1 < 6) {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) - Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) - Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) - Li12 * d(el, node11));
	//									}
	//									else {
	//										Njx1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * b(el, node12) + Li12 * b(el, node11));
	//										Njy1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * c(el, node12) + Li12 * c(el, node11));
	//										Njz1 = l(el, edge_loc1) / pow(6.0 * Ve(el), 2) * (Li11 * d(el, node12) + Li12 * d(el, node11));
	//									}
	//									ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
	//									Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
	//									if (edge_glb1 != 0) {
	//										for (int ii0 = 0; ii0 < 3; ii0++) {
	//											cossin.real(cos(kr_cosphi));
	//											cossin.imag(sin(kr_cosphi));
	//											Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
	//										}
	//									}
	//								}
	//							}
	//						}
	//					}
	//				}
	//				Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
	//				Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
	//				Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
	//				Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
	//				Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
	//				Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
	//				Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
	//				nfar_cnt++;
	//			}
	//		}
	//		double E_max = abs(Efield_far(0));
	//		for (int i = 0; i < num_sample; i++) {
	//			E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
	//		}
	//		ofstream ofs111("Direction_map_90_16_245_phi3.txt");
	//		for (int i = 0; i < num_sample; i++) {
	//			ofs111 << Efield_far(i) / E_max << endl;
	//		}
	//	}
	//}



		

	}


	if (fre_pattern != 1e15) {

		Freq = fre_pattern;
		omega = 2.0 * pi * Freq;

		Solve_E_H_boundary2(el_subdomain, myid);

		int of1, ofn1;
		double nx, ny, nz;
		double  P_rad = 0;
		double U_max = 0;

		time_t start_FETI, end_FETI;
		start_FETI = time(NULL);

		time_t start_pre_pro, end_pre_pro;
		start_pre_pro = time(NULL);

		int num_face_ABC = 0;


		for (int el = 0; el < num_element_subdomain; el++) {
			for (int nn = 0; nn < 4; nn++) {
				of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
				if (ofn1 == -7) {
					num_face_ABC++;
				}
			}
		}
		cout << "myid is " << myid << "    num_face_ABC  is " << num_face_ABC << endl;
		MatrixXd Noperator_free_of_theta_phi;
		MatrixXd Loperator_free_of_theta_phi;

		//Noperator_free_of_theta_phi.resize(0);
		Noperator_free_of_theta_phi.resize(num_face_ABC, 3);
		Loperator_free_of_theta_phi.resize(num_face_ABC, 3);
		Noperator_free_of_theta_phi.setZero();
		Loperator_free_of_theta_phi.setZero();

		MatrixXd Xgs_free_of_theta_phi;
		MatrixXd Ygs_free_of_theta_phi;
		MatrixXd Zgs_free_of_theta_phi;
		Xgs_free_of_theta_phi.resize(num_face_ABC, 6);
		Ygs_free_of_theta_phi.resize(num_face_ABC, 6);
		Zgs_free_of_theta_phi.resize(num_face_ABC, 6);
		Xgs_free_of_theta_phi.setZero();
		Ygs_free_of_theta_phi.setZero();
		Zgs_free_of_theta_phi.setZero();

		MatrixXd JS;
		JS.resize(num_face_ABC, 6 * 3);
		JS.setZero();
		MatrixXd MS;
		MS.resize(num_face_ABC, 6 * 3);
		MS.setZero();


		complex<double>* X_sol_temp = new complex<double>[num_face_ABC * 6];
		complex<double>* X_sol_temp1 = new complex<double>[num_face_ABC * 6];


		int count_face_ABC = -1;
		complex<double>omegaepsil2;
		omegaepsil2.real(0); omegaepsil2.imag(-omega * mur);


		for (int el = 0; el < num_element_subdomain; el++) {
			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
			//double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
			//double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
			//double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
			//double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };

			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
			Vector3d ncrosNj, Js1;
			for (int nn = 0; nn < 4; nn++) {
				of1 = el_subdomain[el].face[nn].opp[0]; ofn1 = el_subdomain[el].face[nn].opp[1];
				nx = el_subdomain[el].face[nn].N_Vector[0]; ny = el_subdomain[el].face[nn].N_Vector[1]; nz = el_subdomain[el].face[nn].N_Vector[2];
				if (ofn1 == -7) {
					count_face_ABC++;
					zb_node1(0) = el_subdomain[el].node[face_node[nn][0] - 1].zb[0]; zb_node1(1) = el_subdomain[el].node[face_node[nn][0] - 1].zb[1]; zb_node1(2) = el_subdomain[el].node[face_node[nn][0] - 1].zb[2];
					zb_node2(0) = el_subdomain[el].node[face_node[nn][1] - 1].zb[0]; zb_node2(1) = el_subdomain[el].node[face_node[nn][1] - 1].zb[1]; zb_node2(2) = el_subdomain[el].node[face_node[nn][1] - 1].zb[2];
					zb_node3(0) = el_subdomain[el].node[face_node[nn][2] - 1].zb[0]; zb_node3(1) = el_subdomain[el].node[face_node[nn][2] - 1].zb[1]; zb_node3(2) = el_subdomain[el].node[face_node[nn][2] - 1].zb[2];
					for (int pp = 0; pp < 6; pp++) {
						Xgs = alpha0[pp] * zb_node1(0) + beta0[pp] * zb_node2(0) + gamma0[pp] * zb_node3(0);
						Ygs = alpha0[pp] * zb_node1(1) + beta0[pp] * zb_node2(1) + gamma0[pp] * zb_node3(1);
						Zgs = alpha0[pp] * zb_node1(2) + beta0[pp] * zb_node2(2) + gamma0[pp] * zb_node3(2);
						Xgs_free_of_theta_phi(count_face_ABC, pp) = Xgs;
						Ygs_free_of_theta_phi(count_face_ABC, pp) = Ygs;
						Zgs_free_of_theta_phi(count_face_ABC, pp) = Zgs;
						for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
							node11 = edge_node_local[edge_cnt][0] - 1; node12 = edge_node_local[edge_cnt][1] - 1;
							edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_cnt]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_cnt + 12];
							Njx = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].c[node11] * el_subdomain[el].d[node12] - el_subdomain[el].d[node11] * el_subdomain[el].c[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							Njy = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].d[node11] * el_subdomain[el].b[node12] - el_subdomain[el].b[node11] * el_subdomain[el].d[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							Njz = 2.0 * el_subdomain[el].length[edge_cnt] * (el_subdomain[el].b[node11] * el_subdomain[el].c[node12] - el_subdomain[el].c[node11] * el_subdomain[el].b[node12]) / (36 * el_subdomain[el].Ve * el_subdomain[el].Ve);
							//omegaepsil.real(0); omegaepsil.imag(-omega * mur);
							//ncrosNj(0) = ((ny * Njz - nz * Njy)) / omegaepsil; ncrosNj(1) = ((nz * Njx - nx * Njz)) / omegaepsil; ncrosNj(2) = ((nx * Njy - ny * Njx)) / omegaepsil;
							ncrosNj(0) = ((ny * Njz - nz * Njy)) / 1; ncrosNj(1) = ((nz * Njx - nx * Njz)) / 1; ncrosNj(2) = ((nx * Njy - ny * Njx)) / 1;
							Js1[0] = ncrosNj(0); Js1[1] = ncrosNj(1); Js1[2] = ncrosNj(2);
							JS(count_face_ABC, edge_cnt * 3 + 0) = edge_sign1 * Js1[0]; JS(count_face_ABC, edge_cnt * 3 + 1) = edge_sign1 * Js1[1]; JS(count_face_ABC, edge_cnt * 3 + 2) = edge_sign1 * Js1[2];
							if (edge_glb1 != 0) {
								X_sol_temp[count_face_ABC * 6 + edge_cnt] = Xsol[edge_glb1 - 1] / omegaepsil2;
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//cossin.real(cos(kr_cosphi));
									//cossin.imag(sin(kr_cosphi));
									//Noperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
									Noperator_free_of_theta_phi(count_face_ABC, ii0) = (el_subdomain[el].face[nn].Area * 1 * 1) * 1 * 1 * 1;

								}
							}
							else {
								X_sol_temp[count_face_ABC * 6 + edge_cnt] = 0.0;
							}
						}
						for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
							edge_loc1 = face_edge[nn][edge_cnt] - 1;
							edge_glb1 = el_subdomain[el].Eedge_GBNO[edge_loc1]; edge_sign1 = el_subdomain[el].Eedge_GBNO[edge_loc1 + 12];
							node11 = edge_node_local[edge_loc1][0] - 1; node12 = edge_node_local[edge_loc1][1] - 1;
							Li11 = (el_subdomain[el].a[node11] + el_subdomain[el].b[node11] * Xgs + el_subdomain[el].c[node11] * Ygs + el_subdomain[el].d[node11] * Zgs);
							Li12 = (el_subdomain[el].a[node12] + el_subdomain[el].b[node12] * Xgs + el_subdomain[el].c[node12] * Ygs + el_subdomain[el].d[node12] * Zgs);
							if (edge_loc1 < 6) {
								Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] - Li12 * el_subdomain[el].b[node11]);
								Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] - Li12 * el_subdomain[el].c[node11]);
								Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] - Li12 * el_subdomain[el].d[node11]);
							}
							else {
								Njx1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].b[node12] + Li12 * el_subdomain[el].b[node11]);
								Njy1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].c[node12] + Li12 * el_subdomain[el].c[node11]);
								Njz1 = el_subdomain[el].length[edge_loc1] / pow(6.0 * el_subdomain[el].Ve, 2) * (Li11 * el_subdomain[el].d[node12] + Li12 * el_subdomain[el].d[node11]);
							}
							ncrosNj1(0) = ny * Njz1 - nz * Njy1; ncrosNj1(1) = nz * Njx1 - nx * Njz1;  ncrosNj1(2) = nx * Njy1 - ny * Njx1;
							Ms1(0) = -ncrosNj1(0); Ms1(1) = -ncrosNj1(1); Ms1(2) = -ncrosNj1(2);
							MS(count_face_ABC, edge_cnt * 3 + 0) = edge_sign1 * Ms1[0]; MS(count_face_ABC, edge_cnt * 3 + 1) = edge_sign1 * Ms1[1]; MS(count_face_ABC, edge_cnt * 3 + 2) = edge_sign1 * Ms1[2];
							if (edge_glb1 != 0) {
								X_sol_temp1[count_face_ABC * 6 + edge_cnt] = Xsol[edge_glb1 - 1];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//cossin.real(cos(kr_cosphi));
									//cossin.imag(sin(kr_cosphi));
									/*Loperator(ii0) += complex<double>(el_subdomain[el].face[nn].Area * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];*/
									Loperator_free_of_theta_phi(count_face_ABC, ii0) = (el_subdomain[el].face[nn].Area * 1 * 1 * 1) * 1 * 1;

								}
							}
							else {

								X_sol_temp1[count_face_ABC * 6 + edge_cnt] = 0.0;
							}
						}
					}
				}
			}
		}
		end_pre_pro = time(NULL);

		double time_pre_pro = (double)(end_pre_pro - start_pre_pro);

		cout << "time_pre_pro cost time is  " << time_pre_pro << endl;



		/*   Serial  */


		if (1) {
			int nQuads = 6;
			double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
			double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
			double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
			double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
			double kwave = 2.0 * pi * Freq * sqrt(epsilon0 * mur);
			int num_sample = 361 * 361;
			double* Efield_far_mine = new double[180 * 360];
			for (int i = 0; i < 180 * 360; i++) {
				Efield_far_mine[i] = 0;
			}
			VectorXd Efield_far;
			Vector3d Ur, zb_node1, zb_node2, zb_node3, zb_gaus, ncrosNj1, Ms1;
			Vector3cd ncrosNj, Js1, Noperator, Loperator;
			Efield_far.resize(num_sample);
			Efield_far.setZero();
			int nfar_cnt = 0;
			int face_count, gaus_count, of1, edge_loc1, edge_loc2, edge_glb1, edge_glb2, edge_sign1, edge_sign2, node11, node12;
			double phi0, theta0, Xgs, Ygs, Zgs, kr_cosphi, Li11, Li12, Njx, Njy, Njz, Njx1, Njy1, Njz1;
			complex<double> cossin, omegaepsil, Ntheta, Nphi, Ltheta, Lphi, Etheta, Ephi;
			int num_of_nphi = 360; int num_of_ntheta = 180;
			for (int nphi = 0; nphi < num_of_nphi; nphi++) {
				for (int ntheta = 0; ntheta < 180; ntheta++) {
					double d_phi = nphi + 0.5; double d_ntheta = ntheta + 1;
					phi0 = d_phi * pi / 180.0; theta0 = d_ntheta * pi / 180.0;
					Ur(0) = sin(theta0) * cos(phi0); Ur(1) = sin(theta0) * sin(phi0); Ur(2) = cos(theta0);
					Noperator(0) = 0.0; Noperator(1) = 0.0; Noperator(2) = 0.0; Loperator(0) = 0.0; Loperator(1) = 0.0; Loperator(2) = 0.0;
					for (int face_ABC_number = 0; face_ABC_number < num_face_ABC; face_ABC_number++) {
						for (int pp = 0; pp < 6; pp++) {
							Xgs = Xgs_free_of_theta_phi(face_ABC_number, pp);
							Ygs = Ygs_free_of_theta_phi(face_ABC_number, pp);
							Zgs = Zgs_free_of_theta_phi(face_ABC_number, pp);
							kr_cosphi = kwave * (Xgs * Ur(0) + Ygs * Ur(1) + Zgs * Ur(2));
							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
								Js1[0] = JS(face_ABC_number, edge_cnt * 3 + 0);
								Js1[1] = JS(face_ABC_number, edge_cnt * 3 + 1);
								Js1[2] = JS(face_ABC_number, edge_cnt * 3 + 2);
								complex<double>X_temp1 = X_sol_temp[face_ABC_number * 6 + edge_cnt];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									cossin.real(cos(kr_cosphi));
									cossin.imag(sin(kr_cosphi));
									//Noperator(ii0) += complex<double>(Area(el, nn) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1] * Js1(ii0);
									Noperator(ii0) += complex<double>(Noperator_free_of_theta_phi(face_ABC_number, ii0) * 1 * weight0[pp]) * cossin * X_temp1 * Js1(ii0);
								}
							}
							//count_face_ABC;
							for (int edge_cnt = 0; edge_cnt < 6; edge_cnt++) {
								Ms1[0] = MS(face_ABC_number, edge_cnt * 3 + 0);
								Ms1[1] = MS(face_ABC_number, edge_cnt * 3 + 1);
								Ms1[2] = MS(face_ABC_number, edge_cnt * 3 + 2);
								complex<double>X_temp1 = X_sol_temp1[face_ABC_number * 6 + edge_cnt];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									cossin.real(cos(kr_cosphi));
									cossin.imag(sin(kr_cosphi));
									//Loperator(ii0) += complex<double>(Area(el, nn) * Ms1(ii0) * edge_sign1 * weight0[pp]) * cossin * Xsol[edge_glb1 - 1];
									Loperator(ii0) += complex<double>(Loperator_free_of_theta_phi(face_ABC_number, ii0) * Ms1(ii0) * 1 * weight0[pp]) * cossin * X_temp1;
								}
							}
						}
					}
					Ntheta = Noperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Noperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Noperator(2) * complex<double>(-sin(theta0));
					Nphi = Noperator(0) * complex<double>(-sin(phi0)) + Noperator(1) * complex<double>(cos(phi0));
					Ltheta = Loperator(0) * complex<double>(cos(theta0) * cos(phi0)) + Loperator(1) * complex<double>(cos(theta0) * sin(phi0)) + Loperator(2) * complex<double>(-sin(theta0));
					Lphi = Loperator(0) * complex<double>(-sin(phi0)) + Loperator(1) * complex<double>(cos(phi0));
					Etheta = (complex<double>(120.0 * pi) * Ntheta + Lphi);
					Ephi = (Ltheta - complex<double>(120.0 * pi) * Nphi);
					complex<double>* Etheta_n = new complex<double>[num_domain];
					complex<double>* Ephi_n = new complex<double>[num_domain];
					MPI_Gather(&Etheta, 1, MPI_DOUBLE_COMPLEX, Etheta_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
					MPI_Gather(&Ephi, 1, MPI_DOUBLE_COMPLEX, Ephi_n, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
					if (myid == 0) {
						Etheta.imag(0.0); Ephi.imag(0.0);
						Etheta.real(0.0); Ephi.real(0.0);
						for (int i = 0; i < num_domain; i++) {
							Etheta.imag(Etheta.imag() + Etheta_n[i].imag());
							Etheta.real(Etheta.real() + Etheta_n[i].real());
							Ephi.imag(Ephi.imag() + Ephi_n[i].imag());
							Ephi.real(Ephi.real() + Ephi_n[i].real());
						}
						Efield_far_mine[nphi + ntheta * num_of_nphi] = (norm(Etheta) + norm(Ephi));
						//P_rad += (norm(Etheta) + norm(Ephi)) * sin(theta0) * 1 / 180 * pi;
						P_rad += (norm(Etheta) + norm(Ephi)) * sin(theta0) * 1 / 180 * pi * 1 / 180 * pi;
						U_max = U_max > (norm(Etheta) + norm(Ephi)) ? U_max : (norm(Etheta) + norm(Ephi));
					}
					Efield_far(nfar_cnt) = sqrt(norm(Etheta) + norm(Ephi));
					nfar_cnt++;
				}
			}
			if (myid == 0) {
				ofstream ofs1224("Direction_map_Efield_fare.txt");
				for (int i = 0; i < num_of_nphi * num_of_ntheta; i++) {
					ofs1224 << Efield_far_mine[i] << endl;
				}
				ofstream ofs1225("Direction_map_P_rad.txt");
				ofs1225 << P_rad << endl;
			}

			double E_max = abs(Efield_far(0));
			for (int i = 0; i < num_sample; i++) {
				E_max = E_max > abs(Efield_far(i)) ? E_max : abs(Efield_far(i));
			}
		}

		if (myid == 0) {
			cout << " 4 * pi * U_max / P_rad =  " << 2 * pi * U_max / P_rad << endl;
			//cout << " 4 * pi * U_max / P_rad =  " << 4 * pi * U_max / P_rad << endl;
			//cout << "4 * pi * U_max / P_rad cost time is  " << time_FETI << endl;

		}
		/*  Serial  */
	}


}


