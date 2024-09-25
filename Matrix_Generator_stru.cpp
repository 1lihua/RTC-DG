
#include "Matrix_Generator_stru.h"

using namespace Eigen;
using namespace std;


double Flabel_stru(int vnode1, int vnode2);
int Result_Out1();
double Ffeeder = 2000;
double Fscrew = 900;
double Afeeder = pi * 20e-3 * 20e-3;
double Awasher = pi / 4.0 * (pow(7e-3, 2) - pow(3.75e-3, 2));


vector<Triplet<double>> TriList_stru;

void* Matrix_Generator1_stru(int myid) {
	int unknown_dm, nn2, final_row, nnz_dm, nrhs, cnt_n1, end1, start1, cnt_n2, end2, start2, final_col, mm2;
	double* bb, * x, * a_pro, * b_pro, * pro;
	int maxfct, mnum, mtype, phase1, phase2, phase3, n, msglvl, error, comm, info;
	int iparm[64]{ 0 };
	int pt[64]{ 0 };;
	int* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
	iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;

	for (int n1 = myid; n1 < myid+1; n1++) {

		unknown_dm = num_unKnown_subdomain_stru[n1][0] + num_unKnown_subdomain_stru[n1][1];
		nnz_dm = num_nzero_Pmatrix_stru[n1];
		final_row = r_dm_stru[n1];
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 };
		mkl_dcsrcoo(job, &unknown_dm, p_pardiso_stru , jp_stru , ip_stru , &nnz_dm, acoo1_stru , rowind1_stru, colind1_stru , &info);
		cnt_n1 = Node_Boundary_stru[n1];
		bb = new double[cnt_n1 * unknown_dm]();
		x = new double[cnt_n1 * unknown_dm];
		nrhs = cnt_n1;
		for (int i = 0; i < cnt_n1; i++) {
			bb[i * unknown_dm + i + unknown_dm - cnt_n1] = 1.0;//Identity matrix
		}
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_stru , ip_stru , jp_stru , perm2, &nrhs, iparm, &msglvl, bb, x, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_stru , ip_stru , jp_stru , perm2, &nrhs, iparm, &msglvl, bb, x, &error);
		delete[]bb;
		end1 = num_unKnown_subdomain_stru[n1][0] + num_unKnown_subdomain_stru[n1][1];
		start1 = end1 - cnt_n1 + 1;
		a_pro = new double[cnt_n1 * cnt_n1]();
#		pragma omp parallel for 
		for (int i = 0; i < cnt_n1 * unknown_dm; i++) {
			if (i % unknown_dm + 1 >= start1 && i % unknown_dm + 1 <= end1 && i / unknown_dm + 1 >= 1 && i / unknown_dm + 1 <= cnt_n1) {
				a_pro[cnt_n1 * (i % unknown_dm + 1 - start1) + i / unknown_dm + 1 - 1] = x[i];
			}
		}
		delete[] x;

		double* fbr_temp = new double[unknown_dm];
		//fbb_stru = new double[cnt_n1];
		nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_stru, ip_stru, jp_stru, perm2, &nrhs, iparm, &msglvl, fbr_stru, fbr_temp, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_stru, ip_stru, jp_stru, perm2, &nrhs, iparm, &msglvl, fbr_stru, fbr_temp, &error);

		for (int i = 0; i < cnt_n1; i++) {
			fbb_stru[i] = fbr_temp[i + unknown_dm - cnt_n1];
		}
		delete[] fbr_temp;

		//cout << "TTTTT1" << endl;
		final_col = 0;
		nn2 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_C_stru[n1][n2] != 0) {
				mm2 = nn2 + nnz_C_stru[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Node_Boundary_stru[n2];
					end2 = num_unKnown_subdomain_stru[n2][0] + num_unKnown_subdomain_stru[n2][1];
					start2 = end2 - cnt_n2 + 1;
					b_pro = new double[cnt_n1 * cnt_n2]();
					pro = new double[cnt_n1 * cnt_n2]();
					for (int i = nn2; i < mm2; i++) {
						if (mrow_stru[i] >= start1 && mrow_stru[i] <= end1 && mcol_stru[i] >= start2 && mcol_stru[i] <= end2) {
							b_pro[cnt_n2 * (mrow_stru[i] - start1) + mcol_stru[i] - start2] = m_stru[i];
						}
					}
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, 1, a_pro, cnt_n1, b_pro, cnt_n2, 0, pro, cnt_n2);
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
					TriList_stru.insert(TriList_stru.end(), Tri_tmp.begin(), Tri_tmp.end());
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Node_Boundary_stru[n2];
		}
	}
}

void* Solver_stru(int myid) {
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
		int unknown_dm = num_unKnown_subdomain_stru[n1][0]+ num_unKnown_subdomain_stru[n1][1];
		int cnt_n2 = 0, mm3, nn3 = 0, unknown_dm2;
		int nn4 = 0, mm4, nn2 = 0;
		int nzero1 = 0, row1 = 0;
		double* b_pro, * pro, * fb_tmp;
		double* fbr_stru_temp = new  double[unknown_dm];

		for (int i = 0; i < unknown_dm; i++) {
			fbr_stru_temp[i] = fbr_stru[i];
		}
		int cnt_n1 = Node_Boundary_stru[(n1)];
		for (int n2 = 0; n2 < num_domain; n2++) {
			cnt_n2 = Node_Boundary_stru[(n2)];
			mm4 = nn4 + cnt_n2;
			if (nnz_C_stru[n1][n2] != 0) {
				mm3 = nn3 + nnz_C_stru[n1][n2];
				unknown_dm2 = num_unKnown_subdomain_stru[n2][0] + num_unKnown_subdomain_stru[n2][1];
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

					for (int i = 0; i < cnt_n2; i++) {
						fb_tmp[i] = xx[nn4 + i];
					}
					for (int i = nn3; i < mm3; i++) {
						if (mrow_stru[i] >= start1 && mrow_stru[i] <= end1 && mcol_stru[i] >= start2 && mcol_stru[i] <= end2) {

							b_pro[cnt_n2 * (mrow_stru[i] - start1) + mcol_stru[i] - start2] = m_stru[i];
						}


						//b_pro[unknown_dm2 * (mrow[i] - 1) + mcol[i] - 1] = m_T_q[i];
					}
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, 1, b_pro, cnt_n2, fb_tmp, 1, 0, pro, 1);
					for (int i = 0; i < cnt_n1; i++) {
						fbr_stru_temp[ i + unknown_dm - cnt_n1] -= pro[i];
					}
					delete[] b_pro;
					delete[] pro;
					delete[] fb_tmp;

					//b_pro = new double[unknown_dm * unknown_dm2]();
					//pro = new double[unknown_dm];
					//fb_tmp = new double[unknown_dm2]();
					//for (int i = 0; i < cnt_n2; i++) {
					//	fb_tmp[unknown_dm2 - cnt_n2 + i] = xx[nn4 + i];
					//}
					//for (int i = nn3; i < mm3; i++) {
					//	b_pro[unknown_dm2 * (mrow_stru[i] - 1) + mcol_stru[i] - 1] = m_stru[i];
					//}
					//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, unknown_dm, 1, unknown_dm2, 1, b_pro, unknown_dm2, fb_tmp, 1, 0, pro, 1);
					//for (int i = 0; i < unknown_dm; i++) {
					//	fbr_stru_temp[nn2 + i] -= pro[i];
					//}
					//delete[] b_pro;
					//delete[] pro;
					//delete[] fb_tmp;

				}
				nn3 = mm3;
			}
			nn4 = mm4;
		}
		int nrhs = 1;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso_stru , ip_stru , jp_stru , perm, &nrhs, iparm, &msglvl, fbr_stru_temp , Xstru , &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso_stru , ip_stru , jp_stru , perm, &nrhs, iparm, &msglvl, fbr_stru_temp , Xstru , &error);
	}
	return 0;
}

void Direct_solution_stru1(Element* element, int myid) {

	fbb_stru = new double[num_unKnown_stb]();
	cout << "myid is " << myid << " num_unKnown_stb = " << num_unKnown_stb << " Node_Boundary_stru[myid] = " << Node_Boundary_stru[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	Matrix_Generator1_stru(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	//cout << "myid is " << myid << "  Volumn solve time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);


	//cout << "tt4" << endl;
	int mat_size = TriList_stru.size();   //  nnz of Zi*Cij (j=1,2,3...,n)

	//cout << "myid is " << myid << "  mat_size is " << mat_size << endl;

	double* m_final = new double[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i] = (TriList_stru[i].value());
		rowm[i] = TriList_stru[i].row() + 1;
		colm[i] = TriList_stru[i].col() + 1;
	}
	TriList_stru.clear();
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
		face_mat_size += nnz_ZC_process[j] + Node_Boundary_stru[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size_stru = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Node_Boundary_stru[j];
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
	double* fbb_temp = new double[Node_Boundary_stru[myid]];
	for (int j = 0; j < Node_Boundary_stru[myid]; j++) {
		fbb_temp[j] = fbb_stru[j];
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Node_Boundary_stru[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Node_Boundary_stru[myid], MPI_DOUBLE, fbb_total, Node_Boundary_stru, address_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
	//cout << "tt5" << endl;

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
		nrhs = 1;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		delete[] m_pardiso; delete[] im; delete[] jm; delete[]  fbb_total;

	}


	end = time(NULL);

	double time1 = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}




	Xstru = new double[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]];
	MPI_Bcast(xx, num_unk_boundary, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	Solver_stru(myid);
	delete[] xx;
	end_out = time(NULL);

	double time_out = (double)(end_out - start_out);
	if (myid == 0) {
		cout << "time_face_solve_stru is " << time1 << endl;
		ofstream ofs1119("data_DG_RTC_waveguide_stru.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unKnown_subdomain_stru[0][0] + num_unKnown_subdomain_stru[0][1] << ',' << num_unk_boundary << ',' << time_out << ',' \
			<< time1 << "s" << ',' << num_nzero_Pmatrix_stru[myid] << ',' << face_mat_size << ',' << Node_Boundary_stru[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << mat_size << endl;

	}


	//if (myid == 0) {
	//	ofstream ofs1119("uvw.csv");
	//	for (int i = 0; i < num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]; i++) {
	//		ofs1119 << Xstru[i] << endl;
	//	}
	//}


	int un_T_temp = 0;
	int* Accumulated_T_temp = new int[num_domain];
	int* unknown_T_temp = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		Accumulated_T_temp[j] = 0;
		unknown_T_temp[j] = num_unKnown_subdomain_stru[j][0];

	}

	for (int j = 1; j < num_domain; j++) {
		Accumulated_T_temp[j] = Accumulated_T_temp[j - 1] + num_unKnown_subdomain_stru[j - 1][0];
	}
	cout << "tt7" << endl;

	for (int j = 0; j < num_domain; j++) {
		un_T_temp += num_unKnown_subdomain_stru[j][0];
	}

	double* X_T = new double[un_T_temp];
	double* X_h_temp = new double[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_u = new int[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_v = new int[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_w = new int[num_unKnown_subdomain_stru[myid][0]];
	for (int j = 0; j < num_unKnown_subdomain_stru[myid][0]; j++) {
		node_global_temp_u[j] = 0;
		node_global_temp_v[j] = 0;
		node_global_temp_w[j] = 0;
	}

	int* V_T_u = new int[un_T_temp];
	int* V_T_v = new int[un_T_temp];
	int* V_T_w = new int[un_T_temp];

	for (int j = 0; j < num_element_subdomain; j++) {
		for (int ii = 0; ii < 4; ii++) {
			if (element[j].node[ii].unknown_uvw[0] != 0) {

				node_global_temp_u[element[j].node[ii].unknown_uvw[0] - 1] = element[j].node[ii].ver;

			}
			if (element[j].node[ii].unknown_uvw[1] != 0) {

				node_global_temp_v[element[j].node[ii].unknown_uvw[1] - 1] = element[j].node[ii].ver;

			}
			if (element[j].node[ii].unknown_uvw[2] != 0) {

				node_global_temp_w[element[j].node[ii].unknown_uvw[2] - 1] = element[j].node[ii].ver;

			}

		}
	}

	for (int j = 0; j < num_unKnown_subdomain_stru[myid][0]; j++) {
		X_h_temp[j] = Xstru[j];
	}

	MPI_Gatherv(X_h_temp, num_unKnown_subdomain_stru[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_u, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_u, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_v, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_v, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_w, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_w, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);


	delete[]Accumulated_T_temp; delete[]unknown_T_temp;  delete[]node_global_temp_u; delete[]node_global_temp_v; delete[]node_global_temp_w; delete[]X_h_temp; delete[]Xstru;




	double* XX1_u = new double[num_node];
	double* XX1_v = new double[num_node];
	double* XX1_w = new double[num_node];
	for (int i = 0; i < num_node; i++) {
		XX1_u[i] = 0;
		XX1_v[i] = 0;
		XX1_w[i] = 0;
	}


	if (myid == 0) {
		for (int i = 0; i < un_T_temp; i++) {
			int temp_ver = V_T_u[i];
			if (temp_ver != 0) {
				XX1_u[temp_ver - 1] = X_T[i];
			}
			int temp_ver1 = V_T_v[i];
			if (temp_ver1 != 0) {
				XX1_v[temp_ver1 - 1] = X_T[i];
			}
			int temp_ver2 = V_T_w[i];
			if (temp_ver2 != 0) {
				XX1_w[temp_ver2 - 1] = X_T[i];
			}

		}

		ofstream ofs12157("Patch_u.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12157 << XX1_u[i] << endl;
		}
		ofstream ofs12158("Patch_v.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12158 << XX1_v[i] << endl;
		}
		ofstream ofs12159("Patch_w.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12159 << XX1_w[i] << endl;
		}



	}


	delete[]X_T; delete[]XX1_u; delete[]V_T_u; delete[]V_T_v; delete[]V_T_w;




//	vector<Triplet<double>> Tri_tmp(num_unKnown_stb);
//	for (int i = 0; i < num_unKnown_stb; i++) {
//		Tri_tmp[i] = Triplet<double>(i, i, 1);
//	}
//	TriList_stru.insert(TriList_stru.end(), Tri_tmp.begin(), Tri_tmp.end());
//	double* bb, * x;
//	int cnt1 = 0, cnt2 = 0, tot_dm;
//	MKL_INT maxfct, mnum, mtype, phase1, phase2, phase3, nrhs, msglvl, error;
//	MKL_INT iparm[64]{ 0 };
//	int pt[64]{ 0 };
//	MKL_INT* perm2 = nullptr;
//	maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
//	iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;
//	nzero1 = 0; row1 = 0;
//	int cnt_n1, cnt_n2, start1, start2, end1, end2, final_row = 0, final_col = 0;
//	xx = new double[num_unKnown_stb];
//	int size = TriList_stru.size();
//	double* m_final = new double[size];
//	double* m_pardiso1 = new double[size];
//	int* rowm = new int[size]; int* colm = new int[size]; int* im1 = new int[num_unKnown_stb + 1]; int* jm1 = new int[size];
//	nnz = size;
//#	pragma omp parallel for num_threads(96)
//	for (int i = 0; i < size; i++) {
//		m_final[i] = TriList_stru[i].value();
//		rowm[i] = TriList_stru[i].row() + 1;
//		colm[i] = TriList_stru[i].col() + 1;
//	}
//	TriList_stru.clear();
//	nrhs = 1;
//	int job[8] = { 2, 1, 1, 0, 0, 0, 0, 0 };
//	job[4] = nnz;
//	mkl_dcsrcoo(job, &num_unKnown_stb, m_pardiso1, jm1, im1, &nnz, m_final, rowm, colm, &info);
//	delete[] m_final;
//	delete[] rowm;
//	delete[] colm;
//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unKnown_stb, m_pardiso1, im1, jm1, perm2, &nrhs, iparm, &msglvl, fbb_stru, xx, &error);
//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unKnown_stb, m_pardiso1, im1, jm1, perm2, &nrhs, iparm, &msglvl, fbb_stru, xx, &error);
//	delete[] m_pardiso1;
//	delete[] im1;
//	delete[] jm1;
//	pthread_t tids[NUM_THREADS_stru];
//	for (int i = 0; i < NUM_THREADS_stru; ++i) {
//		int ret = pthread_create(&tids[i], NULL, Solver_stru, (void*)i);
//	}
//	for (int i = 0; i < NUM_THREADS_stru; i++) {
//		pthread_join(tids[i], NULL);
//	}
//	delete[] xx;

}


void Direct_solution_stru(Element* element, int myid) {

	fbb_stru = new double[num_unKnown_stb]();
	cout << "myid is " << myid << " num_unKnown_stb = " << num_unKnown_stb << " Node_Boundary_stru[myid] = " << Node_Boundary_stru[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);
	int myid1 = myid;
	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	Matrix_Generator1_stru(myid1);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	//cout << "myid is " << myid << "  Volumn solve time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);

	double time_out = (double)(end_out - start_out);

	//cout << "tt4" << endl;
	int mat_size = TriList_stru.size();   //  nnz of Zi*Cij (j=1,2,3...,n)

	//cout << "myid is " << myid << "  mat_size is " << mat_size << endl;

	double* m_final = new double[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i] = (TriList_stru[i].value());
		rowm[i] = TriList_stru[i].row() + 1;
		colm[i] = TriList_stru[i].col() + 1;
	}
	TriList_stru.clear();
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
		face_mat_size += nnz_ZC_process[j] + Node_Boundary_stru[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << " my id is " << myid << endl;
		cout << "face_mat_size_stru = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Node_Boundary_stru[j];
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
	double* fbb_temp = new double[Node_Boundary_stru[myid]];
	for (int j = 0; j < Node_Boundary_stru[myid]; j++) {
		fbb_temp[j] = fbb_stru[j];
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Node_Boundary_stru[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Node_Boundary_stru[myid], MPI_DOUBLE, fbb_total, Node_Boundary_stru, address_offset, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;

	time_t start, end;

	start = time(NULL);
	//delete[] xx;
	
	xx = new double[num_unk_boundary];
	for (int i = 0; i < num_unk_boundary; i++) {
		xx[i] = 0.0;
	}
	cout << "tt5" << endl;

	//omp_set_num_threads(64);

	if (myid == 0) {
		cout << "face_mat_size_stru123 = " << num_unk_boundary << endl;
		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j] = 1.0;
		}
		ofstream ofs11141("m_final_total1.txt");
		ofstream ofs11142("rowm_total1.txt");
		ofstream ofs11143("colm_total1.txt");
		for (int j = 0; j < face_mat_size; j++) {
			ofs11141 << m_final_total[face_mat_size-1-j] << endl;
			ofs11142 << rowm_total[face_mat_size - 1 - j] << endl;
			ofs11143 << colm_total[face_mat_size - 1 - j] << endl;
		}
		ofstream ofs11144("num_unk_boundary.txt");
		for (int j = 0; j < num_unk_boundary; j++) {
			ofs11144 << fbb_total[j] << endl;
		}

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
		int maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
		int iparm[64]{ 0 };
		int pt[64]{ 0 };;
		int* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
		int phase3 = 12, phase4 = 33;
		//MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
		//MKL_INT iparm[64]{ 0 };
		//int pt[64]{ 0 };
		//for (int el = 0; el < 64; el++) {
		//	pt[el] = 0; iparm[el] = 0;
		//}
		//MKL_INT* perm2 = nullptr;
		//maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
		//iparm[0] = 0;
		iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;
		nrhs = 1;
		cout << "test8888" << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase3, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		//cout << "test9999" << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total, xx, &error);
		delete[] m_pardiso; delete[] im; delete[] jm; delete[]  fbb_total;

	}


	end = time(NULL);

	double time = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}



	if (myid == 0) {
		cout << "time_face_solve_stru is " << time << endl;
		ofstream ofs1119("data_DG_RTC_patch_stru.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unKnown_subdomain_stru[0][0] + num_unKnown_subdomain_stru[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix_stru[myid] << ',' << face_mat_size << ',' << Node_Boundary_stru[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << mat_size << endl;

	}
	Xstru = new double[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]];
	MPI_Bcast(xx, num_unk_boundary, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	Solver_stru(myid);
	delete[] xx;
	if (myid == 0) {
		ofstream ofs1119("uvw.csv");
		for (int i = 0; i < num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]; i++) {
			ofs1119 << Xstru[i] << endl;
		}
	}


	int un_T_temp = 0;
	int* Accumulated_T_temp = new int[num_domain];
	int* unknown_T_temp = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		Accumulated_T_temp[j] = 0;
		unknown_T_temp[j] = num_unKnown_subdomain_stru[j][0];

	}

	for (int j = 1; j < num_domain; j++) {
		Accumulated_T_temp[j] = Accumulated_T_temp[j - 1] + num_unKnown_subdomain_stru[j - 1][0];
	}
	cout << "tt7" << endl;

	for (int j = 0; j < num_domain; j++) {
		un_T_temp += num_unKnown_subdomain_stru[j][0];
	}

	double* X_T = new double[un_T_temp];
	double* von_Mises_temp = new double[un_T_temp];//

	double* X_h_temp = new double[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_u = new int[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_v = new int[num_unKnown_subdomain_stru[myid][0]];
	int* node_global_temp_w = new int[num_unKnown_subdomain_stru[myid][0]];
	for (int j = 0; j < num_unKnown_subdomain_stru[myid][0]; j++) {
		node_global_temp_u[j] = 0;
		node_global_temp_v[j] = 0;
		node_global_temp_w[j] = 0;
	}

	int* V_T_u = new int[un_T_temp];
	int* V_T_v = new int[un_T_temp];
	int* V_T_w = new int[un_T_temp];

	double sigma_six[6];//
	double* Xstruvon_Mises_temp = new double[num_unKnown_subdomain_stru[myid][0]];
	//double* von_Mises_temp = new double[num_node];//
	for (int i = 0; i < num_unKnown_subdomain_stru[myid][0]; i++) {
		Xstruvon_Mises_temp[i] = 0;
	}

	for (int j = 0; j < num_element_subdomain; j++) {
		double TT = 0;
		for (int node_ii = 0; node_ii < 4; node_ii++) {
			int unknown_TT = element[j].node[node_ii].unknown_T;
			if (unknown_TT == 0) {
				TT += Td;
			}
			else {
				TT += Xh[unknown_TT - 1];
			}
		}
		TT = TT / 4.0;
		double u_temp0[4]; double v_temp0[4]; double w_temp0[4];//

		for (int ii = 0; ii < 4; ii++) {
			double u_temp = 0; double v_temp = 0; double w_temp = 0;//
			if (element[j].node[ii].unknown_uvw[0] != 0) {
				node_global_temp_u[element[j].node[ii].unknown_uvw[0] - 1] = element[j].node[ii].ver;
				u_temp = Xstru[element[j].node[ii].unknown_uvw[0] - 1];//
				u_temp0[ii] = Xstru[element[j].node[ii].unknown_uvw[0] - 1];//
			}
			else {
				u_temp0[ii] = 0;
			}
			if (element[j].node[ii].unknown_uvw[1] != 0) {
				node_global_temp_v[element[j].node[ii].unknown_uvw[1] - 1] = element[j].node[ii].ver;
				v_temp = Xstru[element[j].node[ii].unknown_uvw[1] - 1];//
				v_temp0[ii] = Xstru[element[j].node[ii].unknown_uvw[1] - 1];//
			}
			else {
				v_temp0[ii] = 0;
			}

			if (element[j].node[ii].unknown_uvw[2] != 0) {

				node_global_temp_w[element[j].node[ii].unknown_uvw[2] - 1] = element[j].node[ii].ver;
				w_temp = Xstru[element[j].node[ii].unknown_uvw[2] - 1];//
				w_temp0[ii] = Xstru[element[j].node[ii].unknown_uvw[2] - 1];//
			}
			else {
				w_temp0[ii] = 0;
			}
			//double b_temp = element[j].b[0] / (6 * element[j].Ve) + element[j].b[1] / (6 * element[j].Ve) + element[j].b[2] / (6 * element[j].Ve) + element[j].b[3] / (6 * element[j].Ve);
			//double c_temp = element[j].c[0] / (6 * element[j].Ve) + element[j].c[1] / (6 * element[j].Ve) + element[j].c[2] / (6 * element[j].Ve) + element[j].c[3] / (6 * element[j].Ve);
			//double d_temp = element[j].d[0] / (6 * element[j].Ve) + element[j].d[1] / (6 * element[j].Ve) + element[j].d[2] / (6 * element[j].Ve) + element[j].d[3] / (6 * element[j].Ve);
			////sigma_six[0] = (element[j].b[ii] / (6 * element[j].Ve) * u_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (element[j].c[ii] / (6 * element[j].Ve) * v_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (element[j].d[ii] / (6 * element[j].Ve) * w_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
			////sigma_six[1] = (element[j].b[ii] / (6 * element[j].Ve) * u_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (element[j].c[ii] / (6 * element[j].Ve) * v_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (element[j].d[ii] / (6 * element[j].Ve) * w_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
			////sigma_six[2] = (element[j].b[ii] / (6 * element[j].Ve) * u_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (element[j].c[ii] / (6 * element[j].Ve) * v_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (element[j].d[ii] / (6 * element[j].Ve) * w_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]);//
			////sigma_six[3] = element[j].c[ii] / (6 * element[j].Ve) * u_temp * (mu[j]) + element[j].b[ii] / (6 * element[j].Ve) * v_temp * (mu[j]);//
			////sigma_six[4] = element[j].d[ii] / (6 * element[j].Ve) * v_temp * (mu[j]) + element[j].c[ii] / (6 * element[j].Ve) * w_temp * (mu[j]);//
			////sigma_six[5] = element[j].d[ii] / (6 * element[j].Ve) * u_temp * (mu[j]) + element[j].b[ii] / (6 * element[j].Ve) * w_temp * (mu[j]);//
			//sigma_six[0] = (b_temp * u_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (c_temp * v_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (d_temp * w_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
			//sigma_six[1] = (b_temp * u_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (c_temp * v_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (d_temp * w_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
			//sigma_six[2] = (b_temp * u_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (c_temp * v_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (d_temp * w_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]);//
			//sigma_six[3] = c_temp * u_temp * (mu[j]) + b_temp * v_temp * (mu[j]);//
			//sigma_six[4] = d_temp * v_temp * (mu[j]) + c_temp * w_temp * (mu[j]);//
			//sigma_six[5] = d_temp * u_temp * (mu[j]) + b_temp * w_temp * (mu[j]);//
			//double sigma_temp = sqrt(1.0 / 2 * ((sigma_six[0] - sigma_six[1]) * (sigma_six[0] - sigma_six[2]) + (sigma_six[2] - sigma_six[1]) * (sigma_six[0] - sigma_six[1]) + (sigma_six[0] - sigma_six[1]) * (sigma_six[0] - sigma_six[1]) + \
			//	3*(sigma_six[3]* sigma_six[3]+ sigma_six[4] * sigma_six[4]+ sigma_six[5] * sigma_six[5])));//
			//element[j].node[ii].von_Mises = sigma_temp;//
			////von_Mises_temp[element[j].node[ii].ver - 1] = sigma_temp;
			//if (element[j].node[ii].unknown_uvw[0] != 0) {
			//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[0] - 1] = sigma_temp;
			//}
			//if (element[j].node[ii].unknown_uvw[1] != 0) {
			//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[1] - 1] = sigma_temp;
			//}
			//if (element[j].node[ii].unknown_uvw[2] != 0) {
			//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[2] - 1] = sigma_temp;
			//}
		}

		double b_temp = element[j].b[0] / (6 * element[j].Ve) * u_temp0[0] + element[j].b[1] / (6 * element[j].Ve) * u_temp0[1] + element[j].b[2] / (6 * element[j].Ve) * u_temp0[2] + element[j].b[3] / (6 * element[j].Ve) * u_temp0[3];
		double c_temp = element[j].c[0] / (6 * element[j].Ve) * v_temp0[0] + element[j].c[1] / (6 * element[j].Ve) * v_temp0[1] + element[j].c[2] / (6 * element[j].Ve) * v_temp0[2] + element[j].c[3] / (6 * element[j].Ve) * v_temp0[3];
		double d_temp = element[j].d[0] / (6 * element[j].Ve) * w_temp0[0] + element[j].d[1] / (6 * element[j].Ve) * w_temp0[1] + element[j].d[2] / (6 * element[j].Ve) * w_temp0[2] + element[j].d[3] / (6 * element[j].Ve) * w_temp0[3];
		sigma_six[0] = (b_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (c_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (d_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
		sigma_six[1] = (b_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (c_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]) + (d_temp - element[j].alpha * (TT - T0)) * (lambda[j]);//
		sigma_six[2] = (b_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (c_temp - element[j].alpha * (TT - T0)) * (lambda[j]) + (d_temp - element[j].alpha * (TT - T0)) * (lambda[j] + 2 * mu[j]);//

		double b_temp0 = element[j].b[0] / (6 * element[j].Ve) * v_temp0[0] + element[j].b[1] / (6 * element[j].Ve) * v_temp0[1] + element[j].b[2] / (6 * element[j].Ve) * v_temp0[2] + element[j].b[3] / (6 * element[j].Ve) * v_temp0[3];
		double c_temp0 = element[j].c[0] / (6 * element[j].Ve) * w_temp0[0] + element[j].c[1] / (6 * element[j].Ve) * w_temp0[1] + element[j].c[2] / (6 * element[j].Ve) * w_temp0[2] + element[j].c[3] / (6 * element[j].Ve) * w_temp0[3];
		double d_temp0 = element[j].d[0] / (6 * element[j].Ve) * u_temp0[0] + element[j].d[1] / (6 * element[j].Ve) * u_temp0[1] + element[j].d[2] / (6 * element[j].Ve) * u_temp0[2] + element[j].d[3] / (6 * element[j].Ve) * u_temp0[3];

		double b_temp1 = element[j].b[0] / (6 * element[j].Ve) * w_temp0[0] + element[j].b[1] / (6 * element[j].Ve) * w_temp0[1] + element[j].b[2] / (6 * element[j].Ve) * w_temp0[2] + element[j].b[3] / (6 * element[j].Ve) * w_temp0[3];
		double c_temp1 = element[j].c[0] / (6 * element[j].Ve) * u_temp0[0] + element[j].c[1] / (6 * element[j].Ve) * u_temp0[1] + element[j].c[2] / (6 * element[j].Ve) * u_temp0[2] + element[j].c[3] / (6 * element[j].Ve) * u_temp0[3];
		double d_temp1 = element[j].d[0] / (6 * element[j].Ve) * v_temp0[0] + element[j].d[1] / (6 * element[j].Ve) * v_temp0[1] + element[j].d[2] / (6 * element[j].Ve) * v_temp0[2] + element[j].d[3] / (6 * element[j].Ve) * v_temp0[3];

		sigma_six[3] = c_temp1 * (mu[j]) + b_temp0 * (mu[j]);//
		sigma_six[4] = d_temp1 * (mu[j]) + c_temp0 * (mu[j]);//
		sigma_six[5] = d_temp0 * (mu[j]) + b_temp1 * (mu[j]);//


		double sigma_temp = sqrt(1.0 / 2 * ((sigma_six[0] - sigma_six[1]) * (sigma_six[0] - sigma_six[1]) + (sigma_six[2] - sigma_six[0]) * (sigma_six[2] - sigma_six[0]) + (sigma_six[2] - sigma_six[1]) * (sigma_six[2] - sigma_six[1]) + \
			3 * (sigma_six[3] * sigma_six[3] + sigma_six[4] * sigma_six[4] + sigma_six[5] * sigma_six[5])));//
		if (element[j].Material == 3 || element[j].Material == 1) {
			for (int ii = 0; ii < 4; ii++) {
				element[j].node[ii].von_Mises = sigma_temp;//
				//von_Mises_temp[element[j].node[ii].ver - 1] = sigma_temp;
				//if (element[j].node[ii].unknown_uvw[0] != 0) {
				//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[0] - 1] = sigma_temp;
				//}
				//if (element[j].node[ii].unknown_uvw[1] != 0) {
				//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[1] - 1] = sigma_temp;
				//}
				//if (element[j].node[ii].unknown_uvw[2] != 0) {
				//	Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[2] - 1] = sigma_temp;
				//}
			}
		}
		for (int ii = 0; ii < 4; ii++) {
			element[j].node[ii].von_Mises = sigma_temp;//
			//von_Mises_temp[element[j].node[ii].ver - 1] = sigma_temp;
			if (element[j].node[ii].unknown_uvw[0] != 0) {
				Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[0] - 1] = sigma_temp;
			}
			if (element[j].node[ii].unknown_uvw[1] != 0) {
				Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[1] - 1] = sigma_temp;
			}
			if (element[j].node[ii].unknown_uvw[2] != 0) {
				Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[2] - 1] = sigma_temp;
			}
		}
		//for (int ii = 0; ii < 4; ii++) {
		//	element[j].node[ii].von_Mises = sigma_temp;//
		//	//von_Mises_temp[element[j].node[ii].ver - 1] = sigma_temp;
		//	if (element[j].node[ii].unknown_uvw[0] != 0) {
		//		Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[0] - 1] = sigma_temp;
		//	}
		//	if (element[j].node[ii].unknown_uvw[1] != 0) {
		//		Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[1] - 1] = sigma_temp;
		//	}
		//	if (element[j].node[ii].unknown_uvw[2] != 0) {
		//		Xstruvon_Mises_temp[element[j].node[ii].unknown_uvw[2] - 1] = sigma_temp;
		//	}
		//}
	}

	//
	double* X_von_Mises_temp = new double[num_node];
	for (int i = 0; i < num_node; i++) {
		X_von_Mises_temp[i] = 0;
	}
	for (int j = 0; j < num_element_subdomain; j++) {
		if (element[j].Material == 3 || element[j].Material == 1) {
			for (int ii = 0; ii < 4; ii++) {
				X_von_Mises_temp[element[j].node[ii].ver - 1] = element[j].node[ii].von_Mises;
			}
		}

		//for (int ii = 0; ii < 4; ii++) {
		//	X_von_Mises_temp[element[j].node[ii].ver - 1] = element[j].node[ii].von_Mises;
		//}
	}

	int un_F_temp = 0;
	int* Accumulated_F_temp = new int[num_domain];
	int* unknown_F_temp = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		Accumulated_F_temp[j] = 0;
		unknown_F_temp[j] = num_node;

	}

	for (int j = 1; j < num_domain; j++) {
		Accumulated_F_temp[j] = Accumulated_F_temp[j - 1] + num_node;
	}
	double* X_F = new double[num_node * num_domain];
	MPI_Gatherv(X_von_Mises_temp, num_node, MPI_DOUBLE, X_F, unknown_F_temp, Accumulated_F_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete[]X_von_Mises_temp;

	double* X_von_Mises_temp2 = new double[num_node];

	if (myid == 0) {
		for (int i = 0; i < num_node * num_domain; i++) {
			int index_temp = i / num_node;
			//if (X_F[i] != 0) {
			if (X_F[i] > 100) {
				X_von_Mises_temp2[index_temp] = X_F[i];
			}
		}
		ofstream ofs12161("Patch_von_Mises1.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12161 << X_von_Mises_temp2[i] << endl;
		}
	}


	//



	for (int j = 0; j < num_unKnown_subdomain_stru[myid][0]; j++) {
		X_h_temp[j] = Xstru[j];
	}

	MPI_Gatherv(Xstruvon_Mises_temp, num_unKnown_subdomain_stru[myid][0], MPI_DOUBLE, von_Mises_temp, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(X_h_temp, num_unKnown_subdomain_stru[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_u, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_u, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_v, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_v, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(node_global_temp_w, num_unKnown_subdomain_stru[myid][0], MPI_INT, V_T_w, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);


	delete[]Accumulated_T_temp; delete[]unknown_T_temp;  delete[]node_global_temp_u; delete[]node_global_temp_v; delete[]node_global_temp_w; delete[]X_h_temp; delete[]Xstru;




	double* XX1_u = new double[num_node];
	double* XX1_v = new double[num_node];
	double* XX1_w = new double[num_node];
	double* Xstruvon_Mises_O = new double[num_node];//

	for (int i = 0; i < num_node; i++) {
		XX1_u[i] = 0;
		XX1_v[i] = 0;
		XX1_w[i] = 0;
	}


	if (myid == 0) {
		for (int i = 0; i < un_T_temp; i++) {
			int temp_ver = V_T_u[i];
			if (temp_ver != 0) {
				XX1_u[temp_ver - 1] = X_T[i];
				Xstruvon_Mises_O[temp_ver - 1] = von_Mises_temp[i];//
			}
			int temp_ver1 = V_T_v[i];
			if (temp_ver1 != 0) {
				XX1_v[temp_ver1 - 1] = X_T[i];
			}
			int temp_ver2 = V_T_w[i];
			if (temp_ver2 != 0) {
				XX1_w[temp_ver2 - 1] = X_T[i];
			}

		}

		ofstream ofs12157("Patch_u.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12157 << XX1_u[i] << endl;
		}
		ofstream ofs12158("Patch_v.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12158 << XX1_v[i] << endl;
		}
		ofstream ofs12159("Patch_w.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12159 << XX1_w[i] << endl;
		}
		ofstream ofs12160("Patch_von_Mises.txt");
		for (int i = 0; i < num_node; i++) {
			ofs12160 << Xstruvon_Mises_O[i] << endl;
		}


	}


	delete[]X_T; delete[]XX1_u; delete[]V_T_u; delete[]V_T_v; delete[]V_T_w; delete[]von_Mises_temp;//




	//	vector<Triplet<double>> Tri_tmp(num_unKnown_stb);
	//	for (int i = 0; i < num_unKnown_stb; i++) {
	//		Tri_tmp[i] = Triplet<double>(i, i, 1);
	//	}
	//	TriList_stru.insert(TriList_stru.end(), Tri_tmp.begin(), Tri_tmp.end());
	//	double* bb, * x;
	//	int cnt1 = 0, cnt2 = 0, tot_dm;
	//	MKL_INT maxfct, mnum, mtype, phase1, phase2, phase3, nrhs, msglvl, error;
	//	MKL_INT iparm[64]{ 0 };
	//	int pt[64]{ 0 };
	//	MKL_INT* perm2 = nullptr;
	//	maxfct = 1; mnum = 1; mtype = 11; phase1 = 13; phase2 = -1; msglvl = 0;
	//	iparm[0] = 1; iparm[1] = 2; iparm[9] = 13; iparm[10] = 1; iparm[12] = 1; iparm[17] = -1; iparm[18] = -1; iparm[24] = 1;
	//	nzero1 = 0; row1 = 0;
	//	int cnt_n1, cnt_n2, start1, start2, end1, end2, final_row = 0, final_col = 0;
	//	xx = new double[num_unKnown_stb];
	//	int size = TriList_stru.size();
	//	double* m_final = new double[size];
	//	double* m_pardiso1 = new double[size];
	//	int* rowm = new int[size]; int* colm = new int[size]; int* im1 = new int[num_unKnown_stb + 1]; int* jm1 = new int[size];
	//	nnz = size;
	//#	pragma omp parallel for num_threads(96)
	//	for (int i = 0; i < size; i++) {
	//		m_final[i] = TriList_stru[i].value();
	//		rowm[i] = TriList_stru[i].row() + 1;
	//		colm[i] = TriList_stru[i].col() + 1;
	//	}
	//	TriList_stru.clear();
	//	nrhs = 1;
	//	int job[8] = { 2, 1, 1, 0, 0, 0, 0, 0 };
	//	job[4] = nnz;
	//	mkl_dcsrcoo(job, &num_unKnown_stb, m_pardiso1, jm1, im1, &nnz, m_final, rowm, colm, &info);
	//	delete[] m_final;
	//	delete[] rowm;
	//	delete[] colm;
	//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unKnown_stb, m_pardiso1, im1, jm1, perm2, &nrhs, iparm, &msglvl, fbb_stru, xx, &error);
	//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unKnown_stb, m_pardiso1, im1, jm1, perm2, &nrhs, iparm, &msglvl, fbb_stru, xx, &error);
	//	delete[] m_pardiso1;
	//	delete[] im1;
	//	delete[] jm1;
	//	pthread_t tids[NUM_THREADS_stru];
	//	for (int i = 0; i < NUM_THREADS_stru; ++i) {
	//		int ret = pthread_create(&tids[i], NULL, Solver_stru, (void*)i);
	//	}
	//	for (int i = 0; i < NUM_THREADS_stru; i++) {
	//		pthread_join(tids[i], NULL);
	//	}
	//	delete[] xx;

}

int Matrix_Generator_stru(Element* element, Element* oppo_element, int myid) {
	unordered_map<long long, double> mat_full2;
	ofstream ofs;
	cout.precision(16);
	int material, ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	//double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;
	Vector3d zb1, zb2, zb3, zb4, zbc;
	//Vector4d Phi, T;
	long long temp;
	//fbr_stru.resize(num_unKnown_stru);
	//fbr_stru.setZero();

	fbr_stru = new double[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]];
	for (int i = 0; i < num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]; i++) {
		fbr_stru[i] = 0.0;
	}

	lambda = new double[num_element_subdomain];
	mu = new double[num_element_subdomain];
	for (el = 0; el < num_element_subdomain; el++) {
		lambda[el] = element[el].nu * element[el].E / (1 + element[el].nu) / (1 - 2 * element[el].nu);
		mu[el] = element[el].E / 2.0 / (1 + element[el].nu);

	}

	for (el = 0; el < num_element_subdomain; el++) {

		for (ii = 0; ii < 4; ii++) {
			//node_ii_glb = ver(el, ii);
			//ai = a(el, ii); bi = b(el, ii); ci = c(el, ii); di = d(el, ii);
			ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

			for (jj = 0; jj < 4; jj++) {
				//node_jj_glb = ver(el, jj);
				//aj = a(el, jj); bj = b(el, jj); cj = c(el, jj); dj = d(el, jj);
				aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];
				//node_ii_loc = ver(el, ii + 8);
				//node_jj_loc = ver(el, jj + 8);
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];

				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * bi * bj + mu[el] * ci * cj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * cj + mu[el] * ci * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * dj + mu[el] * di * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * bj + mu[el] * bi * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * ci * cj + mu[el] * bi * bj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * dj + mu[el] * di * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * bj + mu[el] * bi * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * cj + mu[el] * ci * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * di * dj + mu[el] * bi * bj + mu[el] * ci * cj) / (36.0 * element[el].Ve);
				}
			}
			if (element[el].node[ii].unknown_uvw[0] != 0 && element[el].node[ii].unknown_uvw[1] != 0 && element[el].node[ii].unknown_uvw[2] != 0) {
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				double TT = 0;
				for (int node_ii = 0; node_ii < 4; node_ii++) {
					int unknown_TT = element[el].node[node_ii].unknown_T;
					if (unknown_TT == 0) {
						TT += Td;
					}
					else {
						TT += Xh[unknown_TT - 1];
					}
				}
				//if (el < 3) { 
				//	cout << TT / 4.0 << endl;
				//}
				TT = TT / 4.0;
				if (TT < 293.15) {
					//cout <<"error TT = "<< TT << endl;
					//TT = 293.15;

				}
				//TT = 310.0;
				//fbr_stru[ver(el, ii + 8) - 1]  += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				//fbr_stru[ver(el, ii + 12) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				//fbr_stru[ver(el, ii + 16) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[0] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[1] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[2] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
			}
		}
		//for (nn = 0; nn < 4; nn++) {
		//	ofn1 = opp_stru[el][nn + 4];
		//	if (ofn1 == -71) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -72) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += -Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += -3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -73) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 12);
		//			//fbr_stru[node_ii_loc - 1] += Ffeeder / Afeeder * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 1.6e0 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//}
	}
	//cout << "tt1" << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (nn = 0; nn < 4; nn++) {
			//op1 = opp_stru[el][nn]; ofn1 = opp_stru[el][nn + 4];
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];

			if (op1 != 0) {
				//ndm_op = domain[op1 - 1];

				//if (ndm != ndm_op) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					//if (count == 10) {
					//	cout << "yess" << endl;
					//}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}

					//double factor = (E[el] * (1 - nu[el]) / ((1 + nu[el]) * (1 - 2 * nu[el])) + E[op1 - 1] * (1 - nu[op1 - 1]) / ((1 + nu[op1 - 1]) * (1 - 2 * nu[op1 - 1]))) / 2.0;
					double lambda_oppo = oppo_element[count_oppo_element - 1].nu * oppo_element[count_oppo_element - 1].E / (1 + oppo_element[count_oppo_element - 1].nu) / (1 - 2 * oppo_element[count_oppo_element - 1].nu);
					double mu_oppo = oppo_element[count_oppo_element - 1].E / 2.0 / (1 + oppo_element[count_oppo_element - 1].nu);
					double factor_temp = (element[el].E * (1 - element[el].nu) / ((1 + element[el].nu) * (1 - 2 * element[el].nu)) + oppo_element[count_oppo_element - 1].E * (1 - oppo_element[count_oppo_element - 1].nu) / ((1 + oppo_element[count_oppo_element - 1].nu) * (1 - 2 * oppo_element[count_oppo_element - 1].nu))) / 2.0;


					double factor = (lambda[el] + 2 * mu[el] + lambda_oppo + 2 * mu_oppo) / 2.0;
					double factor1 = factor;
					double factor2 = factor;

					//double factor = (lambda[el] + 2 * mu[el] + lambda_oppo + 2 * mu_oppo) / 2.0;
					//double factor1 = lambda[el] + 2 * mu[el];
					//double factor2 = lambda_oppo + 2 * mu_oppo;


					if (count == 10) {
						cout << "factor_temp = " << factor_temp << endl;
						cout << "factor = " << factor << endl;
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);
						for (jj = 0; jj < 3; jj++) {
							//node_jj_glb = ver(el, face_node[nn][jj] - 1);
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							//node_jj1_loc = ver(el, face_node[nn][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm - 1);
							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;//T T
							//mat_full2[temp] -= 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));


							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							//mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}


							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += 0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));



							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							//mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][1];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}



							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc= Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							//mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(el).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							//mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
								mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}


						}
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);

						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;

						for (jj = 0; jj < 3; jj++) {

							//node_jj_glb = ver(op1 - 1, face_node[ofn1 - 1][jj] - 1);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm_op - 1);
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							//mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 11);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							if (node_ii1_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += -0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							//mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							if (node_ii1_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}


							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += -1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							//mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							if (node_ii2_loc != 0 && node_jj1_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(op1 - 1).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							//temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							//mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							if (node_ii2_loc != 0 && node_jj2_loc != 0) {
								temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
								mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							}
						}
					}
				}
			}
		}
	}

	//cout << "tt2" << endl;
	int matsize = mat_full2.size();
	iRow_stru = new int[matsize];
	jCol_stru = new int[matsize];
	double* Me_full2 = new double[matsize];
	int nnz = 0;
	for (auto full : mat_full2) {
		Me_full2[nnz] = full.second;
		iRow_stru[nnz] = full.first / num_unKnown_stru + 1;
		jCol_stru[nnz++] = full.first % num_unKnown_stru + 1;
	}
	cout << nnz << endl;

	num_nzero_Pmatrix_stru = new int[num_domain];
	acoo1_stru = new double[nnz];
	rowind1_stru = new int[nnz];
	colind1_stru = new int[nnz];
	m_stru = new double[nnz];
	mrow_stru = new int[nnz];
	mcol_stru = new int[nnz];
	int nn1 = 0, nn2 = 0, nn3 = 0, mm1, mm2, mm3, cc = 0;;
	nnz_dm = 0; nnz_tot = 0;

	nnz_C_stru = new int* [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_C_stru[i] = new int[num_domain];

	for (int i = 0; i < num_domain; ++i) {
		for (int j = 0; j < num_domain; ++j) {
			nnz_C_stru[i][j] = 0;
		}
	}


	for (int n1 = myid; n1 < myid + 1; n1++) {
		mm1 = nn1 + num_unKnown_subdomain_stru[n1][0] + num_unKnown_subdomain_stru[n1][1];
		nn2 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2 = nn2 + num_unKnown_subdomain_stru[n2][0] + num_unKnown_subdomain_stru[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow_stru[j] > nn1 && iRow_stru[j] <= mm1 && jCol_stru[j] > nn2 && jCol_stru[j] <= mm2) {
#						pragma omp critical
					{
						++nnz_dm;
						m_stru[cc] = Me_full2[j];
						mrow_stru[cc] = iRow_stru[j] - nn1;
						mcol_stru[cc] = jCol_stru[j] - nn2;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo1_stru[nnz_tot] = Me_full2[j];
							rowind1_stru[nnz_tot] = iRow_stru[j] - nn1;
							colind1_stru[nnz_tot] = jCol_stru[j] - nn2;
							++nnz_tot;
						}
					}
				}
			}
			nn2 = mm2;
			nnz_C_stru[n1][n2] = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_stru[n1] = nnz_dm;
		}
		nn1 = mm1;
	}
	delete[] Me_full2;
	delete[] iRow_stru;
	delete[] jCol_stru;
	p_pardiso_stru = new double[nnz_tot];
	ip_stru = new int[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1] + 1];
	jp_stru = new int[nnz_tot];
	//cout << "tt3" << endl;

	r_dm_stru = new int[num_domain];
	r_dm_stru[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_stru[i] = r_dm_stru[i - 1] + Node_Boundary_stru[(i - 1)];
	}
	cout << " my 123id is " << myid << endl;
	int myid1 = myid;

	Direct_solution_stru(element, myid1);
	//Result_Out1();
	return 0;
}

int Matrix_Generator_stru2(Element* element, Element* oppo_element, int myid) {
	unordered_map<long long, double> mat_full2;
	ofstream ofs;
	cout.precision(16);
	int material, ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	//double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;
	Vector3d zb1, zb2, zb3, zb4, zbc;
	//Vector4d Phi, T;
	long long temp;
	//fbr_stru.resize(num_unKnown_stru);
	//fbr_stru.setZero();

	fbr_stru = new double[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]];
	for (int i = 0; i < num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]; i++) {
		fbr_stru[i] = 0.0;
	}

	lambda = new double[num_element_subdomain];
	mu = new double[num_element_subdomain];
	for (el = 0; el < num_element_subdomain; el++) {
		lambda[el] = element[el].nu * element[el].E / (1 + element[el].nu) / (1 - 2 * element[el].nu);
		mu[el] = element[el].E / 2.0 / (1 + element[el].nu);

	}

	for (el = 0; el < num_element_subdomain; el++) {

		for (ii = 0; ii < 4; ii++) {
			//node_ii_glb = ver(el, ii);
			//ai = a(el, ii); bi = b(el, ii); ci = c(el, ii); di = d(el, ii);
			ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

			for (jj = 0; jj < 4; jj++) {
				//node_jj_glb = ver(el, jj);
				//aj = a(el, jj); bj = b(el, jj); cj = c(el, jj); dj = d(el, jj);
				aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];
				//node_ii_loc = ver(el, ii + 8);
				//node_jj_loc = ver(el, jj + 8);
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];

				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * bi * bj + mu[el] * ci * cj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * cj + mu[el] * ci * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * dj + mu[el] * di * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * bj + mu[el] * bi * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * ci * cj + mu[el] * bi * bj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * dj + mu[el] * di * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * bj + mu[el] * bi * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * cj + mu[el] * ci * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * di * dj + mu[el] * bi * bj + mu[el] * ci * cj) / (36.0 * element[el].Ve);
				}
			}
			if (element[el].node[ii].unknown_uvw[0] != 0 && element[el].node[ii].unknown_uvw[1] != 0 && element[el].node[ii].unknown_uvw[2] != 0) {
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				double TT = 0;
				for (int node_ii = 0; node_ii < 4; node_ii++) {
					int unknown_TT = element[el].node[node_ii].unknown_T;
					if (unknown_TT == 0) {
						TT += T0;
					}
					else {
						TT += Xh[unknown_TT - 1];
					}
				}
				//if (el < 3) { 
				//	cout << TT / 4.0 << endl;
				//}
				TT = TT / 4.0;
				if (TT < 293.15) {
					//cout <<"error TT = "<< TT << endl;
					//TT = 293.15;

				}
				//TT = 310.0;
				//fbr_stru[ver(el, ii + 8) - 1]  += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				//fbr_stru[ver(el, ii + 12) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				//fbr_stru[ver(el, ii + 16) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[0] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[1] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[2] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
			}
		}
		//for (nn = 0; nn < 4; nn++) {
		//	ofn1 = opp_stru[el][nn + 4];
		//	if (ofn1 == -71) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -72) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += -Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += -3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -73) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 12);
		//			//fbr_stru[node_ii_loc - 1] += Ffeeder / Afeeder * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 1.6e0 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//}
	}
	//cout << "tt1" << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (nn = 0; nn < 4; nn++) {
			//op1 = opp_stru[el][nn]; ofn1 = opp_stru[el][nn + 4];
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];

			if (op1 != 0) {
				//ndm_op = domain[op1 - 1];

				//if (ndm != ndm_op) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					//if (count == 10) {
					//	cout << "yess" << endl;
					//}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}

					//double factor = (E[el] * (1 - nu[el]) / ((1 + nu[el]) * (1 - 2 * nu[el])) + E[op1 - 1] * (1 - nu[op1 - 1]) / ((1 + nu[op1 - 1]) * (1 - 2 * nu[op1 - 1]))) / 2.0;
					double lambda_oppo = oppo_element[count_oppo_element - 1].nu * oppo_element[count_oppo_element - 1].E / (1 + oppo_element[count_oppo_element - 1].nu) / (1 - 2 * oppo_element[count_oppo_element - 1].nu);
					double mu_oppo = oppo_element[count_oppo_element - 1].E / 2.0 / (1 + oppo_element[count_oppo_element - 1].nu);
					double factor_temp = (element[el].E * (1 - element[el].nu) / ((1 + element[el].nu) * (1 - 2 * element[el].nu)) + oppo_element[count_oppo_element - 1].E * (1 - oppo_element[count_oppo_element - 1].nu) / ((1 + oppo_element[count_oppo_element - 1].nu) * (1 - 2 * oppo_element[count_oppo_element - 1].nu))) / 2.0;


					//double factor = (lambda[el] + 2 * mu[el] + lambda_oppo + 2 * mu_oppo) / 2.0;
					//double factor1 = factor;
					//double factor2 = factor;

					double factor = (lambda[el] + 2 * mu[el] + lambda_oppo + 2 * mu_oppo) / 2.0;
					double factor1 = lambda[el] + 2 * mu[el];
					double factor2 = lambda_oppo + 2 * mu_oppo;


					if (count == 10) {
						cout << "factor_temp = " << factor_temp << endl;
						cout << "factor = " << factor << endl;
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);
						for (jj = 0; jj < 3; jj++) {
							//node_jj_glb = ver(el, face_node[nn][jj] - 1);
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							//node_jj1_loc = ver(el, face_node[nn][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm - 1);
							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;//T T
							//mat_full2[temp] -= 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));


							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += 0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));



							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][1];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));



							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc= Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(el).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


						}
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);

						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;

						for (jj = 0; jj < 3; jj++) {

							//node_jj_glb = ver(op1 - 1, face_node[ofn1 - 1][jj] - 1);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm_op - 1);
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 11);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += -0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];

							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii1_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += -1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(op1 - 1).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii2_loc - 1) * num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
						}
					}
				}
			}
		}
	}

	//cout << "tt2" << endl;
	int matsize = mat_full2.size();
	iRow_stru = new int[matsize];
	jCol_stru = new int[matsize];
	double* Me_full2 = new double[matsize];
	int nnz = 0;
	for (auto full : mat_full2) {
		Me_full2[nnz] = full.second;
		iRow_stru[nnz] = full.first / num_unKnown_stru + 1;
		jCol_stru[nnz++] = full.first % num_unKnown_stru + 1;
	}
	cout << nnz << endl;

	num_nzero_Pmatrix_stru = new int[num_domain];
	acoo1_stru = new double[nnz];
	rowind1_stru = new int[nnz];
	colind1_stru = new int[nnz];
	m_stru = new double[nnz];
	mrow_stru = new int[nnz];
	mcol_stru = new int[nnz];
	int nn1 = 0, nn2 = 0, nn3 = 0, mm1, mm2, mm3, cc = 0;;
	nnz_dm = 0; nnz_tot = 0;

	nnz_C_stru = new int* [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_C_stru[i] = new int[num_domain];

	for (int i = 0; i < num_domain; ++i) {
		for (int j = 0; j < num_domain; ++j) {
			nnz_C_stru[i][j] = 0;
		}
	}


	for (int n1 = myid; n1 < myid + 1; n1++) {
		mm1 = nn1 + num_unKnown_subdomain_stru[n1][0] + num_unKnown_subdomain_stru[n1][1];
		nn2 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2 = nn2 + num_unKnown_subdomain_stru[n2][0] + num_unKnown_subdomain_stru[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow_stru[j] > nn1 && iRow_stru[j] <= mm1 && jCol_stru[j] > nn2 && jCol_stru[j] <= mm2) {
#						pragma omp critical
					{
						++nnz_dm;
						m_stru[cc] = Me_full2[j];
						mrow_stru[cc] = iRow_stru[j] - nn1;
						mcol_stru[cc] = jCol_stru[j] - nn2;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo1_stru[nnz_tot] = Me_full2[j];
							rowind1_stru[nnz_tot] = iRow_stru[j] - nn1;
							colind1_stru[nnz_tot] = jCol_stru[j] - nn2;
							++nnz_tot;
						}
					}
				}
			}
			nn2 = mm2;
			nnz_C_stru[n1][n2] = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_stru[n1] = nnz_dm;
		}
		nn1 = mm1;
	}
	delete[] Me_full2;
	delete[] iRow_stru;
	delete[] jCol_stru;
	p_pardiso_stru = new double[nnz_tot];
	ip_stru = new int[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1] + 1];
	jp_stru = new int[nnz_tot];
	//cout << "tt3" << endl;

	r_dm_stru = new int[num_domain];
	r_dm_stru[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_stru[i] = r_dm_stru[i - 1] + Node_Boundary_stru[(i - 1)];
	}
	Direct_solution_stru(element, myid);
	//Result_Out1();
	return 0;
}


int Matrix_Generator_stru_error(Element* element, Element* oppo_element, int myid) {
	unordered_map<long long, double> mat_full2;
	ofstream ofs;
	cout.precision(16);
	int material, ndm, el, ii, jj, node_ii, node_jj, nn, mm, op1, ofn1, node_ii1, node_jj1, ndm_op, node_ii_glb, node_ii1_loc, node_ii2_loc, node_jj_glb, node_jj1_loc, node_jj2_loc, node_jj_loc, node_ii_loc;
	int nodeNo, nzero, nzero1, nzero2, nnz_dm, nnz_tot, unknown_dm, row1, row2, pp, qq, flag1, flag2, kk1, kk2, node1, node2, node3, node4, nGauss, num_nnz_est, info;
	double ai, aj, bj, cj, dj, bi, ci, di, xc, yc, zc, weight, cnt, Ni, Nj, N1, N2, N3, N4, Resitivity, Esqu, k;
	//double Q, Q0, Q1, T1, T2, T3, T4, Phi1, Phi2, Phi3, Phi4, factor;
	Vector3d zb1, zb2, zb3, zb4, zbc;
	//Vector4d Phi, T;
	long long temp;
	//fbr_stru.resize(num_unKnown_stru);
	//fbr_stru.setZero();

	fbr_stru = new double[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1]];
	lambda = new double[num_element_subdomain];
	mu = new double[num_element_subdomain];
	for (el = 0; el < num_element_subdomain; el++) {
		lambda[el] = element[el].nu * element[el].E / (1 + element[el].nu) / (1 - 2 * element[el].nu);
		mu[el] = element[el].E / 2.0 / (1 + element[el].nu);

	}

	for (el = 0; el < num_element_subdomain; el++) {

		for (ii = 0; ii < 4; ii++) {
			//node_ii_glb = ver(el, ii);
			//ai = a(el, ii); bi = b(el, ii); ci = c(el, ii); di = d(el, ii);
			ai = element[el].a[ii]; bi = element[el].b[ii]; ci = element[el].c[ii]; di = element[el].d[ii];

			for (jj = 0; jj < 4; jj++) {
				//node_jj_glb = ver(el, jj);
				//aj = a(el, jj); bj = b(el, jj); cj = c(el, jj); dj = d(el, jj);
				aj = element[el].a[jj]; bj = element[el].b[jj]; cj = element[el].c[jj]; dj = element[el].d[jj];
				//node_ii_loc = ver(el, ii + 8);
				//node_jj_loc = ver(el, jj + 8);
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];

				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * bi * bj + mu[el] * ci * cj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * cj + mu[el] * ci * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[0];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * bi * dj + mu[el] * di * bj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * bj + mu[el] * bi * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * ci * cj + mu[el] * bi * bj + mu[el] * di * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[1];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * ci * dj + mu[el] * di * cj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[0];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * bj + mu[el] * bi * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[1];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += (lambda[el] * di * cj + mu[el] * ci * dj) / (36.0 * element[el].Ve);
				}
				node_ii_loc = element[el].node[ii].unknown_uvw[2];
				node_jj_loc = element[el].node[jj].unknown_uvw[2];
				if ((node_ii_loc != 0) && (node_jj_loc != 0)) {
					temp = (long long)(node_ii_loc - 1) * num_unKnown_stru + node_jj_loc - 1 + Accumulated_unknowns_stru[myid];
					mat_full2[temp] += ((lambda[el] + 2 * mu[el]) * di * dj + mu[el] * bi * bj + mu[el] * ci * cj) / (36.0 * element[el].Ve);
				}
			}
			if (element[el].node[ii].unknown_uvw[0] != 0 && element[el].node[ii].unknown_uvw[1] != 0 && element[el].node[ii].unknown_uvw[2] != 0) {
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				//double TT = (Xtemp[ver(el, 0) - 1] + Xtemp[ver(el, 1) - 1] + Xtemp[ver(el, 2) - 1] + Xtemp[ver(el, 3) - 1]) / 4.0;
				double TT = 0;
				for (int node_ii = 0; node_ii < 4; node_ii++) {
					int unknown_TT = element[el].node[node_ii].unknown_T;
					if (unknown_TT == 0) {
						TT += T0;
					}
					else {
						TT += Xh[unknown_TT - 1];
					}
				}
				//if (el < 3) { 
				//	cout << TT / 4.0 << endl;
				//}
				TT = TT / 4.0;
				if (TT < 293.15) {
					//cout <<"error TT = "<< TT << endl;
					//TT = 293.15;

				}
				//fbr_stru[ver(el, ii + 8) - 1]  += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				//fbr_stru[ver(el, ii + 12) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				//fbr_stru[ver(el, ii + 16) - 1] += alpha[el] * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[0] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].b[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[1] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].c[ii] / 6;
				fbr_stru[element[el].node[ii].unknown_uvw[2] - 1] += element[el].alpha * (TT - T0) * (3 * lambda[el] + 2 * mu[el]) * element[el].d[ii] / 6;
			}
		}
		//for (nn = 0; nn < 4; nn++) {
		//	ofn1 = opp_stru[el][nn + 4];
		//	if (ofn1 == -71) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -72) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 8);
		//			//fbr_stru[node_ii_loc - 1] += -Fscrew / Awasher * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += -3.0e1 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//	if (ofn1 == -73) {
		//		for (ii = 0; ii < 3; ii++) {
		//			node_ii = face_node[nn][ii] - 1; node_ii_loc = ver(el, node_ii + 12);
		//			//fbr_stru[node_ii_loc - 1] += Ffeeder / Afeeder * 1.0 / 3 * Area(el, nn);
		//			fbr_stru[node_ii_loc - 1] += 1.6e0 * 1.0 / 3 * Area(el, nn);
		//		}
		//	}
		//}
	}
	cout << "tt1" << endl;

	int count = 0;
	int count_oppo_element = 0;
	for (el = 0; el < num_element_subdomain; el++) {
		ndm = element[el].domain;

		for (nn = 0; nn < 4; nn++) {
			//op1 = opp_stru[el][nn]; ofn1 = opp_stru[el][nn + 4];
			op1 = element[el].face[nn].opp[0]; ofn1 = element[el].face[nn].opp[1];

			if (op1 != 0) {
				//ndm_op = domain[op1 - 1];

				//if (ndm != ndm_op) {
				if (element[el].face[nn].whether_boundary) {
					count++;
					//if (count == 10) {
					//	cout << "yess" << endl;
					//}
					if (oppo_element[count_oppo_element].Global_num == op1) {
						op1 = count_oppo_element + 1;
						count_oppo_element++;
						//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
					}
					else {
						cout << "wrong" << endl;
					}

					//double factor = (E[el] * (1 - nu[el]) / ((1 + nu[el]) * (1 - 2 * nu[el])) + E[op1 - 1] * (1 - nu[op1 - 1]) / ((1 + nu[op1 - 1]) * (1 - 2 * nu[op1 - 1]))) / 2.0;
					double lambda_oppo = oppo_element[count_oppo_element - 1].nu * oppo_element[count_oppo_element - 1].E / (1 + oppo_element[count_oppo_element - 1].nu) / (1 - 2 * oppo_element[count_oppo_element - 1].nu);
					double mu_oppo = oppo_element[count_oppo_element - 1].E / 2.0 / (1 + oppo_element[count_oppo_element - 1].nu);
					double factor_temp = (element[el].E * (1 - element[el].nu) / ((1 + element[el].nu) * (1 - 2 * element[el].nu)) + oppo_element[count_oppo_element - 1].E * (1 - oppo_element[count_oppo_element - 1].nu) / ((1 + oppo_element[count_oppo_element - 1].nu) * (1 - 2 * oppo_element[count_oppo_element - 1].nu))) / 2.0;


					double factor = (lambda[el] + 2 * mu[el] + lambda_oppo + 2 * mu_oppo) / 2.0;
					double factor1 = factor;
					double factor2 = factor;

					//double factor = (lambda[el] + 3 * mu[el] + lambda_oppo + 3 * mu_oppo) / 2.0;
					//double factor1 = lambda[el] + 3 * mu[el];
					//double factor2 = lambda_oppo + 3 * mu_oppo;


					if (count == 10) {
						cout << "factor_temp = " << factor_temp << endl;
						cout << "factor = " << factor << endl;
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);
						for (jj = 0; jj < 3; jj++) {
							//node_jj_glb = ver(el, face_node[nn][jj] - 1);
							node_jj_glb = element[el].node[face_node[nn][jj] - 1].ver;
							//node_jj1_loc = ver(el, face_node[nn][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm - 1);
							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;//T T
							//mat_full2[temp] -= 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));


							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];
							
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += 0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));



							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];

							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][1];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] -= (factor1) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));



							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc= Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(el, face_node[nn][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[0];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = element[el].node[face_node[nn][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(el).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[0];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = element[el].node[face_node[nn][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[myid];
							mat_full2[temp] += factor1 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


						}
					}
					for (ii = 0; ii < 3; ii++) {
						//node_ii_glb = ver(el, face_node[nn][ii] - 1);
						//node_ii1_loc = ver(el, face_node[nn][ii] + 3);
						//node_ii2_loc = Node_RTC_Thermal(node_ii_glb - 1, ndm - 1);

						node_ii_glb = element[el].node[face_node[nn][ii] - 1].ver;

						for (jj = 0; jj < 3; jj++) {

							//node_jj_glb = ver(op1 - 1, face_node[ofn1 - 1][jj] - 1);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 3);
							//node_jj2_loc = Node_RTC_Thermal(node_jj_glb - 1, ndm_op - 1);
							node_jj_glb = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].ver;



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += 4.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];

							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 11);
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 11);
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] -= (4.0) * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));



							//temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += -0.5 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));
							//node_ii1_loc = ver(el, face_node[nn][ii] + 7);
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];

							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii1_loc = element[el].node[face_node[nn][ii] - 1].unknown_uvw[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii1_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += (factor2) * 0.5 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));


							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj1_loc - 1;
							//mat_full2[temp] += -1.0 * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj1_loc = ver(op1 - 1, face_node[ofn1 - 1][jj] + 7);
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[0];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[1];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj1_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_uvw[2];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj1_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += -factor * 1.0 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));

							//temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_Thermal + node_jj2_loc - 1;
							//mat_full2[temp] += Kh.row(op1 - 1).mean() * Area(el, nn) / 12 * (1 + Flabel(node_jj_glb, node_ii_glb));

							//node_ii2_loc = Node_RTC_stru[node_ii_glb - 1][ndm - 1][0];
							//node_jj2_loc = Node_RTC_stru[node_jj_glb - 1][ndm_op - 1][0];
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[0];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[0];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[1];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[1];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
							node_ii2_loc = element[el].node[face_node[nn][ii] - 1].unknown_abc[2];
							node_jj2_loc = oppo_element[count_oppo_element - 1].node[face_node[ofn1 - 1][jj] - 1].unknown_abc[2];
							temp = (long long)(node_ii2_loc - 1) * (long long)num_unKnown_stru + node_jj2_loc - 1 + Accumulated_unknowns_stru[oppo_element[count_oppo_element - 1].domain - 1];
							mat_full2[temp] += factor2 * element[el].face[nn].Area / 12 * (1 + Flabel_stru(node_jj_glb, node_ii_glb));
						}
					}
				}
			}
		}
	}

	//cout << "tt2" << endl;
	int matsize = mat_full2.size();
	iRow_stru = new int[matsize];
	jCol_stru = new int[matsize];
	double* Me_full2 = new double[matsize];
	int nnz = 0;
	for (auto full : mat_full2) {
		Me_full2[nnz] = full.second;
		iRow_stru[nnz] = full.first / num_unKnown_stru + 1;
		jCol_stru[nnz++] = full.first % num_unKnown_stru + 1;
	}
	cout << nnz << endl;

	num_nzero_Pmatrix_stru = new int[num_domain];
	acoo1_stru = new double[nnz];
	rowind1_stru = new int[nnz];
	colind1_stru = new int[nnz];
	m_stru = new double[nnz];
	mrow_stru = new int[nnz];
	mcol_stru = new int[nnz];
	int nn1 = 0, nn2 = 0, nn3 = 0, mm1, mm2, mm3, cc = 0;;
	nnz_dm = 0; nnz_tot = 0;

	nnz_C_stru = new int * [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_C_stru[i] = new int[num_domain];

	for (int i = 0; i < num_domain; ++i) {
		for (int j = 0; j < num_domain; ++j) {
			nnz_C_stru[i][j] = 0;
		}
	}


	for (int n1 = myid; n1 < myid+1; n1++) {
		mm1 = nn1 + num_unKnown_subdomain_stru[n1][0] + num_unKnown_subdomain_stru[n1][1];
		nn2 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			nnz_dm = 0;
			mm2 = nn2 + num_unKnown_subdomain_stru[n2][0] + num_unKnown_subdomain_stru[n2][1];
#			pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow_stru[j] > nn1 && iRow_stru[j] <= mm1 && jCol_stru[j] > nn2 && jCol_stru[j] <= mm2) {
#						pragma omp critical
					{
						++nnz_dm;
						m_stru[cc] = Me_full2[j];
						mrow_stru[cc] = iRow_stru[j] - nn1;
						mcol_stru[cc] = jCol_stru[j] - nn2;
						++cc;
					}
					if (n1 == n2) {
#						pragma omp critical
						{
							acoo1_stru[nnz_tot] = Me_full2[j];
							rowind1_stru[nnz_tot] = iRow_stru[j] - nn1;
							colind1_stru[nnz_tot] = jCol_stru[j] - nn2;
							++nnz_tot;
						}
					}
				}
			}
			nn2 = mm2;
			nnz_C_stru[n1][n2] = nnz_dm;
			if (n1 == n2) num_nzero_Pmatrix_stru[n1] = nnz_dm;
		}
		nn1 = mm1;
	}
	delete[] Me_full2;
	delete[] iRow_stru;
	delete[] jCol_stru;
	p_pardiso_stru = new double[nnz_tot];
	ip_stru = new int[num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1] + 1];
	jp_stru = new int[nnz_tot];
	//cout << "tt3" << endl;

	r_dm_stru = new int[num_domain];
	r_dm_stru[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		r_dm_stru[i] = r_dm_stru[i - 1] + Node_Boundary_stru[(i - 1)];
	}
	Direct_solution_stru(element,myid);
	//Result_Out1();
	return 0;
}


int Result_Out1() {
	double* Xtemp_stru = new double[num_node];
	double* Xtemp_stru_V = new double[num_node];
	double* Xtemp_stru_W = new double[num_node];
	int el, ii, kk, nodeii;
	ofstream ofs2;
	ofs2.precision(16);
	for (int el = 0; el < num_element_subdomain; el++) {
		for (int ii = 0; ii < 4; ii++) {
			if (ver(el, ii + 8) == 0 && ver(el, ii + 12) == 0 && ver(el, ii + 16) == 0) {
				Xtemp_stru[ver(el, ii) - 1] = 0;
				Xtemp_stru_V[ver(el, ii) - 1] = 0;
				Xtemp_stru_W[ver(el, ii) - 1] = 0;
			}
			else {
				Xtemp_stru[ver(el, ii) - 1] = Xstru[ver(el, ii + 8) - 1] * 1e6;
				Xtemp_stru_V[ver(el, ii) - 1] = Xstru[ver(el, ii + 12) - 1] * 1e6;
				Xtemp_stru_W[ver(el, ii) - 1] = Xstru[ver(el, ii + 16) - 1] * 1e6;
			}
		}
	}
	string   T = "u_RTC.txt";
	string   T1 = "v_RTC.txt";
	string   T2 = "w_RTC.txt";
	ofs2.open(T, ios::out);
	for (kk = 0; kk < num_node; kk++) {
		ofs2 << Xtemp_stru[kk] << '\n';
	}
	ofs2.close();
	ofs2.open(T1, ios::out);
	for (kk = 0; kk < num_node; kk++) {
		ofs2 << Xtemp_stru_V[kk] << '\n';
	}
	ofs2.close();
	ofs2.open(T2, ios::out);
	for (kk = 0; kk < num_node; kk++) {
		ofs2 << Xtemp_stru_W[kk] << '\n';
	}
	ofs2.close();
	return 0;
}

double Flabel_stru(int vnode1, int vnode2) {
	double Flable0;
	if (vnode1 == vnode2) {
		Flable0 = 1.0;
	}
	else {
		Flable0 = 0.0;
	}
	return Flable0;
}