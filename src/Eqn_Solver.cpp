#include<ctime>
#include<algorithm>
#include "Element.h"
#include "pthread.h"
#include"Global_Data.h"
#include<Eigen>
#include<unordered_map>
#include<complex>
#include<fstream>
#include"Matrix_Module.h"
#include <stdlib.h>
#include<mpi.h>
#include<omp.h>
#include"mkl_scalapack.h"
using namespace std;
using namespace Eigen;
using Eigen::Triplet;

double sign_judge(int edgeNo1, int edgeNo2);
double Flabel(int vnode1, int vnode2);
long long  matsize = 0;
vector<Eigen::Triplet<complex<double>>> TriList;

void FETI_like_procedure_iter1(int myid) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = 0;


		/*        Columnar calculation      */
		//
		//cout << "test2" << endl;
		//cnt_n1 = Edge_Boundary[n1];
		//int num_col_solution = 100;
		//int num_solution;
		//num_solution = cnt_n1 / num_col_solution;
		//Mat_I = new complex<double>[unknown_dm * num_col_solution]();
		//Inv_M_A = (complex<double>*)malloc(unknown_dm * num_col_solution * sizeof(complex<double>));
		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		//nrhs = num_col_solution;
		//cout << "test1" << endl;
		//int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		//mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		//for (int i = 0; i < num_solution; i++) {
		//	if (i > 0) {
		//		for (int j = 0; j < num_col_solution; j++) {
		//			Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
		//		}
		//	}
		//	for (int j = 0; j < num_col_solution; j++) {
		//		Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		//	}

		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		//	for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
		//		for (int j = 0; j < num_col_solution; j++) {
		//			int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
		//			int InvAindex = k + j * unknown_dm;
		//			mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
		//		}
		//	}
		//}

		//delete[] Mat_I;
		//Mat_I = nullptr;
		//free(Inv_M_A);
		//int Remainder_column = cnt_n1 - num_col_solution * num_solution;
		//if (Remainder_column != 0) {
		//	Mat_I = new complex<double>[unknown_dm * Remainder_column]();
		//	Inv_M_A = (complex<double>*)malloc(unknown_dm * Remainder_column * sizeof(complex<double>));
		//	nrhs = Remainder_column;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		//	for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
		//		for (int j = 0; j < Remainder_column; j++) {
		//			int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + num_solution * num_col_solution;
		//			int InvAindex = k + j * unknown_dm;
		//			mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
		//		}
		//	}
		//	delete[] Mat_I;
		//	Mat_I = nullptr;
		//	free(Inv_M_A);
		//}

		/*    Columnar calculation      */




		cnt_n1 = Edge_Boundary[n1];
		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
		nrhs = cnt_n1;
		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
		for (int i = 0; i < cnt_n1; i++) {
			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		}
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);

		//for (int i = 0; i < num_domain; i++) {
		//	MPI_Barrier(MPI_COMM_WORLD);
		//	if (i == myid) {
		//		cout << "myid is " << myid << endl;
		//		cout << "start pardiso "  << endl;
		//		clock_t beginTime = clock();
		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//		if (error != 0)cout << "Pardiso error " << error << endl;
		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		//		clock_t EndTime = clock();
		//	}
		//}
		//cout << "start pardiso" << endl;
		clock_t beginTime = clock();
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		clock_t EndTime = clock();


		delete[] Mat_I;
		Mat_I = nullptr;

		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
#		pragma omp parallel for 
		for (int i = 0; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int InvAindex = i + j * unknown_dm;
				int ReIndex = cnt_n1 * i + j;
				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
			}
		}
		free(Inv_M_A);

		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
#		pragma omp parallel for  
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
				int InvAindex = cnt_n1 * i + j;
				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
			}
		}
		free(R_product_InvA);






		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		//complex<double>* yii = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
		//for (int i = 0; i < unknown_dm; i++) {
		//	yii[i].real(fbr[i].real);
		//	yii[i].imag(fbr[i].imag);
		//}
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure_iter(int myid) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error, phase3, phase4;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0; phase3 = 12; phase4 = 33;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = 0;


		/*        Columnar calculation      */

		cout << "test2" << endl;
		cnt_n1 = Edge_Boundary[n1];
		int num_col_solution = 200;
		cout << "num_col_solution = " << num_col_solution << endl;
		int num_solution;
		num_solution = cnt_n1 / num_col_solution;
		Mat_I = new complex<double>[unknown_dm * num_col_solution]();
		Inv_M_A = (complex<double>*)malloc(unknown_dm * num_col_solution * sizeof(complex<double>));
		long long cnt_n1_long = cnt_n1;
		mat_zm = new complex<double>[cnt_n1_long * cnt_n1_long];
		nrhs = num_col_solution;
		cout << "test1" << endl;
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);

		pardiso(pt, &maxfct, &mnum, &mtype, &phase3, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）

		for (int i = 0; i < num_solution; i++) {
			if (i > 0) {
				for (int j = 0; j < num_col_solution; j++) {
					Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
				}
			}
			for (int j = 0; j < num_col_solution; j++) {
				Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			}

			pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
				for (int j = 0; j < num_col_solution; j++) {
					long long long_temp1 = k - unknown_dm + cnt_n1;
					long long long_temp2 = j + i * num_col_solution;

					long long Zmindex = (long_temp1) * cnt_n1_long + long_temp2;
					int InvAindex = k + j * unknown_dm;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
		}
		cout << "test8" << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		delete[] Mat_I;
		Mat_I = nullptr;
		free(Inv_M_A);
		int Remainder_column = cnt_n1 - num_col_solution * num_solution;
		if (Remainder_column != 0) {
			Mat_I = new complex<double>[unknown_dm * Remainder_column]();
			for (int j = 0; j < Remainder_column; j++) {
				Mat_I[j * unknown_dm + (num_solution)*num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			}
			Inv_M_A = (complex<double>*)malloc(unknown_dm * Remainder_column * sizeof(complex<double>));
			nrhs = Remainder_column;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
			for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
				for (int j = 0; j < Remainder_column; j++) {

					long long long_temp1 = k - unknown_dm + cnt_n1;
					long long long_temp2 = j + num_solution * num_col_solution;

					long long Zmindex = (long_temp1) * cnt_n1_long + long_temp2;
					int InvAindex = k + j * unknown_dm;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
			delete[] Mat_I;
			Mat_I = nullptr;
			free(Inv_M_A);
		}

		/*    Columnar calculation      */




//		cnt_n1 = Edge_Boundary[n1];
//		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
//		nrhs = cnt_n1;
//		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
//		for (int i = 0; i < cnt_n1; i++) {
//			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
//		}
//		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
//		mkl_zcsrcoo(job, &unknown_dm, p_pardiso , jp , ip , &nnz_dm, acoo , rowind , colind , &info);
//
//		//for (int i = 0; i < num_domain; i++) {
//		//	MPI_Barrier(MPI_COMM_WORLD);
//		//	if (i == myid) {
//		//		cout << "myid is " << myid << endl;
//		//		cout << "start pardiso "  << endl;
//		//		clock_t beginTime = clock();
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		//		if (error != 0)cout << "Pardiso error " << error << endl;
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		//		clock_t EndTime = clock();
//		//	}
//		//}
//		//cout << "start pardiso" << endl;
//		clock_t beginTime = clock();
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		if (error != 0)cout << "Pardiso error " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		clock_t EndTime = clock();
//
//
//		delete[] Mat_I;
//		Mat_I = nullptr;
//
//		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
//#		pragma omp parallel for 
//		for (int i = 0; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int InvAindex = i + j * unknown_dm;
//				int ReIndex = cnt_n1 * i + j;
//				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
//			}
//		}
//		free(Inv_M_A);
//
//		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
//#		pragma omp parallel for  
//		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
//				int InvAindex = cnt_n1 * i + j;
//				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
//			}
//		}
//		free(R_product_InvA);






		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		//complex<double>* yii = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
		//for (int i = 0; i < unknown_dm; i++) {
		//	yii[i].real(fbr[i].real);
		//	yii[i].imag(fbr[i].imag);
		//}
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					long long cnt_n1_long = cnt_n1;
					long long cnt_n2_long = cnt_n2;
					long long matsize = cnt_n1_long * cnt_n2_long;
					//int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (long long i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							long long long_temp1 = mrow[i] - start1;
							long long long_temp2 = mcol[i] - start2;

							b_pro[cnt_n2_long * (long_temp1) + long_temp2].real(m[i].real);
							b_pro[cnt_n2_long * (long_temp1) + long_temp2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					//size_t List_step = 0;

					long long List_step = 0;
					for (long long i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure(int myid){
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid+1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1,  nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = r_dm[myid];


		/*        Columnar calculation      */
		//
		//cout << "test2" << endl;
		//cnt_n1 = Edge_Boundary[n1];
		//int num_col_solution = 100;
		//int num_solution;
		//num_solution = cnt_n1 / num_col_solution;
		//Mat_I = new complex<double>[unknown_dm * num_col_solution]();
		//Inv_M_A = (complex<double>*)malloc(unknown_dm * num_col_solution * sizeof(complex<double>));
		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		//nrhs = num_col_solution;
		//cout << "test1" << endl;
		//int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		//mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		//for (int i = 0; i < num_solution; i++) {
		//	if (i > 0) {
		//		for (int j = 0; j < num_col_solution; j++) {
		//			Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
		//		}
		//	}
		//	for (int j = 0; j < num_col_solution; j++) {
		//		Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		//	}

		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		//	for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
		//		for (int j = 0; j < num_col_solution; j++) {
		//			int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
		//			int InvAindex = k + j * unknown_dm;
		//			mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
		//		}
		//	}
		//}

		//delete[] Mat_I;
		//Mat_I = nullptr;
		//free(Inv_M_A);
		//int Remainder_column = cnt_n1 - num_col_solution * num_solution;
		//if (Remainder_column != 0) {
		//	Mat_I = new complex<double>[unknown_dm * Remainder_column]();
		//	Inv_M_A = (complex<double>*)malloc(unknown_dm * Remainder_column * sizeof(complex<double>));
		//	nrhs = Remainder_column;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		//	for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
		//		for (int j = 0; j < Remainder_column; j++) {
		//			int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + num_solution * num_col_solution;
		//			int InvAindex = k + j * unknown_dm;
		//			mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
		//		}
		//	}
		//	delete[] Mat_I;
		//	Mat_I = nullptr;
		//	free(Inv_M_A);
		//}

		/*    Columnar calculation      */




		cnt_n1 = Edge_Boundary[n1];
		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
		nrhs = cnt_n1;
		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
		for (int i = 0; i < cnt_n1; i++) {
			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		}
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso , jp , ip , &nnz_dm, acoo , rowind , colind , &info);

		//for (int i = 0; i < num_domain; i++) {
		//	MPI_Barrier(MPI_COMM_WORLD);
		//	if (i == myid) {
		//		cout << "myid is " << myid << endl;
		//		cout << "start pardiso "  << endl;
		//		clock_t beginTime = clock();
		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//		if (error != 0)cout << "Pardiso error " << error << endl;
		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		//		clock_t EndTime = clock();
		//	}
		//}
		//cout << "start pardiso" << endl;
		clock_t beginTime = clock();
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		clock_t EndTime = clock();


		delete[] Mat_I;
		Mat_I = nullptr;

		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
#		pragma omp parallel for 
		for (int i = 0; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int InvAindex = i + j * unknown_dm;
				int ReIndex = cnt_n1 * i + j;
				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
			}
		}
		free(Inv_M_A);

		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
#		pragma omp parallel for  
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
				int InvAindex = cnt_n1 * i + j;
				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
			}
		}
		free(R_product_InvA);
		





		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		//complex<double>* yii = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
		//for (int i = 0; i < unknown_dm; i++) {
		//	yii[i].real(fbr[i].real);
		//	yii[i].imag(fbr[i].imag);
		//}
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
					#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure_E_T(int myid) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;
	int startdm = myid;
	int enddm = myid + 1;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];  
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = r_dm[myid];
		cnt_n1 = Edge_Boundary[n1];
		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
		nrhs = cnt_n1;
		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
		for (int i = 0; i < cnt_n1; i++) {
			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		}
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		clock_t beginTime = clock();
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		clock_t EndTime = clock();
		delete[] Mat_I;
		Mat_I = nullptr;
		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
#		pragma omp parallel for 
		for (int i = 0; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int InvAindex = i + j * unknown_dm;
				int ReIndex = cnt_n1 * i + j;
				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
			}
		}
		free(Inv_M_A);

		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
#		pragma omp parallel for  
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
				int InvAindex = cnt_n1 * i + j;
				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
			}
		}
		free(R_product_InvA);
		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		free(pro);
		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;
					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}
		delete[] mat_zm;
	}
}

void FETI_like_procedure1(int myid) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error, phase3, phase4;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0; phase3 = 12, phase4 = 33;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = r_dm[myid];


		/*        Columnar calculation      */

		cout << "test2" << endl;
		cnt_n1 = Edge_Boundary[n1];
		int num_col_solution = 400;
		cout << "num_col_solution = " << num_col_solution << endl;
		int num_solution;
		num_solution = cnt_n1 / num_col_solution;
		Mat_I = new complex<double>[unknown_dm * num_col_solution]();
		Inv_M_A = (complex<double>*)malloc(unknown_dm * num_col_solution * sizeof(complex<double>));
		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		nrhs = num_col_solution;
		cout << "test1" << endl;
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase3, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		for (int i = 0; i < num_solution; i++) {
			time_t start, end;

			start = time(NULL);
			if (i > 0) {
				for (int j = 0; j < num_col_solution; j++) {
					Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
				}
			}
			for (int j = 0; j < num_col_solution; j++) {
				Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			}

			pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
				for (int j = 0; j < num_col_solution; j++) {
					int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
					int InvAindex = k + j * unknown_dm;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
			end = time(NULL);
			double time = (double)(end - start);
			if (myid == 0 && i == 0) {
				cout << "column solve time is" << time << endl;
			}



		}
		//cout << "test8" << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		delete[] Mat_I;
		Mat_I = nullptr;
		free(Inv_M_A);
		int Remainder_column = cnt_n1 - num_col_solution * num_solution;
		if (Remainder_column != 0) {
			Mat_I = new complex<double>[unknown_dm * Remainder_column]();
			for (int j = 0; j < Remainder_column; j++) {
				Mat_I[j * unknown_dm + (num_solution)*num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			}
			Inv_M_A = (complex<double>*)malloc(unknown_dm * Remainder_column * sizeof(complex<double>));
			nrhs = Remainder_column;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
			for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
				for (int j = 0; j < Remainder_column; j++) {
					int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + num_solution * num_col_solution;
					int InvAindex = k + j * unknown_dm;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
			delete[] Mat_I;
			Mat_I = nullptr;
			free(Inv_M_A);
		}

		/*    Columnar calculation      */




//		cnt_n1 = Edge_Boundary[n1];
//		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
//		nrhs = cnt_n1;
//		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
//		for (int i = 0; i < cnt_n1; i++) {
//			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
//		}
//		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
//		mkl_zcsrcoo(job, &unknown_dm, p_pardiso , jp , ip , &nnz_dm, acoo , rowind , colind , &info);
//
//		//for (int i = 0; i < num_domain; i++) {
//		//	MPI_Barrier(MPI_COMM_WORLD);
//		//	if (i == myid) {
//		//		cout << "myid is " << myid << endl;
//		//		cout << "start pardiso "  << endl;
//		//		clock_t beginTime = clock();
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		//		if (error != 0)cout << "Pardiso error " << error << endl;
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		//		clock_t EndTime = clock();
//		//	}
//		//}
//		//cout << "start pardiso" << endl;
//		clock_t beginTime = clock();
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		if (error != 0)cout << "Pardiso error " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		clock_t EndTime = clock();
//
//
//		delete[] Mat_I;
//		Mat_I = nullptr;
//
//		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
//#		pragma omp parallel for 
//		for (int i = 0; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int InvAindex = i + j * unknown_dm;
//				int ReIndex = cnt_n1 * i + j;
//				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
//			}
//		}
//		free(Inv_M_A);
//
//		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
//#		pragma omp parallel for  
//		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
//				int InvAindex = cnt_n1 * i + j;
//				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
//			}
//		}
//		free(R_product_InvA);






		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		//complex<double>* yii = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
		//for (int i = 0; i < unknown_dm; i++) {
		//	yii[i].real(fbr[i].real);
		//	yii[i].imag(fbr[i].imag);
		//}
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure2(int myid) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		nnz_dm = num_nzero_Pmatrix[n1];
		nn2 = 0;
		final_row = r_dm[myid];

		long long unknown_dm_temp = unknown_dm;
		/*        Columnar calculation      */

		cout << "test2" << endl;
		cnt_n1 = Edge_Boundary[n1];
		long long cnt_n1_temp = cnt_n1;
		long long num_col_solution = 1600;
		long long num_solution;
		num_solution = cnt_n1_temp / num_col_solution;
		Mat_I = new complex<double>[unknown_dm_temp * num_col_solution]();

		Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * num_col_solution * sizeof(complex<double>));
		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		nrhs = num_col_solution;
		cout << "test1" << endl;
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		for (long long i = 0; i < num_solution; i++) {
			if (i > 0) {
				for (long long j = 0; j < num_col_solution; j++) {
					Mat_I[j * unknown_dm_temp + (i - 1) * num_col_solution + j + unknown_dm_temp - cnt_n1_temp].real(0.0);//Identity matrix
				}
			}
			for (long long j = 0; j < num_col_solution; j++) {
				Mat_I[j * unknown_dm_temp + i * num_col_solution + j + unknown_dm_temp - cnt_n1_temp].real(1.0);//Identity matrix
			}
			//cout << "test9" << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (long long k = unknown_dm_temp - cnt_n1_temp; k < unknown_dm_temp; k++) {
				for (long long j = 0; j < num_col_solution; j++) {
					long long Zmindex = (k - unknown_dm_temp + cnt_n1_temp) * cnt_n1_temp + j + i * num_col_solution;
					long long InvAindex = k + j * unknown_dm_temp;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
		}
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		cout << "test3" << endl;
		delete[] Mat_I;
		Mat_I = nullptr;
		free(Inv_M_A);
		long long Remainder_column = cnt_n1_temp - num_col_solution * num_solution;
		if (Remainder_column != 0) {
			Mat_I = new complex<double>[unknown_dm_temp * Remainder_column]();
			for (int j = 0; j < Remainder_column; j++) {
				Mat_I[j * unknown_dm_temp + (num_solution)*num_col_solution + j + unknown_dm_temp - cnt_n1_temp].real(1.0);//Identity matrix
			}
			Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * Remainder_column * sizeof(complex<double>));
			nrhs = Remainder_column;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
			for (long long k = unknown_dm - cnt_n1_temp; k < unknown_dm; k++) {
				for (long long j = 0; j < Remainder_column; j++) {
					long long Zmindex = (k - unknown_dm_temp + cnt_n1_temp) * cnt_n1_temp + j + num_solution * num_col_solution;
					long long InvAindex = k + j * unknown_dm_temp;
					mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
				}
			}
			delete[] Mat_I;
			Mat_I = nullptr;
			free(Inv_M_A);
		}

		/*    Columnar calculation      */
		cout << "test4" << endl;



//		cnt_n1 = Edge_Boundary[n1];
//		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
//		nrhs = cnt_n1;
//		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
//		for (int i = 0; i < cnt_n1; i++) {
//			Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
//		}
//		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
//		mkl_zcsrcoo(job, &unknown_dm, p_pardiso , jp , ip , &nnz_dm, acoo , rowind , colind , &info);
//
//		//for (int i = 0; i < num_domain; i++) {
//		//	MPI_Barrier(MPI_COMM_WORLD);
//		//	if (i == myid) {
//		//		cout << "myid is " << myid << endl;
//		//		cout << "start pardiso "  << endl;
//		//		clock_t beginTime = clock();
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		//		if (error != 0)cout << "Pardiso error " << error << endl;
//		//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		//		clock_t EndTime = clock();
//		//	}
//		//}
//		//cout << "start pardiso" << endl;
//		clock_t beginTime = clock();
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		if (error != 0)cout << "Pardiso error " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		clock_t EndTime = clock();
//
//
//		delete[] Mat_I;
//		Mat_I = nullptr;
//
//		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
//#		pragma omp parallel for 
//		for (int i = 0; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int InvAindex = i + j * unknown_dm;
//				int ReIndex = cnt_n1 * i + j;
//				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
//			}
//		}
//		free(Inv_M_A);
//
//		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
//#		pragma omp parallel for  
//		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
//				int InvAindex = cnt_n1 * i + j;
//				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
//			}
//		}
//		free(R_product_InvA);






		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		//complex<double>* yii = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
		//for (int i = 0; i < unknown_dm; i++) {
		//	yii[i].real(fbr[i].real);
		//	yii[i].imag(fbr[i].imag);
		//}
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure_method2(int myid, MKL_Complex16* acoo_temp_1225, int* num_nzero_Pmatrix_E_inter, int* rowind_temp_1225, int* colind_temp_1225) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		nnz_dm = num_nzero_Pmatrix_E_inter[n1];
		cout << "myid is " << myid << "   nnz_dm is " << nnz_dm << endl;
		nn2 = 0;
		final_row = r_dm[myid];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;
		p_pardiso = new MKL_Complex16[nnz_dm];
		ip = new int[unknown_dm + 1];
		jp = new int[nnz_dm];

		int unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];

		int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];

		cnt_n1 = Edge_Boundary[n1] - num_unknown_subdomain[myid][1];
		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
		nrhs = cnt_n1;
		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));

		int nnz_temp_1225 = num_nzero_Pmatrix[myid];
		cout << "myid is " << myid << "   nnz_temp_1225 is " << nnz_temp_1225 << endl;
		cout << "myid is " << myid << "   cnt_n1 is " << cnt_n1 << endl;
		cout << "myid is " << myid << "   unknown_dm is " << unknown_dm << endl;
		int nnz_of_AIS = 0;
		//#		pragma omp parallel for  
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] <= unknown_E_inter && colind[i] > unknown_E_inter) {
				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
				//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].real(acoo[i].real);
				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].imag(acoo[i].imag);
				nnz_of_AIS++;
			}
		}

		//for (int i = 0; i < cnt_n1; i++) {
		//	Mat_I[i * unknown_dm + i + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		//}

		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		delete[] Mat_I;
		Mat_I = nullptr;



		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
#		pragma omp parallel for 
		for (int i = 0; i < unknown_dm; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int InvAindex = i + j * unknown_dm;
				int ReIndex = cnt_n1 * i + j;
				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
			}
		}
		free(Inv_M_A);

		Inv_M_A = nullptr;




		complex<double>* ASI = (complex<double>*)malloc(Edge_Boundary[myid] * unknown_dm * sizeof(complex<double>));
		int nnz_of_ASI = 0;

		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].real(acoo[i].real);
				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].imag(acoo[i].imag);
				nnz_of_ASI++;
			}
		}
		complex<double>* ASI_multi_R_product_InvA = new complex<double>[Edge_Boundary[myid] * cnt_n1];


		MKL_Complex16 alpa1, beta1; alpa1.real = 1.0; alpa1.imag = 0; beta1.real = 0.0; beta1.imag = 0;
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Edge_Boundary[myid], cnt_n1, unknown_dm, &alpa1, ASI, unknown_dm, R_product_InvA, cnt_n1, &beta1, ASI_multi_R_product_InvA, cnt_n1);
		//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn

		free(R_product_InvA);




		//complex<double>* ASS = (complex<double>*)malloc(Edge_Boundary[myid] * Edge_Boundary[myid] * sizeof(complex<double>));
		MKL_Complex16* ASS = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];

		for (int i = 0; i < Edge_Boundary[myid] * Edge_Boundary[myid]; i++) {
			ASS[i].real = 0.0;
			ASS[i].imag = 0.0;
		}
		int nnz_of_ASS = 0;
		//#		pragma omp parallel for  
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] > unknown_E_inter) {
				//if (colind[i] > unknown_E_inter + unknown_E_boundary) {
				//	cout << "error" << endl;
				//}
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].real = acoo[i].real;
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].imag = acoo[i].imag;
				nnz_of_ASS++;
				//acoo_temp_1225[count_nnz_E_inter] = acoo[i];
				//rowind_temp_1225[count_nnz_E_inter] = rowind[i];
				//colind_temp_1225[count_nnz_E_inter] = colind[i];
				//count_nnz_E_inter++;
			}
		}

		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_multi_R_product_InvA[i * cnt_n1 + j].imag();
				ASS[i * Edge_Boundary[myid] + j].real = ASS[i * Edge_Boundary[myid] + j].real - ASI_multi_R_product_InvA[i * cnt_n1 + j].real();
			}

		}


		int* ASS_row = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		int* ASS_col = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			for (int j = 0; j < Edge_Boundary[myid]; j++) {
				ASS_row[i * Edge_Boundary[myid] + j] = i + 1;
				ASS_col[i * Edge_Boundary[myid] + j] = j + 1;
			}
		}

		//cout << "myid is " << myid << "   test6" << endl;

		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;


		p_pardiso = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];
		ip = new int[Edge_Boundary[myid] + 1];
		jp = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];

		nnz_dm = Edge_Boundary[myid] * Edge_Boundary[myid];
		mkl_zcsrcoo(job, &Edge_Boundary[myid], p_pardiso, jp, ip, &nnz_dm, ASS, ASS_row, ASS_col, &info);

		cout << "myid is " << myid << "   test7 " << endl;
		Mat_I = new complex<double>[Edge_Boundary[myid] * Edge_Boundary[myid]]();

		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			Mat_I[i * Edge_Boundary[myid] + i].real(1.0);//Identity matrix
		}

		Inv_M_A = (complex<double>*)malloc(Edge_Boundary[myid] * Edge_Boundary[myid] * sizeof(complex<double>));
		nrhs = Edge_Boundary[myid];
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &Edge_Boundary[myid], p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &Edge_Boundary[myid], p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		//cout << "myid is " << myid << "   test2 " << endl;
		delete[] Mat_I;
		Mat_I = nullptr;

		//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn


		cnt_n1 = Edge_Boundary[myid];
		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		//#		pragma omp parallel for  
		for (int i = 0; i < cnt_n1; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int Zmindex = i * cnt_n1 + j;
				int InvAindex = i + cnt_n1 * j;
				mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
			}
		}
		free(Inv_M_A); Inv_M_A = nullptr;

		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;

		unknown_dm = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];

		int nnz_tot = num_nzero_Pmatrix[myid];
		int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
		p_pardiso = new MKL_Complex16[nnz_tot];
		ip = new int[number_unknown_subdomain + 1];
		jp = new int[nnz_tot];


		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_tot, acoo, rowind, colind, &info);


		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure_method21(int myid, MKL_Complex16* acoo_temp_1225, int* num_nzero_Pmatrix_E_inter, int* rowind_temp_1225, int* colind_temp_1225) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	int startdm = myid;
	int enddm = myid + 1;
	//cout << "myid is " << myid << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	//cout << "Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {
		cnt1 = 0;
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		nnz_dm = num_nzero_Pmatrix_E_inter[n1];
		cout << "myid is " << myid << "   nnz_dm is " << nnz_dm << endl;
		nn2 = 0;
		final_row = r_dm[myid];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;
		p_pardiso = new MKL_Complex16[nnz_dm];
		ip = new int[unknown_dm + 1];
		jp = new int[nnz_dm];

		int unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];

		int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];

		cnt_n1 = Edge_Boundary[n1] - num_unknown_subdomain[myid][1];







		/*        Columnar calculation      */


		//int nnz_temp_1225 = num_nzero_Pmatrix[myid];
		//int nnz_of_ASI = 0;
		//for (int i = 0; i < nnz_temp_1225; i++) {
		//	if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
		//		if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
		//			cout << "error" << endl;
		//		}
		//		nnz_of_ASI++;
		//	}
		//}
		//complex<double>* ASI = (complex<double>*)malloc(nnz_of_ASI * sizeof(complex<double>));
		//int* ASI_row = new int[nnz_of_ASI];
		//int* ASI_col = new int[nnz_of_ASI];
		//nnz_of_ASI = 0;
		//for (int i = 0; i < nnz_temp_1225; i++) {
		//	if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
		//		if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
		//			cout << "error" << endl;
		//		}
		//		ASI[nnz_of_ASI].real(acoo[i].real);
		//		ASI[nnz_of_ASI].imag(acoo[i].imag);
		//		ASI_row[nnz_of_ASI] = rowind[i] - unknown_E_inter;
		//		ASI_col[nnz_of_ASI] = colind[i];
		//		nnz_of_ASI++;
		//	}
		//}
		//MKL_Complex16* ASI_mult_Y = new MKL_Complex16[unknown_E_boundary * unknown_E_boundary];
		//cout << "test2" << endl;
		////cnt_n1 = Edge_Boundary[n1];
		//int num_col_solution = 100;
		//int num_solution;
		//num_solution = unknown_E_boundary / num_col_solution;
		//Mat_I = new complex<double>[unknown_dm * num_col_solution]();
		//Inv_M_A = (complex<double>*)malloc(unknown_dm * num_col_solution * sizeof(complex<double>));
		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		//nrhs = num_col_solution;
		//cout << "test1" << endl;
		//int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		//mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		//int pre_num_solution = num_solution - 1;
		//for (int i = 0; i < pre_num_solution; i++) {
		//	//if (i > 0) {
		//	//	for (int j = 0; j < num_col_solution; j++) {
		//	//		Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
		//	//	}
		//	//}
		//	for (int j = 0; j < unknown_dm * num_col_solution; j++) {
		//		Mat_I[j] = 0.0;
		//	}

		//	for (int j = 0; j < nnz_temp_1225; j++) {
		//		if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + i * num_col_solution && colind[j] <= unknown_E_inter + (i + 1) * num_col_solution) {
		//			if (colind[j] > unknown_E_inter + unknown_E_boundary) {
		//				cout << "error" << endl;
		//			}
		//			//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
		//			//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
		//			Mat_I[(rowind[j] - 1) + (colind[j] - 1 - unknown_E_inter - i * num_col_solution) * unknown_E_inter].real(acoo[j].real);
		//			Mat_I[(rowind[j] - 1) + (colind[j] - 1 - unknown_E_inter - i * num_col_solution) * unknown_E_inter].imag(acoo[j].imag);
		//		}
		//	}

		//	//for (int j = 0; j < num_col_solution; j++) {
		//	//	Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
		//	//}

		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		//	for (int j = 0; j < unknown_dm; j++) {
		//		for (int k = 0; k < num_col_solution; k++) {
		//			int InvAindex = k * unknown_dm + j;
		//			int Aindex = k + j * num_col_solution;
		//			Mat_I[Aindex] = Inv_M_A[InvAindex];

		//		}
		//	}

		//	for (int j = 0; j < nnz_of_ASI; j++) {
		//		int row_temp = ASI_row[j];
		//		int col_temp = ASI_col[j];
		//		for (int k = 0; k < num_col_solution; k++) {
		//			int index_temp = row_temp * unknown_E_boundary + k + i * num_col_solution;
		//			ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
		//			ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();

		//		}


		//	}




		//	//for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
		//	//	for (int j = 0; j < num_col_solution; j++) {
		//	//		int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
		//	//		int InvAindex = k + j * unknown_dm;
		//	//		mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
		//	//	}
		//	//}
		//}

		//delete[] Mat_I;
		//Mat_I = nullptr;
		//free(Inv_M_A);
		//int Remainder_column = cnt_n1 - num_col_solution * pre_num_solution;
		//if (Remainder_column != 0) {
		//	Mat_I = new complex<double>[unknown_dm * Remainder_column]();

		//	for (int j = 0; j < unknown_dm * Remainder_column; j++) {
		//		Mat_I[j] = 0.0;
		//	}

		//	for (int j = 0; j < nnz_temp_1225; j++) {
		//		if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + pre_num_solution * num_col_solution) {
		//			//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
		//			//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
		//			Mat_I[(rowind[j] - 1) + (colind[j] - 1 - unknown_E_inter - pre_num_solution * num_col_solution) * unknown_E_inter].real(acoo[j].real);
		//			Mat_I[(rowind[j] - 1) + (colind[j] - 1 - unknown_E_inter - pre_num_solution * num_col_solution) * unknown_E_inter].imag(acoo[j].imag);
		//		}
		//	}


		//	Inv_M_A = (complex<double>*)malloc(unknown_dm * Remainder_column * sizeof(complex<double>));
		//	nrhs = Remainder_column;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		//	if (error != 0)cout << "Pardiso error " << error << endl;
		//	pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？


		//	for (int j = 0; j < nnz_of_ASI; j++) {
		//		int row_temp = ASI_row[j];
		//		int col_temp = ASI_col[j];
		//		for (int k = 0; k < Remainder_column; k++) {
		//			int index_temp = row_temp * unknown_E_boundary + k + pre_num_solution * num_col_solution;
		//			ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
		//			ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();

		//		}


		//	}
		//	delete[] Mat_I;
		//	Mat_I = nullptr;
		//	free(Inv_M_A);
		//}


		int nnz_temp_1225 = num_nzero_Pmatrix[myid];
		int nnz_of_ASI = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				nnz_of_ASI++;
			}
		}
		complex<double>* ASI = (complex<double>*)malloc(nnz_of_ASI * sizeof(complex<double>));
		int* ASI_row = new int[nnz_of_ASI];
		int* ASI_col = new int[nnz_of_ASI];
		nnz_of_ASI = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				ASI[nnz_of_ASI].real(acoo[i].real);
				ASI[nnz_of_ASI].imag(acoo[i].imag);
				ASI_row[nnz_of_ASI] = rowind[i] - unknown_E_inter;
				ASI_col[nnz_of_ASI] = colind[i];
				nnz_of_ASI++;
			}
		}
		MKL_Complex16* ASI_mult_Y = new MKL_Complex16[unknown_E_boundary * unknown_E_boundary];
		for (int i = 0; i < unknown_E_boundary * unknown_E_boundary; i++) {
			ASI_mult_Y[i].imag = 0.0;
			ASI_mult_Y[i].real = 0.0;
		}

		long long unknown_dm_temp = unknown_dm;

		//cout << "test2" << endl;
		//cnt_n1 = Edge_Boundary[n1];
		long long  num_col_solution = 200;
		int num_solution;
		num_solution = unknown_E_boundary / num_col_solution;
		Mat_I = new complex<double>[unknown_dm_temp * num_col_solution]();
		Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * num_col_solution * sizeof(complex<double>));
		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		nrhs = num_col_solution;
		cout << "test1" << endl;
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		//mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
		int pre_num_solution = num_solution - 1;
		for (int i = 0; i < pre_num_solution; i++) {
			//if (i > 0) {
			//	for (int j = 0; j < num_col_solution; j++) {
			//		Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
			//	}
			//}
			for (long long j = 0; j < unknown_dm_temp * num_col_solution; j++) {
				Mat_I[j] = 0.0;
			}

			for (int j = 0; j < nnz_temp_1225; j++) {
				if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + i * num_col_solution && colind[j] <= unknown_E_inter + (i + 1) * num_col_solution) {
					if (colind[j] > unknown_E_inter + unknown_E_boundary) {
						cout << "error" << endl;
					}
					long long index_mat_I_row = (rowind[j] - 1);
					long long index_mat_I_col = colind[j] - 1 - unknown_E_inter - i * num_col_solution;


					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].real(acoo[j].real);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].imag(acoo[j].imag);
				}
			}

			//for (int j = 0; j < num_col_solution; j++) {
			//	Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			//}

			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (long long j = 0; j < unknown_dm_temp; j++) {
				for (long long k = 0; k < num_col_solution; k++) {
					long long InvAindex = k * unknown_dm_temp + j;
					long long Aindex = k + j * num_col_solution;
					Mat_I[Aindex] = Inv_M_A[InvAindex];

				}
			}

			for (long long j = 0; j < nnz_of_ASI; j++) {
				long long row_temp = ASI_row[j] - 1;//re
				long long col_temp = ASI_col[j] - 1;//re
				for (long long k = 0; k < num_col_solution; k++) {
					long long index_temp = row_temp * unknown_E_boundary + k + i * num_col_solution;
					ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
					ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();
				}
			}


			//for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
			//	for (int j = 0; j < num_col_solution; j++) {
			//		int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
			//		int InvAindex = k + j * unknown_dm;
			//		mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
			//	}
			//}
		}

		//cout << "test3" << endl;


		delete[] Mat_I;
		Mat_I = nullptr;
		free(Inv_M_A);
		long long Remainder_column = cnt_n1 - num_col_solution * pre_num_solution;
		if (Remainder_column != 0) {
			Mat_I = new complex<double>[unknown_dm_temp * Remainder_column]();

			for (long long j = 0; j < unknown_dm_temp * Remainder_column; j++) {
				Mat_I[j] = 0.0;
			}

			for (int j = 0; j < nnz_temp_1225; j++) {
				if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + pre_num_solution * num_col_solution) {
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
					long long index_mat_I_row = (rowind[j] - 1);
					long long index_mat_I_col = colind[j] - 1 - unknown_E_inter - pre_num_solution * num_col_solution;

					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].real(acoo[j].real);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].imag(acoo[j].imag);
				}
			}


			Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * Remainder_column * sizeof(complex<double>));
			nrhs = Remainder_column;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (long long j = 0; j < unknown_dm_temp; j++) {
				for (long long k = 0; k < Remainder_column; k++) {
					long long InvAindex = k * unknown_dm_temp + j;
					long long Aindex = k + j * Remainder_column;
					Mat_I[Aindex] = Inv_M_A[InvAindex];

				}
			}


			for (int j = 0; j < nnz_of_ASI; j++) {
				long long row_temp = ASI_row[j] - 1;
				long long col_temp = ASI_col[j] - 1;
				for (long long k = 0; k < Remainder_column; k++) {
					int index_temp = row_temp * unknown_E_boundary + k + pre_num_solution * num_col_solution;
					ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * Remainder_column + k].real() - ASI[j].imag() * Mat_I[col_temp * Remainder_column + k].imag();
					ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * Remainder_column + k].imag() + ASI[j].imag() * Mat_I[col_temp * Remainder_column + k].real();

				}
			}
			//for (int j = 0; j < nnz_of_ASI; j++) {
			//	long long row_temp = ASI_row[j] - 1;
			//	long long col_temp = ASI_col[j] - 1;
			//	for (long long k = 0; k < Remainder_column; k++) {
			//		int index_temp = row_temp * unknown_E_boundary + k + pre_num_solution * num_col_solution;
			//		ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
			//		ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();

			//	}
			//}

			delete[] Mat_I;
			Mat_I = nullptr;
			free(Inv_M_A);
		}

		/*    Columnar calculation      */

















//		Mat_I = new complex<double>[unknown_dm * cnt_n1]();
//		nrhs = cnt_n1;
//		Inv_M_A = (complex<double>*)malloc(unknown_dm * cnt_n1 * sizeof(complex<double>));
//
//		int nnz_temp_1225 = num_nzero_Pmatrix[myid];
//		cout << "myid is " << myid << "   nnz_temp_1225 is " << nnz_temp_1225 << endl;
//		cout << "myid is " << myid << "   cnt_n1 is " << cnt_n1 << endl;
//		cout << "myid is " << myid << "   unknown_dm is " << unknown_dm << endl;
//		int nnz_of_AIS = 0;
//		//#		pragma omp parallel for  
//		for (int i = 0; i < nnz_temp_1225; i++) {
//			if (rowind[i] <= unknown_E_inter && colind[i] > unknown_E_inter) {
//				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].real(acoo[i].real);
//				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].imag(acoo[i].imag);
//				nnz_of_AIS++;
//			}
//		}
//
//		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
//		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		if (error != 0)cout << "Pardiso error " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//
//		delete[] Mat_I;
//		Mat_I = nullptr;
//
//		//cout << "myid is " << myid << "   Edge_Boundary[myid] is " << Edge_Boundary[myid] << endl;
//
//		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1 * unknown_dm * sizeof(complex<double>));
//#		pragma omp parallel for 
//		for (int i = 0; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int InvAindex = i + j * unknown_dm;
//				int ReIndex = cnt_n1 * i + j;
//				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
//			}
//		}
//		free(Inv_M_A);
//		Inv_M_A = nullptr;
//		complex<double>* ASI = (complex<double>*)malloc(Edge_Boundary[myid] * unknown_dm * sizeof(complex<double>));
//		int nnz_of_ASI = 0;
//		for (int i = 0; i < nnz_temp_1225; i++) {
//			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
//				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
//					cout << "error" << endl;
//				}
//				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].real(acoo[i].real);
//				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].imag(acoo[i].imag);
//				nnz_of_ASI++;
//				//acoo_temp_1225[count_nnz_E_inter] = acoo[i];
//				//rowind_temp_1225[count_nnz_E_inter] = rowind[i];
//				//colind_temp_1225[count_nnz_E_inter] = colind[i];
//				//count_nnz_E_inter++;
//			}
//		}
//		complex<double>* ASI_multi_R_product_InvA = new complex<double>[Edge_Boundary[myid] * cnt_n1];
//		MKL_Complex16 alpa1, beta1; alpa1.real = 1.0; alpa1.imag = 0; beta1.real = 0.0; beta1.imag = 0;
//		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Edge_Boundary[myid], cnt_n1, unknown_dm, &alpa1, ASI, unknown_dm, R_product_InvA, cnt_n1, &beta1, ASI_multi_R_product_InvA, cnt_n1);
//		free(R_product_InvA);






		//complex<double>* ASS = (complex<double>*)malloc(Edge_Boundary[myid] * Edge_Boundary[myid] * sizeof(complex<double>));
		MKL_Complex16* ASS = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];

		for (int i = 0; i < Edge_Boundary[myid] * Edge_Boundary[myid]; i++) {
			ASS[i].real = 0.0;
			ASS[i].imag = 0.0;
		}
		int nnz_of_ASS = 0;
		//#		pragma omp parallel for  
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] > unknown_E_inter) {
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].real = acoo[i].real;
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].imag = acoo[i].imag;
				nnz_of_ASS++;
			}
		}

		for (int i = 0; i < cnt_n1; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_mult_Y[i * cnt_n1 + j].imag;
				ASS[i * Edge_Boundary[myid] + j].real = ASS[i * Edge_Boundary[myid] + j].real - ASI_mult_Y[i * cnt_n1 + j].real;
			}

		}


		int* ASS_row = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		int* ASS_col = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			for (int j = 0; j < Edge_Boundary[myid]; j++) {
				ASS_row[i * Edge_Boundary[myid] + j] = i + 1;
				ASS_col[i * Edge_Boundary[myid] + j] = j + 1;
			}
		}

		//cout << "myid is " << myid << "   test6" << endl;

		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;


		p_pardiso = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];
		ip = new int[Edge_Boundary[myid] + 1];
		jp = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];

		nnz_dm = Edge_Boundary[myid] * Edge_Boundary[myid];
		mkl_zcsrcoo(job, &Edge_Boundary[myid], p_pardiso, jp, ip, &nnz_dm, ASS, ASS_row, ASS_col, &info);

		//cout << "myid is " << myid << "   test7 " << endl;
		Mat_I = new complex<double>[Edge_Boundary[myid] * Edge_Boundary[myid]]();

		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			Mat_I[i * Edge_Boundary[myid] + i].real(1.0);//Identity matrix
		}

		Inv_M_A = (complex<double>*)malloc(Edge_Boundary[myid] * Edge_Boundary[myid] * sizeof(complex<double>));
		nrhs = Edge_Boundary[myid];
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &Edge_Boundary[myid], p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &Edge_Boundary[myid], p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

		//cout << "myid is " << myid << "   test2 " << endl;
		delete[] Mat_I;
		Mat_I = nullptr;

		//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn


		cnt_n1 = Edge_Boundary[myid];
		mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		//#		pragma omp parallel for  
		for (int i = 0; i < cnt_n1; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				int Zmindex = i * cnt_n1 + j;
				int InvAindex = i + cnt_n1 * j;
				mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
			}
		}
		free(Inv_M_A); Inv_M_A = nullptr;

		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;

		unknown_dm = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];

		int nnz_tot = num_nzero_Pmatrix[myid];
		int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
		p_pardiso = new MKL_Complex16[nnz_tot];
		ip = new int[number_unknown_subdomain + 1];
		jp = new int[nnz_tot];


		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_tot, acoo, rowind, colind, &info);


		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		nrhs = 1;


		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
			fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
			fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
		}
		//free(yii);
		free(pro);
		//if (myid == num_domain-1) {
		//	ofstream ofs239("unknownUifbb123.csv", ios::app);
		//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
		//		ofs239 << fbb[i].real << ',' << fbb[i].imag << endl;
		//	}
		//}




		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize];
					pro = new complex<double>[matsize];
#		pragma omp parallel for 
					for (int i = 0; i < matsize; i++) {
						b_pro[i].real(0);
						b_pro[i].imag(0);
						pro[i].real(0);
						pro[i].imag(0);
					}
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
						}
					}
					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;


					for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
						if (pro[i] != (complex<double>)0) {
							int mr = final_row + i / cnt_n2;
							int mc = final_col + i % cnt_n2;
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, pro[i]);
						}
					}
					Tri_tmp.resize(List_step);
					//pthread_mutex_lock(&lock);
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//pthread_mutex_unlock(&lock);
					delete[] b_pro;
					delete[] pro;
				}
				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}

		//final_row += Edge_Boundary[n1];
		delete[] mat_zm;
	}
}

void FETI_like_procedure_method3(int myid, MKL_Complex16* acoo_temp_1225, int* num_nzero_Pmatrix_E_inter, int* rowind_temp_1225, int* colind_temp_1225) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;
	int startdm = myid;
	int enddm = myid + 1;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	//int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int  row1, nn2, nnz_dm, final_row, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int unknown_dm, cnt_n1;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {

		//unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		//nnz_dm = num_nzero_Pmatrix[n1];
		final_row = r_dm[myid];

		int unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];
		cnt_n1 = Edge_Boundary[n1] - num_unknown_subdomain[myid][1];

		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		nnz_dm = num_nzero_Pmatrix_E_inter[n1];
		cout << "myid is " << myid << "   nnz_dm is " << nnz_dm << endl;
		nn2 = 0;
		final_row = r_dm[myid];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;
		p_pardiso = new MKL_Complex16[nnz_dm];
		ip = new int[unknown_dm + 1];
		jp = new int[nnz_dm];



		long long unknown_dm_temp = unknown_dm;
		long long cnt_n1_temp = cnt_n1;
		Mat_I = new complex<double>[unknown_dm_temp * cnt_n1_temp]();
		nrhs = cnt_n1;
		Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * cnt_n1_temp * sizeof(complex<double>));

		int nnz_temp_1225 = num_nzero_Pmatrix[myid];
		cout << "myid is " << myid << "   nnz_temp_1225 is " << nnz_temp_1225 << endl;
		cout << "myid is " << myid << "   cnt_n1 is " << cnt_n1 << endl;
		cout << "myid is " << myid << "   unknown_dm is " << unknown_dm << endl;
		int nnz_of_AIS = 0;
		//#		pragma omp parallel for  
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] <= unknown_E_inter && colind[i] > unknown_E_inter) {
				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
				//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].real(acoo[i].real);
				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].imag(acoo[i].imag);
				nnz_of_AIS++;
			}
		}
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
		if (error != 0)cout << "Pardiso error " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
		delete[] Mat_I;
		Mat_I = nullptr;



		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1_temp * unknown_dm_temp * sizeof(complex<double>));
#		pragma omp parallel for 
		for (long long i = 0; i < unknown_dm_temp; i++) {
			for (long long j = 0; j < cnt_n1_temp; j++) {
				long long InvAindex = i + j * unknown_dm_temp;
				long long ReIndex = cnt_n1_temp * i + j;
				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
			}
		}
		free(Inv_M_A);

		complex<double>* ASI = (complex<double>*)malloc(cnt_n1_temp * unknown_dm_temp * sizeof(complex<double>));
		int nnz_of_ASI = 0;

		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].real(acoo[i].real);
				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].imag(acoo[i].imag);
				nnz_of_ASI++;
			}
		}
		complex<double>* ASI_multi_R_product_InvA = new complex<double>[cnt_n1_temp * cnt_n1_temp];


		MKL_Complex16 alpa1, beta1; alpa1.real = 1.0; alpa1.imag = 0; beta1.real = 0.0; beta1.imag = 0;
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n1, unknown_dm, &alpa1, ASI, unknown_dm, R_product_InvA, cnt_n1, &beta1, ASI_multi_R_product_InvA, cnt_n1);
		free(R_product_InvA);




		MKL_Complex16* ASS = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];

		for (int i = 0; i < Edge_Boundary[myid] * Edge_Boundary[myid]; i++) {
			ASS[i].real = 0.0;
			ASS[i].imag = 0.0;
		}
		int nnz_of_ASS = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] > unknown_E_inter) {
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].real = acoo[i].real;
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].imag = acoo[i].imag;
				if (rowind[i] > unknown_E_inter + cnt_n1 || colind[i] > unknown_E_inter + cnt_n1) {
					nnz_of_ASS++;
				}


			}
		}
		for (int i = 0; i < cnt_n1; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_multi_R_product_InvA[i * cnt_n1 + j].imag();
				ASS[i * Edge_Boundary[myid] + j].real = ASS[i * Edge_Boundary[myid] + j].real - ASI_multi_R_product_InvA[i * cnt_n1 + j].real();
			}

		}

		//if (myid == 1) {
		//	ofstream ofs1224("ASI_multi_R_product_InvA.txt");
		//	for (int i = 0; i < cnt_n1; i++) {
		//		for (int j = 0; j < cnt_n1; j++) {
		//			//ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_mult_Y[i * cnt_n1 + j].imag;
		//			ofs1224 << ASI_multi_R_product_InvA[i * cnt_n1 + j].imag() << ',' << ASI_multi_R_product_InvA[i * cnt_n1 + j].real() << endl;
		//		}
		//	}
		//}

		//int* ASS_row = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		//int* ASS_col = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];


		complex<double>* ASS_temp = new complex<double>[nnz_of_ASS + cnt_n1 * cnt_n1];
		int* ASS_row_temp = new int[nnz_of_ASS + cnt_n1 * cnt_n1];
		int* ASS_col_temp = new int[nnz_of_ASS + cnt_n1 * cnt_n1];

		cout << "nnz_of_ASS + cnt_n1 * cnt_n1 = " << nnz_of_ASS + cnt_n1 * cnt_n1 << endl;

		int count_nnz_ASS = 0;
		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			for (int j = 0; j < Edge_Boundary[myid]; j++) {
				if (ASS[i * Edge_Boundary[myid] + j].imag != 0 || ASS[i * Edge_Boundary[myid] + j].real != 0) {
					ASS_temp[count_nnz_ASS].imag(ASS[i * Edge_Boundary[myid] + j].imag);
					ASS_temp[count_nnz_ASS].real(ASS[i * Edge_Boundary[myid] + j].real);
					ASS_row_temp[count_nnz_ASS] = i + 1;
					ASS_col_temp[count_nnz_ASS] = j + 1;
					count_nnz_ASS++;
				}
				//ASS_row[i * Edge_Boundary[myid] + j] = i + 1;
				//ASS_col[i * Edge_Boundary[myid] + j] = j + 1;

			}
		}

		cout << "count_nnz_ASS = " << count_nnz_ASS << endl;

		delete[]ASS;





		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
//#		pragma omp parallel for  
//		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
//				int InvAindex = cnt_n1 * i + j;
//				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
//			}
//		}
		//free(R_product_InvA);



		MKL_Complex16* fbr_temp = new MKL_Complex16[unknown_dm];
		for (int i = 0; i < unknown_dm; i++) {
			fbr_temp[i].real = fbr[i].real;
			fbr_temp[i].imag = fbr[i].imag;// store boundy fi	
		}



		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr_temp, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr_temp, pro, &error);



		complex<double>* pro2 = (std::complex<double>*)malloc(cnt_n1 * sizeof(std::complex<double>));
		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, unknown_dm, &alpa1, ASI, unknown_dm, pro, 1, &beta1, pro2, 1);

		//if (myid == 1) {
		//	ofstream ofs1224("pro2r.txt");
		//	for (int i = 0; i < cnt_n1; i++) {
		//		ofs1224 << pro2[i].imag() << ',' << pro2[i].real() << endl;
		//	}
		//}

		free(pro);
		delete[]fbr_temp;


		for (int i = 0; i < cnt_n1; i++) {
			fbb[i].real = fbr[i + unknown_dm].real - pro2[i].real();
			fbb[i].imag = fbr[i + unknown_dm].imag - pro2[i].imag();// store boundy fi	
		}
		free(pro2);
		for (int i = cnt_n1; i < Edge_Boundary[myid]; i++) {
			fbb[i].real = fbr[i + unknown_dm].real;
			fbb[i].imag = fbr[i + unknown_dm].imag;// store boundy fi	
		}

		cnt_n1 = Edge_Boundary[myid];
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;


		int nnz_tot = num_nzero_Pmatrix[myid];
		int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
		p_pardiso = new MKL_Complex16[nnz_tot];
		ip = new int[number_unknown_subdomain + 1];
		jp = new int[nnz_tot];
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_tot, acoo, rowind, colind, &info);

		delete[]acoo; delete[]rowind; delete[]colind;


		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = mm2 - nn2;
					//					b_pro = new complex<double>[matsize];
					//					pro = new complex<double>[matsize];
					//#		pragma omp parallel for 
					//					for (int i = 0; i < matsize; i++) {
					//						b_pro[i].real(0);
					//						b_pro[i].imag(0);
					//						pro[i].real(0);
					//						pro[i].imag(0);
					//					}
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							int mr = final_row + mrow[i] - start1 + 1;
							int mc = final_col + mcol[i] - start2 + 1;
							//b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							//b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
							complex<double> comp_temp;
							comp_temp.imag(m[i].imag);
							comp_temp.real(m[i].real);
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, comp_temp);
						}
					}




					//for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
					//	if (pro[i] != (complex<double>)0) {
					//		int mr = final_row + i / cnt_n2;
					//		int mc = final_col + i % cnt_n2;
					//		
					//	}
					//}
					Tri_tmp.resize(List_step);

					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//delete[] b_pro;
					//delete[] pro;
				}
				else {

					vector<Eigen::Triplet<complex<double>>> Tri_tmp(count_nnz_ASS);
					for (int i = 0; i < count_nnz_ASS; i++) {
						Tri_tmp[i] = Triplet<complex<double>>(ASS_row_temp[i] + final_row, ASS_col_temp[i] + final_col, ASS_temp[i]);
					}
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());

				}



				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}
		delete[] ASS_row_temp; delete[] ASS_col_temp; delete[] ASS_temp;
		//final_row += Edge_Boundary[n1];
		//delete[] mat_zm;
	}
}

void FETI_like_procedure_method31(int myid, MKL_Complex16* acoo_temp_1225, int* num_nzero_Pmatrix_E_inter, int* rowind_temp_1225, int* colind_temp_1225) {
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;
	int startdm = myid;
	int enddm = myid + 1;
	fbb = new MKL_Complex16[Edge_Boundary[myid]];
	//int cnt1, unknown_dm, row1, nn2, nnz_dm, final_row, cnt_n1, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int  row1, nn2, nnz_dm, final_row, end1, start1, final_col, mm2, cnt_n2, end2, start2;
	int unknown_dm, cnt_n1;
	int* PointE_p = nullptr;
	int* col_index_invA, * row_index_invA;
	sparse_matrix_t csrP = nullptr, destP = nullptr; sparse_index_base_t indexP; int size_pcsr;
	complex<double>* Mat_I, * Inv_M_A_temp, * Inv_M_A, * mat_zm;
	complex<double>* b_pro, * pro;
	for (int n1 = startdm; n1 < enddm; n1++) {

		//unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		//nnz_dm = num_nzero_Pmatrix[n1];
		final_row = r_dm[myid];

		long long unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];
		cnt_n1 = Edge_Boundary[n1] - num_unknown_subdomain[myid][1];

		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];
		nnz_dm = num_nzero_Pmatrix_E_inter[n1];
		cout << "myid is " << myid << "   nnz_dm is " << nnz_dm << endl;
		nn2 = 0;
		final_row = r_dm[myid];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;
		p_pardiso = new MKL_Complex16[nnz_dm];
		ip = new int[unknown_dm + 1];
		jp = new int[nnz_dm];
		long long unknown_dm_temp = unknown_dm;
		long long cnt_n1_temp = cnt_n1;

		//Mat_I = new complex<double>[unknown_dm_temp * cnt_n1_temp]();
		//nrhs = cnt_n1;
		//Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * cnt_n1_temp * sizeof(complex<double>));
		//int nnz_temp_1225 = num_nzero_Pmatrix[myid];


		/*        Columnar calculation      */
		int nnz_temp_1225 = num_nzero_Pmatrix[myid];
		int nnz_of_ASI = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				nnz_of_ASI++;
			}
		}
		complex<double>* ASI = (complex<double>*)malloc(nnz_of_ASI * sizeof(complex<double>));
		int* ASI_row = new int[nnz_of_ASI];
		int* ASI_col = new int[nnz_of_ASI];
		nnz_of_ASI = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
				if (rowind[i] > unknown_E_inter + unknown_E_boundary) {
					cout << "error" << endl;
				}
				ASI[nnz_of_ASI].real(acoo[i].real);
				ASI[nnz_of_ASI].imag(acoo[i].imag);
				ASI_row[nnz_of_ASI] = rowind[i] - unknown_E_inter;
				ASI_col[nnz_of_ASI] = colind[i];
				nnz_of_ASI++;
			}
		}
		MKL_Complex16* ASI_mult_Y = new MKL_Complex16[unknown_E_boundary * unknown_E_boundary];
		for (int i = 0; i < unknown_E_boundary * unknown_E_boundary; i++) {
			ASI_mult_Y[i].imag = 0.0;
			ASI_mult_Y[i].real = 0.0;
		}



		//cout << "test2" << endl;
		//cnt_n1 = Edge_Boundary[n1];
		long long  num_col_solution = 2000;
		int num_solution;
		num_solution = unknown_E_boundary / num_col_solution;
		Mat_I = new complex<double>[unknown_dm_temp * num_col_solution]();
		Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * num_col_solution * sizeof(complex<double>));
		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
		nrhs = num_col_solution;
		cout << "test1" << endl;
		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
		//mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
		int pre_num_solution = num_solution - 1;
		for (int i = 0; i < pre_num_solution; i++) {
			//if (i > 0) {
			//	for (int j = 0; j < num_col_solution; j++) {
			//		Mat_I[j * unknown_dm + (i - 1) * num_col_solution + j + unknown_dm - cnt_n1].real(0.0);//Identity matrix
			//	}
			//}
			for (long long j = 0; j < unknown_dm_temp * num_col_solution; j++) {
				Mat_I[j] = 0.0;
			}

			for (int j = 0; j < nnz_temp_1225; j++) {
				if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + i * num_col_solution && colind[j] <= unknown_E_inter + (i + 1) * num_col_solution) {
					if (colind[j] > unknown_E_inter + unknown_E_boundary) {
						cout << "error" << endl;
					}
					long long index_mat_I_row = (rowind[j] - 1);
					long long index_mat_I_col = colind[j] - 1 - unknown_E_inter - i * num_col_solution;


					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].real(acoo[j].real);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].imag(acoo[j].imag);
				}
			}

			//for (int j = 0; j < num_col_solution; j++) {
			//	Mat_I[j * unknown_dm + i * num_col_solution + j + unknown_dm - cnt_n1].real(1.0);//Identity matrix
			//}

			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (long long j = 0; j < unknown_dm_temp; j++) {
				for (long long k = 0; k < num_col_solution; k++) {
					long long InvAindex = k * unknown_dm_temp + j;
					long long Aindex = k + j * num_col_solution;
					Mat_I[Aindex] = Inv_M_A[InvAindex];

				}
			}

			for (long long j = 0; j < nnz_of_ASI; j++) {
				long long row_temp = ASI_row[j] - 1;//re
				long long col_temp = ASI_col[j] - 1;//re
				for (long long k = 0; k < num_col_solution; k++) {
					long long index_temp = row_temp * unknown_E_boundary + k + i * num_col_solution;
					ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
					ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();
				}
			}


			//for (int k = unknown_dm - cnt_n1; k < unknown_dm; k++) {
			//	for (int j = 0; j < num_col_solution; j++) {
			//		int Zmindex = (k - unknown_dm + cnt_n1) * cnt_n1 + j + i * num_col_solution;
			//		int InvAindex = k + j * unknown_dm;
			//		mat_zm[Zmindex] = Inv_M_A[InvAindex];  //gain Ri * A(-1) *RiT
			//	}
			//}
		}

		//cout << "test3" << endl;


		delete[] Mat_I;
		Mat_I = nullptr;
		free(Inv_M_A);
		long long Remainder_column = cnt_n1 - num_col_solution * pre_num_solution;
		if (Remainder_column != 0) {
			Mat_I = new complex<double>[unknown_dm_temp * Remainder_column]();

			for (long long j = 0; j < unknown_dm_temp * Remainder_column; j++) {
				Mat_I[j] = 0.0;
			}

			for (int j = 0; j < nnz_temp_1225; j++) {
				if (rowind[j] <= unknown_E_inter && colind[j] > unknown_E_inter + pre_num_solution * num_col_solution) {
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].real(acoo[i].real);
					//Mat_I[(rowind[i] - 1) * cnt_n1 + colind[i] - 1 - unknown_E_inter].imag(acoo[i].imag);
					long long index_mat_I_row = (rowind[j] - 1);
					long long index_mat_I_col = colind[j] - 1 - unknown_E_inter - pre_num_solution * num_col_solution;

					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].real(acoo[j].real);
					Mat_I[index_mat_I_row + index_mat_I_col * unknown_E_inter].imag(acoo[j].imag);
				}
			}


			Inv_M_A = (complex<double>*)malloc(unknown_dm_temp * Remainder_column * sizeof(complex<double>));
			nrhs = Remainder_column;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
			if (error != 0)cout << "Pardiso error " << error << endl;
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？

			for (long long j = 0; j < unknown_dm_temp; j++) {
				for (long long k = 0; k < Remainder_column; k++) {
					long long InvAindex = k * unknown_dm_temp + j;
					long long Aindex = k + j * Remainder_column;
					Mat_I[Aindex] = Inv_M_A[InvAindex];

				}
			}


			for (int j = 0; j < nnz_of_ASI; j++) {
				long long row_temp = ASI_row[j] - 1;
				long long col_temp = ASI_col[j] - 1;
				for (long long k = 0; k < Remainder_column; k++) {
					int index_temp = row_temp * unknown_E_boundary + k + pre_num_solution * num_col_solution;
					ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * Remainder_column + k].real() - ASI[j].imag() * Mat_I[col_temp * Remainder_column + k].imag();
					ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * Remainder_column + k].imag() + ASI[j].imag() * Mat_I[col_temp * Remainder_column + k].real();

				}
			}
			//for (int j = 0; j < nnz_of_ASI; j++) {
			//	long long row_temp = ASI_row[j] - 1;
			//	long long col_temp = ASI_col[j] - 1;
			//	for (long long k = 0; k < Remainder_column; k++) {
			//		int index_temp = row_temp * unknown_E_boundary + k + pre_num_solution * num_col_solution;
			//		ASI_mult_Y[index_temp].real += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].real() - ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].imag();
			//		ASI_mult_Y[index_temp].imag += ASI[j].real() * Mat_I[col_temp * num_col_solution + k].imag() + ASI[j].imag() * Mat_I[col_temp * num_col_solution + k].real();

			//	}
			//}

			delete[] Mat_I;
			Mat_I = nullptr;
			free(Inv_M_A);
		}

		/*    Columnar calculation      */




//		int nnz_of_AIS = 0;
//		//#		pragma omp parallel for  
//		for (int i = 0; i < nnz_temp_1225; i++) {
//			if (rowind[i] <= unknown_E_inter && colind[i] > unknown_E_inter) {
//				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
//					cout << "error" << endl;
//				}
//				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].real(acoo[i].real);
//				Mat_I[(rowind[i] - 1) + (colind[i] - 1 - unknown_E_inter) * unknown_E_inter].imag(acoo[i].imag);
//				nnz_of_AIS++;
//			}
//		}
//
//
//		int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
//		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo_temp_1225, rowind_temp_1225, colind_temp_1225, &info);
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error); // x 为 A（-1）
//		if (error != 0)cout << "Pardiso error " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, Mat_I, Inv_M_A, &error);   //释放内存（存储内部格式的矩阵）？
//		delete[] Mat_I;
//		Mat_I = nullptr;
//
//
//
//		complex<double>* R_product_InvA = (complex<double>*)malloc(cnt_n1_temp * unknown_dm_temp * sizeof(complex<double>));
//#		pragma omp parallel for 
//		for (long long i = 0; i < unknown_dm_temp; i++) {
//			for (long long j = 0; j < cnt_n1_temp; j++) {
//				long long InvAindex = i + j * unknown_dm_temp;
//				long long ReIndex = cnt_n1_temp * i + j;
//				R_product_InvA[ReIndex] = Inv_M_A[InvAindex];  //得到R1 * A(-1) 
//			}
//		}
//		free(Inv_M_A);
//
//		complex<double>* ASI = (complex<double>*)malloc(cnt_n1_temp * unknown_dm_temp * sizeof(complex<double>));
//		int nnz_of_ASI = 0;
//
//		for (int i = 0; i < nnz_temp_1225; i++) {
//			if (rowind[i] > unknown_E_inter && colind[i] <= unknown_E_inter) {
//				if (colind[i] > unknown_E_inter + unknown_E_boundary) {
//					cout << "error" << endl;
//				}
//				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].real(acoo[i].real);
//				ASI[(rowind[i] - 1 - unknown_E_inter) * unknown_dm + colind[i] - 1].imag(acoo[i].imag);
//				nnz_of_ASI++;
//			}
//		}
//		complex<double>* ASI_multi_R_product_InvA = new complex<double>[cnt_n1_temp * cnt_n1_temp];
//
//		//sparse_status_t mkl_sparse_z_mm(const sparse_operation_t operation, const
//		//	MKL_Complex16 alpha, const sparse_matrix_t A, const struct matrix_descr descr, const
//		//	sparse_layout_t layout, const MKL_Complex16 * B, const MKL_INT columns, const MKL_INT
//		//	ldb, const MKL_Complex16 beta, MKL_Complex16 * C, const MKL_INT ldc);
//
//
//
//
//		MKL_Complex16 alpa1, beta1; alpa1.real = 1.0; alpa1.imag = 0; beta1.real = 0.0; beta1.imag = 0;
//		cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n1, unknown_dm, &alpa1, ASI, unknown_dm, R_product_InvA, cnt_n1, &beta1, ASI_multi_R_product_InvA, cnt_n1);
//		free(R_product_InvA);




		MKL_Complex16* ASS = new MKL_Complex16[Edge_Boundary[myid] * Edge_Boundary[myid]];

		for (int i = 0; i < Edge_Boundary[myid] * Edge_Boundary[myid]; i++) {
			ASS[i].real = 0.0;
			ASS[i].imag = 0.0;
		}
		int nnz_of_ASS = 0;
		for (int i = 0; i < nnz_temp_1225; i++) {
			if (rowind[i] > unknown_E_inter && colind[i] > unknown_E_inter) {
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].real = acoo[i].real;
				ASS[(rowind[i] - 1 - unknown_E_inter) * Edge_Boundary[myid] + colind[i] - 1 - unknown_E_inter].imag = acoo[i].imag;
				if (rowind[i] > unknown_E_inter + cnt_n1 || colind[i] > unknown_E_inter + cnt_n1) {
					nnz_of_ASS++;
				}


			}
		}


		for (int i = 0; i < cnt_n1; i++) {
			for (int j = 0; j < cnt_n1; j++) {
				ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_mult_Y[i * cnt_n1 + j].imag;
				ASS[i * Edge_Boundary[myid] + j].real = ASS[i * Edge_Boundary[myid] + j].real - ASI_mult_Y[i * cnt_n1 + j].real;
			}

		}





		//int* ASS_row = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];
		//int* ASS_col = new int[Edge_Boundary[myid] * Edge_Boundary[myid]];


		complex<double>* ASS_temp = new complex<double>[nnz_of_ASS + cnt_n1 * cnt_n1];
		int* ASS_row_temp = new int[nnz_of_ASS + cnt_n1 * cnt_n1];
		int* ASS_col_temp = new int[nnz_of_ASS + cnt_n1 * cnt_n1];

		cout << "nnz_of_ASS + cnt_n1 * cnt_n1 = " << nnz_of_ASS + cnt_n1 * cnt_n1 << endl;

		int count_nnz_ASS = 0;
		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			for (int j = 0; j < Edge_Boundary[myid]; j++) {
				if (ASS[i * Edge_Boundary[myid] + j].imag != 0 || ASS[i * Edge_Boundary[myid] + j].real != 0) {
					ASS_temp[count_nnz_ASS].imag(ASS[i * Edge_Boundary[myid] + j].imag);
					ASS_temp[count_nnz_ASS].real(ASS[i * Edge_Boundary[myid] + j].real);
					ASS_row_temp[count_nnz_ASS] = i + 1;
					ASS_col_temp[count_nnz_ASS] = j + 1;
					count_nnz_ASS++;
				}
				//ASS_row[i * Edge_Boundary[myid] + j] = i + 1;
				//ASS_col[i * Edge_Boundary[myid] + j] = j + 1;

			}
		}

		cout << "count_nnz_ASS = " << count_nnz_ASS << endl;

		delete[]ASS;


		if (myid == 1) {
			ofstream ofs1224("ASI_mult_Y.txt");
			for (int i = 0; i < cnt_n1; i++) {
				for (int j = 0; j < cnt_n1; j++) {
					//ASS[i * Edge_Boundary[myid] + j].imag = ASS[i * Edge_Boundary[myid] + j].imag - ASI_mult_Y[i * cnt_n1 + j].imag;
					ofs1224 << ASI_mult_Y[i * cnt_n1 + j].imag << ',' << ASI_mult_Y[i * cnt_n1 + j].real << endl;
				}

			}
			//for (int i = 0; i < cnt_n1; i++) {
			//	ofs1224 << pro2[i].imag() << ',' << pro2[i].real() << endl;
			//}
		}


		//mat_zm = new complex<double>[cnt_n1 * cnt_n1];
//#		pragma omp parallel for  
//		for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
//			for (int j = 0; j < cnt_n1; j++) {
//				int Zmindex = (i - unknown_dm + cnt_n1) * cnt_n1 + j;
//				int InvAindex = cnt_n1 * i + j;
//				mat_zm[Zmindex] = R_product_InvA[InvAindex];  //gain Ri * A(-1) *RiT
//			}
//		}
		//free(R_product_InvA);



		MKL_Complex16* fbr_temp = new MKL_Complex16[unknown_dm];
		for (int i = 0; i < unknown_dm; i++) {
			fbr_temp[i].real = fbr[i].real;
			fbr_temp[i].imag = fbr[i].imag;// store boundary fi	
		}



		pro = (std::complex<double>*)malloc(unknown_dm * sizeof(std::complex<double>));
		nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr_temp, pro, &error);// pro 为 A（-1）*fi
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr_temp, pro, &error);


		//MKL_Complex16 alpa1, beta1; alpa1.real = 1.0; alpa1.imag = 0; beta1.real = 0.0; beta1.imag = 0;
		complex<double>* pro2 = (std::complex<double>*)malloc(cnt_n1 * sizeof(std::complex<double>));

		for (int j = 0; j < cnt_n1; j++) {
			pro2[j].real(0.0);
			pro2[j].imag(0.0);
		}



		//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, unknown_dm, &alpa1, ASI, unknown_dm, pro, 1, &beta1, pro2, 1);//?
		for (int j = 0; j < nnz_of_ASI; j++) {
			long long row_temp = ASI_row[j] - 1;
			long long col_temp = ASI_col[j] - 1;
			pro2[row_temp].real(pro2[row_temp].real() + ASI[j].real() * pro[col_temp].real() - ASI[j].imag() * pro[col_temp].imag());
			pro2[row_temp].imag(pro2[row_temp].imag() + ASI[j].real() * pro[col_temp].imag() + ASI[j].imag() * pro[col_temp].real());
		}

		free(ASI);
		free(pro);
		delete[]fbr_temp;
		//if (myid == 1) {
		//	ofstream ofs1224("pro2.txt");
		//	for (int i = 0; i < cnt_n1; i++) {
		//		ofs1224 << pro2[i].imag() << ',' << pro2[i].real() << endl;
		//	}
		//}



		//cout << "test4" << endl;
		for (int i = 0; i < cnt_n1; i++) {
			fbb[i].real = fbr[i + unknown_dm].real - pro2[i].real();
			fbb[i].imag = fbr[i + unknown_dm].imag - pro2[i].imag();// store boundary fi	
		}
		free(pro2);
		for (int i = cnt_n1; i < Edge_Boundary[myid]; i++) {
			fbb[i].real = fbr[i + unknown_dm].real;
			fbb[i].imag = fbr[i + unknown_dm].imag;// store boundary fi	
		}

		cnt_n1 = Edge_Boundary[myid];
		unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		delete[] p_pardiso; delete[] ip; delete[] jp;
		p_pardiso = nullptr;
		ip = nullptr;
		jp = nullptr;


		int nnz_tot = num_nzero_Pmatrix[myid];
		int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
		p_pardiso = new MKL_Complex16[nnz_tot];
		ip = new int[number_unknown_subdomain + 1];
		jp = new int[nnz_tot];
		mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_tot, acoo, rowind, colind, &info);

		delete[]acoo; delete[]rowind; delete[]colind;

		//cout << "test5" << endl;

		end1 = unknown_dm;
		start1 = end1 - cnt_n1 + 1;
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				mm2 = nn2 + nnz_c[n1][n2];
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];;
					start2 = end2 - cnt_n2 + 1;
					int matsize = mm2 - nn2;
					//					b_pro = new complex<double>[matsize];
					//					pro = new complex<double>[matsize];
					//#		pragma omp parallel for 
					//					for (int i = 0; i < matsize; i++) {
					//						b_pro[i].real(0);
					//						b_pro[i].imag(0);
					//						pro[i].real(0);
					//						pro[i].imag(0);
					//					}
					vector<Eigen::Triplet<complex<double>>> Tri_tmp(matsize);
					size_t List_step = 0;
					//#		pragma omp parallel for 
					for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
							int mr = final_row + mrow[i] - start1 + 1;
							int mc = final_col + mcol[i] - start2 + 1;
							//b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].real(m[i].real);
							//b_pro[cnt_n2 * (mrow[i] - start1) + mcol[i] - start2].imag(m[i].imag);
							complex<double> comp_temp;
							comp_temp.imag(m[i].imag);
							comp_temp.real(m[i].real);
							Tri_tmp[List_step++] = Triplet<complex<double>>(mr, mc, comp_temp);
						}
					}




					//for (int i = 0; i < matsize; i++) {//将Zm*Cij的内容存储下来
					//	if (pro[i] != (complex<double>)0) {
					//		int mr = final_row + i / cnt_n2;
					//		int mc = final_col + i % cnt_n2;
					//		
					//	}
					//}
					Tri_tmp.resize(List_step);

					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());
					//delete[] b_pro;
					//delete[] pro;
				}
				else {

					vector<Eigen::Triplet<complex<double>>> Tri_tmp(count_nnz_ASS);
					for (int i = 0; i < count_nnz_ASS; i++) {
						Tri_tmp[i] = Triplet<complex<double>>(ASS_row_temp[i] + final_row, ASS_col_temp[i] + final_col, ASS_temp[i]);
					}
					TriList.insert(TriList.end(), Tri_tmp.begin(), Tri_tmp.end());

				}



				nn2 = mm2;
			}
			final_col += Edge_Boundary[n2];
		}
		delete[] ASS_row_temp; delete[] ASS_col_temp; delete[] ASS_temp;
		//final_row += Edge_Boundary[n1];
		//delete[] mat_zm;
	}
}

void* Matrix_Generator_JSIM6(int* nnz_of_C, int myid, MKL_Complex16* CC, complex<double>* fbb_temp, complex<double>* fbb_total) {
	for (int n1 = myid; n1 < myid + 1; n1++) {
		int cnt_n1 = Edge_Boundary[n1];
		int mm4 = 0;
		int nn4 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			mm4 = nn4 + Edge_Boundary[n2];
			if (nnz_c[n1][n2] != 0) {
				if (n1 != n2) {
					int cnt_n2 = Edge_Boundary[n2];
					complex<double>* pro1 = new complex<double>[cnt_n1];
					for (int i = 0; i < cnt_n1; i++) {
						pro1[i].real(0.0);
						pro1[i].imag(0.0);
					}

					complex<double>* Un_temp = new complex<double>[cnt_n2]();
					for (int i = 0; i < cnt_n2; i++) {
						Un_temp[i].real(fbb_total[nn4 + i].real());
						Un_temp[i].imag(fbb_total[nn4 + i].imag());
					}

					int  start_of_C = nnz_of_C[n2];
					//#					pragma omp parallel for 
					//					for (int i = 0; i < cnt_n1; i++) {
					//						for (int j = 0; j < cnt_n2; j++) {
					//							pro1[i].real(pro1[i].real() + CC[start_of_C + i * cnt_n2 + j].real * Un_temp[j].real() - CC[start_of_C + i * cnt_n2 + j].imag * Un_temp[j].imag());
					//							pro1[i].imag(pro1[i].imag() + CC[start_of_C + i * cnt_n2 + j].real * Un_temp[j].imag() + CC[start_of_C + i * cnt_n2 + j].imag * Un_temp[j].real());
					//
					//						}
					//					}

					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, &alpa, CC + start_of_C, cnt_n2, fbb_total + nn4, 1, &beta, pro1, 1);
					for (int i = 0; i < cnt_n1; i++) {   // ym - Cmn * Rt * un
						fbb_temp[i].real(fbb_temp[i].real() - pro1[i].real());
						fbb_temp[i].imag(fbb_temp[i].imag() - pro1[i].imag());
					}
					delete[] pro1;
				}
			}
			nn4 = mm4;
		}
	}

	return nullptr;
}

void* Matrix_Generator_CY2(int matsize_of_C, int* nnz_of_C, int myid, MKL_Complex16* CC) {

	int final_row, cnt_n1, final_col, cnt_n2;
	int mat_size = TriList.size();
	int count = 0;
	for (int n1 = myid; n1 < myid + 1; n1++) {
		final_row = 0;
		cnt_n1 = Edge_Boundary[n1];
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					int rowmm = 0;
					int colmm = 0;
					int j = 0;
					int k = 0;
					int start_of_C = nnz_of_C[n2];
					for (int i = 0; i < mat_size; i++) {
						rowmm = TriList[i].row() + 1;
						colmm = TriList[i].col() + 1;
						if (colmm - final_col < cnt_n2 + 1 && colmm - final_col>0) {
							int location_of_C = start_of_C + (rowmm - 1) * cnt_n2 + colmm - final_col - 1;
							CC[location_of_C].imag = TriList[i].value().imag();
							CC[location_of_C].real = TriList[i].value().real();
							count++;
						}

					}
				}
			}
			final_col += Edge_Boundary[n2];
		}
		//final_row += Edge_Boundary[n1];
	}
	//cout << "mat_size = " << mat_size << endl;
	//cout << "count = " << count << endl;
	return nullptr;
}

void* Matrix_Generator_JSIM_long(long long* nnz_of_C, int myid, MKL_Complex16* CC, complex<double>* fbb_temp, complex<double>* fbb_total) {
	for (int n1 = myid; n1 < myid + 1; n1++) {
		int cnt_n1 = Edge_Boundary[n1];
		int mm4 = 0;
		int nn4 = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			mm4 = nn4 + Edge_Boundary[n2];
			if (nnz_c[n1][n2] != 0) {
				if (n1 != n2) {
					int cnt_n2 = Edge_Boundary[n2];
					complex<double>* pro1 = new complex<double>[cnt_n1];
					for (int i = 0; i < cnt_n1; i++) {
						pro1[i].real(0.0);
						pro1[i].imag(0.0);
					}

					complex<double>* Un_temp = new complex<double>[cnt_n2]();
					for (int i = 0; i < cnt_n2; i++) {
						Un_temp[i].real(fbb_total[nn4 + i].real());
						Un_temp[i].imag(fbb_total[nn4 + i].imag());
					}

					long long  start_of_C = nnz_of_C[n2];
#					pragma omp parallel for 
					for (int i = 0; i < cnt_n1; i++) {
						for (int j = 0; j < cnt_n2; j++) {
							long long index_of_C = start_of_C + (i * cnt_n2 + j);
							pro1[i].real(pro1[i].real() + CC[index_of_C].real * Un_temp[j].real() - CC[index_of_C].imag * Un_temp[j].imag());
							pro1[i].imag(pro1[i].imag() + CC[index_of_C].real * Un_temp[j].imag() + CC[index_of_C].imag * Un_temp[j].real());

						}
					}

					//MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, &alpa, CC + start_of_C, cnt_n2, fbb_total + nn4, 1, &beta, pro1, 1);
					for (int i = 0; i < cnt_n1; i++) {   // ym - Cmn * Rt * un
						fbb_temp[i].real(fbb_temp[i].real() - pro1[i].real());
						fbb_temp[i].imag(fbb_temp[i].imag() - pro1[i].imag());
					}
					delete[] pro1;
				}
			}
			nn4 = mm4;
		}
	}

	return nullptr;
}

void* Matrix_Generator_CY_long(long long matsize_of_C, long long* nnz_of_C, int myid, MKL_Complex16* CC) {

	int final_row, cnt_n1, final_col, cnt_n2;
	long long mat_size = TriList.size();
	cout << "long long mat_size = " << mat_size << endl;

	//int count = 0;
	for (int n1 = myid; n1 < myid + 1; n1++) {
		final_row = 0;
		cnt_n1 = Edge_Boundary[n1];
		final_col = 0;
		for (int n2 = 0; n2 < num_domain; n2++) {
			if (nnz_c[n1][n2] != 0) {
				if (n1 != n2) {
					cnt_n2 = Edge_Boundary[n2];
					long long cnt_n2_long = Edge_Boundary[n2];
					long long rowmm = 0;
					long long colmm = 0;
					int j = 0;
					int k = 0;
					long long start_of_C = nnz_of_C[n2];
					for (long long i = 0; i < mat_size; i++) {
						rowmm = TriList[i].row() + 1;
						colmm = TriList[i].col() + 1;
						if (colmm - final_col < cnt_n2 + 1 && colmm - final_col>0) {

							long long location_of_C = start_of_C + (rowmm - 1) * cnt_n2_long + colmm - final_col - 1;
							CC[location_of_C].imag = TriList[i].value().imag();
							CC[location_of_C].real = TriList[i].value().real();
							//count++;
						}

					}
				}
			}
			final_col += Edge_Boundary[n2];
		}
		//final_row += Edge_Boundary[n1];
	}
	//cout << "mat_size = " << mat_size << endl;
	//cout << "count = " << count << endl;
	return nullptr;
}

void* Solver_mpi(int myid) {


	
	int startdm = myid;
	int enddm = myid+1;

	int maxfct, mnum, mtype, phase1, phase2, n, nrhs, msglvl, error;
	int* perm = nullptr;
	int iparm[64];
	int pt[64];
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;

	int nn4 = 0, mm4, unknown_dm2;
	//MKL_Complex16* Un_temp;

	int  nn3 = 0;
	int nn2 = 0;
	int nzero1 = 0; int  row1 = 0; nrhs = 1;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}

	iparm[0] = 1;
	iparm[1] = 3;
	//iparm[2] = 16;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 2;


	for (int n1 = startdm; n1 < enddm; ++n1) {
		int mm3 = 0;
		mm4 = 0; nn4 = 0;
		int unknown_dm = num_unknown_subdomain[n1][0] + num_unknown_subdomain[n1][1];
		int nnz_dm = num_nzero_Pmatrix[n1];
		//int tot_dm = unknown_dm * unknown_dm;
		nn3 = 0;
		int cnt_n2 = 0;
		int cnt_n1 = Edge_Boundary[n1];//
		complex<double>* b_pro, * pro, * Un_temp;

		for (int ndm2 = 0; ndm2 != num_domain; ++ndm2) {
			cnt_n2 = Edge_Boundary[ndm2];
			mm4 = nn4 + cnt_n2;
			if (nnz_c[n1][ndm2] != 0) {
				mm3 = nn3 + nnz_c[n1][ndm2];
				unknown_dm2 = num_unknown_subdomain[ndm2][0] + num_unknown_subdomain[ndm2][1];
				if (n1 != ndm2) {

					long long cnt_n1_long = cnt_n1;
					long long cnt_n2_long = cnt_n2;
					long  long matsize = cnt_n1_long * cnt_n2_long;

					//int matsize = cnt_n1 * cnt_n2;
					b_pro = new complex<double>[matsize]();
					pro = new complex<double>[cnt_n1]();
					Un_temp = new complex<double>[cnt_n2]();
					for (int i = 0; i < cnt_n2; i++) {
						Un_temp[i].real(unKnownUi[nn4 + i].real);//R*u
						Un_temp[i].imag(unKnownUi[nn4 + i].imag);
					}

					for (int i = nn3; i < mm3; i++) {
						if (mrow[i] > unknown_dm - cnt_n1 && mcol[i] > unknown_dm2 - cnt_n2) {
							long long long_temp1 = mrow[i] - unknown_dm + cnt_n1 - 1;
							long long long_temp2 = mcol[i] - unknown_dm2 + cnt_n2 - 1;

							b_pro[cnt_n2_long * (long_temp1)+long_temp2].real(m[i].real);
							b_pro[cnt_n2_long * (long_temp1)+long_temp2].imag(m[i].imag);
						}

					}

					MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
					cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, &alpa, b_pro, cnt_n2, Un_temp, 1, &beta, pro, 1);   //Cmn * Rt * un
					
					for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {   // ym - Cmn * Rt * un
						fbr[i].real -= pro[i - unknown_dm + cnt_n1].real();
						fbr[i].imag -= pro[i - unknown_dm + cnt_n1].imag();
					}
					delete[] b_pro;
					delete[] pro;
					delete[] Un_temp;

				}
				nn3 = mm3;
			}
			nn4 = mm4;
		}
		//cout << "inter Pardiso" << endl;
		int nrhs = 1;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso , ip , jp , perm, &nrhs, iparm, &msglvl, fbr , Xsol , &error);
		if (error != 0)	cout << "eror = " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso , ip , jp , perm, &nrhs, iparm, &msglvl, fbr , Xsol , &error);
	}
	//return 0;
	return nullptr;//gai
}

int Fbr_Setting(Element* el, int myid) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	//double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();
	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	fbr = (MKL_Complex16*)malloc(sizeof(MKL_Complex16) * number_unknown_subdomain);
	for (int i = 0; i < number_unknown_subdomain; ++i) {
		fbr[i].real = 0.0;
		fbr[i].imag = 0.0;
	}

	double disx, disy, r;
	double Kwave[3];
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;
	cout << " *********************************" << powerful << endl;
	cout << " powerful = " << powerful << endl;
	cout << " ********************************* " << powerful << endl;
	double yc1;
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);
	for (int ith = 0; ith < num_element_subdomain; ith++) {

		//if (el[ith].Material == 1) {
		//	continue;
		//}
		double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
		//if (el[ith].Material != wave_material) {
		if (zt < wave_zzz) {
			continue;
		}


		Kwave[2] = omega * sqrt(el[ith].epsil * mur);
		for (nn = 0; nn < 4; nn++) {
			op1 = el[ith].face[nn].opp[0]; ofn1 = el[ith].face[nn].opp[1];
			if ((ofn1 == -200)) {
				xrod1 = x_wave[0]; yrod1 = y_wave[0];
				//xrod1 = 6.6e-3; yrod1 = 0;
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
							disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
							r = sqrt(pow(disx, 2) + pow(disy, 2));
							//value of test basis function
							Li1[0] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].b[node_ii2 - 1];
							Li1[1] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].c[node_ii2 - 1];
							Li1[2] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].d[node_ii2 - 1];

							Li2[0] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].b[node_ii1 - 1];
							Li2[1] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].c[node_ii1 - 1];
							Li2[2] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].d[node_ii1 - 1];

							//Incident E has Z component;
							EincX.real((cos(Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
							EincX.imag((-sin(Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
							EincY.real((cos(Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
							EincY.imag((-sin(Kwave[2] * rs[2])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
							EincZ.real(0.0);
							EincZ.imag(0.0);

							Einc[0].real(EincX.real() / el[ith].eta); Einc[0].imag(EincX.imag() / el[ith].eta);
							Einc[1].real(EincY.real() / el[ith].eta); Einc[1].imag(EincY.imag() / el[ith].eta);
							Einc[2].real(EincZ.real() / el[ith].eta); Einc[2].imag(EincZ.imag() / el[ith].eta);
							ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
							ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
							ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];
							ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
							ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
							ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];
							for (int ii0 = 0; ii0 < 3; ii0++) {
								//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
								//fbr[edgeii_E - 1].imag += 2 * omega* signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].real() + ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
								//fbr[edgeii_E - 1].imag -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].imag());
								//fbr[edgeii_E - 1].real -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].real());
								//fbr[edgeii_E - 1].real+= - 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].imag()) * el[ith].face[nn].Area * weight0[nGauss];
								//fbr[edgeii_E - 1].imag+= + 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
								fbr[edgeii_E - 1].real = (fbr[edgeii_E - 1].real - powerful * 2 * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].imag() * el[ith].face[nn].Area * weight0[nGauss]);
								fbr[edgeii_E - 1].imag = (fbr[edgeii_E - 1].imag + powerful * 2 * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].real() * el[ith].face[nn].Area * weight0[nGauss]);
							}
						}
					}
				}
			}

		}
	}
	//cout << "omega * mur * w_wg / pi  = " << omega * mur * w_wg / pi << endl;
	return 0;
}

int Fbr_Setting2(Element* el, int myid) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	//double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();
	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	fbr = (MKL_Complex16*)malloc(sizeof(MKL_Complex16) * number_unknown_subdomain);
	for (int i = 0; i < number_unknown_subdomain; ++i) {
		fbr[i].real = 0.0;
		fbr[i].imag = 0.0;
	}

	double disx, disy, r;
	double Kwave[3];
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;


	double yc1;
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);


	for (int i = 0; i < num_wave; i++) {
		for (int ith = 0; ith < num_element_subdomain; ith++) {

			double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
			//if (el[ith].Material != wave_material) {
			if (zt < wave_zzz) {
				continue;
			}
			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
			for (nn = 0; nn < 4; nn++) {
				op1 = el[ith].face[nn].opp[0]; ofn1 = el[ith].face[nn].opp[1];
				if (ofn1 == (-200 - i)) {
					xrod1 = x_wave[i]; yrod1 = y_wave[i];
					//xrod1 = 6.6e-3; yrod1 = 0;
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
								disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
								r = sqrt(pow(disx, 2) + pow(disy, 2));
								//value of test basis function
								Li1[0] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].b[node_ii2 - 1];
								Li1[1] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].c[node_ii2 - 1];
								Li1[2] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].d[node_ii2 - 1];

								Li2[0] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].b[node_ii1 - 1];
								Li2[1] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].c[node_ii1 - 1];
								Li2[2] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].d[node_ii1 - 1];

								//Incident E has Z component;
								EincX.real((cos(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
								EincX.imag((-sin(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
								EincY.real((cos(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
								EincY.imag((-sin(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
								EincZ.real(0.0);
								EincZ.imag(0.0);

								Einc[0].real(EincX.real() / el[ith].eta); Einc[0].imag(EincX.imag() / el[ith].eta);
								Einc[1].real(EincY.real() / el[ith].eta); Einc[1].imag(EincY.imag() / el[ith].eta);
								Einc[2].real(EincZ.real() / el[ith].eta); Einc[2].imag(EincZ.imag() / el[ith].eta);
								ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
								ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
								ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];
								ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
								ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
								ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
									//fbr[edgeii_E - 1].imag += 2 * omega* signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].real() + ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
									//fbr[edgeii_E - 1].imag -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].imag());
									//fbr[edgeii_E - 1].real -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].real());
									//fbr[edgeii_E - 1].real+= - 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].imag()) * el[ith].face[nn].Area * weight0[nGauss];
									//fbr[edgeii_E - 1].imag+= + 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
									fbr[edgeii_E - 1].real = (fbr[edgeii_E - 1].real - 2* Amplitude_wave[i] * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].imag() * el[ith].face[nn].Area * weight0[nGauss]);
									fbr[edgeii_E - 1].imag = (fbr[edgeii_E - 1].imag + 2 * Amplitude_wave[i] * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].real() * el[ith].face[nn].Area * weight0[nGauss]);
								}
							}
						}
					}
				}

			}
		}
	}

	return 0;
}


int Fbr_Setting22(Element* el, int myid) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	//double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();
	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	fbr = (MKL_Complex16*)malloc(sizeof(MKL_Complex16) * number_unknown_subdomain);
	for (int i = 0; i < number_unknown_subdomain; ++i) {
		fbr[i].real = 0.0;
		fbr[i].imag = 0.0;
	}

	double disx, disy, r;
	double Kwave[3];
	Kwave[0] = 0.0; Kwave[2] = omega * sqrt(epsilon0 * mur); Kwave[1] = 0.0;


	double yc1;
	double A10 = 2 * pi * sqrt(Waveguide_P / omega / mur / w_wg / w_wg / w_wg / h_wg);


	for (int i = 0; i < num_wave; i++) {
		for (int ith = 0; ith < num_element_subdomain; ith++) {

			double zt = (el[ith].node[0].zb[2] + el[ith].node[1].zb[2] + el[ith].node[2].zb[2] + el[ith].node[3].zb[2]) / 4.0;
			//if (el[ith].Material != wave_material) {
			if (zt < wave_zzz) {
				continue;
			}
			Kwave[2] = omega * sqrt(el[ith].epsil * mur);
			for (nn = 0; nn < 4; nn++) {
				op1 = el[ith].face[nn].opp[0]; ofn1 = el[ith].face[nn].opp[1];
				if (ofn1 == (-200 - i)) {
					xrod1 = x_wave[i]; yrod1 = y_wave[i];
					//xrod1 = 6.6e-3; yrod1 = 0;
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
								disx = (rs[0] - xrod1); disy = (rs[1] - yrod1);
								r = sqrt(pow(disx, 2) + pow(disy, 2));
								//value of test basis function
								Li1[0] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].b[node_ii2 - 1];
								Li1[1] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].c[node_ii2 - 1];
								Li1[2] = el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii1 - 1] + el[ith].b[node_ii1 - 1] * rs[0] + el[ith].c[node_ii1 - 1] * rs[1] + el[ith].d[node_ii1 - 1] * rs[2]) * el[ith].d[node_ii2 - 1];

								Li2[0] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].b[node_ii1 - 1];
								Li2[1] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].c[node_ii1 - 1];
								Li2[2] = sign_judge(edgeii_loc, 6) * el[ith].length[edgeii_loc - 1] / (36.0 * el[ith].Ve * el[ith].Ve) * (el[ith].a[node_ii2 - 1] + el[ith].b[node_ii2 - 1] * rs[0] + el[ith].c[node_ii2 - 1] * rs[1] + el[ith].d[node_ii2 - 1] * rs[2]) * el[ith].d[node_ii1 - 1];

								//Incident E has Z component;
								EincX.real((cos(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
								EincX.imag((-sin(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disx / r));
								EincY.real((cos(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
								EincY.imag((-sin(Kwave[2] * rs[2] + phase_wave[i])) * ((1 / r / log(rb_wave / ra_wave)) * disy / r));
								EincZ.real(0.0);
								EincZ.imag(0.0);

								Einc[0].real(EincX.real() / el[ith].eta); Einc[0].imag(EincX.imag() / el[ith].eta);
								Einc[1].real(EincY.real() / el[ith].eta); Einc[1].imag(EincY.imag() / el[ith].eta);
								Einc[2].real(EincZ.real() / el[ith].eta); Einc[2].imag(EincZ.imag() / el[ith].eta);
								ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
								ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
								ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];
								ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
								ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
								ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];
								for (int ii0 = 0; ii0 < 3; ii0++) {
									//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
									//fbr[edgeii_E - 1].imag += 2 * omega* signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].real() + ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
									//fbr[edgeii_E - 1].imag -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].imag());
									//fbr[edgeii_E - 1].real -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].real());
									//fbr[edgeii_E - 1].real+= - 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].imag()) * el[ith].face[nn].Area * weight0[nGauss];
									//fbr[edgeii_E - 1].imag+= + 2.0 * omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
									fbr[edgeii_E - 1].real = (fbr[edgeii_E - 1].real - 2 * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].imag() * el[ith].face[nn].Area * weight0[nGauss]);
									fbr[edgeii_E - 1].imag = (fbr[edgeii_E - 1].imag + 2 * omega * signE_Ni * (Li1[ii0] + Li2[ii0]) * ncrosEcrosn[ii0].real() * el[ith].face[nn].Area * weight0[nGauss]);
								}
							}
						}
					}
				}

			}
		}
	}

	return 0;
}

int Fbr_Setting_LP(Element* el,int myid) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	//double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();
	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	fbr = (MKL_Complex16*)malloc(sizeof(MKL_Complex16)* number_unknown_subdomain);
	for (int i = 0; i < number_unknown_subdomain; ++i) {
		fbr[i].real = 0.0;
		fbr[i].imag = 0.0;
	}
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
							EincX.real(0.0); EincX.imag(0.0);
							EincY.real(0.0); EincY.imag(0.0);
							EincZ.real(-factor / Hsubstrate); EincZ.imag(0.0);//E = V0/d

							Einc[0].real(EincX.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[0].imag(EincX.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							Einc[1].real(EincY.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[1].imag(EincY.imag() / (Z0 * Wmicrostrip / Hsubstrate));
							Einc[2].real(EincZ.real() / (Z0 * Wmicrostrip / Hsubstrate)); Einc[2].imag(EincZ.imag() / (Z0 * Wmicrostrip / Hsubstrate));

							ncrosE[0] = nVec[1] * Einc[2] - nVec[2] * Einc[1];
							ncrosE[1] = nVec[2] * Einc[0] - nVec[0] * Einc[2];
							ncrosE[2] = nVec[0] * Einc[1] - nVec[1] * Einc[0];

							ncrosEcrosn[0] = ncrosE[1] * nVec[2] - ncrosE[2] * nVec[1];
							ncrosEcrosn[1] = ncrosE[2] * nVec[0] - ncrosE[0] * nVec[2];
							ncrosEcrosn[2] = ncrosE[0] * nVec[1] - ncrosE[1] * nVec[0];

							//incident H has x component
							//HincX.real(-factor / thickness);
							HincX.real(0.0);
							HincX.imag(0.0);
							HincY.real(0.0); HincY.imag(0.0);
							HincZ.real(0.0); HincZ.imag(0.0);
							Hinc[0].real(HincX.real() / el[ith].eta); Hinc[0].imag(HincX.imag() / el[ith].eta);
							Hinc[1].real(HincY.real() / el[ith].eta); Hinc[1].imag(HincY.imag() / el[ith].eta);
							Hinc[2].real(HincZ.real() / el[ith].eta); Hinc[2].imag(HincZ.imag() / el[ith].eta);
							ncrosH[0] = nVec[1] * Hinc[2] - nVec[2] * Hinc[1];
							ncrosH[1] = nVec[2] * Hinc[0] - nVec[0] * Hinc[2];
							ncrosH[2] = nVec[0] * Hinc[1] - nVec[1] * Hinc[0];
							for (int ii0 = 0; ii0 < 3; ii0++) {
								//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
								fbr[edgeii_E - 1].imag += 2 * omega* signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].real() + ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
							}
						}
					}
				}
			}

		}
	}


	return 0;
}

int Fbr_Setting_WAVE(Element* el, int myid) {
	int op1, ofn1, edgeii_loc, node_ii1, node_ii2, edgeii_E, signE_Ni, nn, ii, jj, nGauss, nCount, nCount1, node_glb, node1, node2, node3, num_count;
	complex<double> nVec[3], DL[3], V1, V2, I1, I2, J1, J2, widthWc, S11, S12, S21, S22;
	double zb1[3], zb2[3], zb3[3], rs[3], Li1[3], Li2[3], zb0[3], zbx, zby, zbz, Jinc1[3], Jinc2[3];
	int nQuads = 6;
	//double weight0[6] = { {0.22338158},{0.22338158},{0.22338158},{0.10995174},{0.10995174},{0.10995174} };
	//double alpha0[6] = { {0.10810301},{0.44594849},{0.44594849},{0.81684757},{0.09157621},{0.09157621} };
	//double beta0[6] = { {0.44594849},{0.10810301},{0.44594849},{0.09157621},{0.81684757},{0.09157621} };
	//double gamma0[6] = { {0.44594849},{0.44594849},{0.10810301},{0.09157621},{0.09157621},{0.81684757} };
	double weight0[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
	double alpha0[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
	double beta0[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
	double gamma0[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };
	complex <double> Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	complex<double> EincX, EincY, EincZ, HincX, HincY, HincZ;
	//fbr.setZero();
	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	fbr = (MKL_Complex16*)malloc(sizeof(MKL_Complex16) * number_unknown_subdomain);
	for (int i = 0; i < number_unknown_subdomain; ++i) {
		fbr[i].real = 0.0;
		fbr[i].imag = 0.0;
	}
	double yc1;
	double A10=2*pi*sqrt(Waveguide_P/omega/mur/w_wg/ w_wg/ w_wg/h_wg);
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
							EincX.real(0.0); EincX.imag(0.0);
							EincY.real(0.0); EincY.imag(0.0);
							//EincZ.real(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * epsil[el] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
							//EincZ.imag(omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-cos((sqrt(omega * omega * epsil[el] * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - rs[0]))));
							EincZ.real(A10 *omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (-sin((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
							EincZ.imag(A10 *omega * mur * w_wg / pi * sin(pi * (rs[1] - y_o) / w_wg) * (cos((sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg))) * (rs[0] - x_o))));
							Einc[0].real(0.0);
							Einc[0].imag(0.0);
							Einc[1].real(0.0);
							Einc[1].imag(0.0);
							Einc[2].real(2 * (sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg)) / mur) * EincZ.imag());
							Einc[2].imag(-2 * (sqrt(omega * omega * el[ith].epsil * mur - (pi / w_wg) * (pi / w_wg)) / mur) * EincZ.real());
							for (int ii0 = 0; ii0 < 3; ii0++) {
								//fbr[edgeii_E - 1].real(fbr[edgeii_E - 1].real() - omega * signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].imag() + ncrosEcrosn[ii0].imag()) * Area[el][nn] * weight0[nGauss]);
								//fbr[edgeii_E - 1].imag += 2 * omega* signE_Ni * 1.0 * (Li1[ii0] + Li2[ii0]) * (ncrosH[ii0].real() + ncrosEcrosn[ii0].real()) * el[ith].face[nn].Area * weight0[nGauss];
								fbr[edgeii_E - 1].imag -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].imag());
								fbr[edgeii_E - 1].real -= (signE_Ni * (Li1[ii0] + Li2[ii0]) * el[ith].face[nn].Area * weight0[nGauss]) * (Einc[ii0].real());
							}
						}
					}
				}
			}

		}
	}
	cout <<"omega * mur * w_wg / pi  = "<< omega * mur * w_wg / pi << endl;

	return 0;
}

int Solve_E_H_boundary_before(Element* el, int myid) {

	Fbr_Setting(el, myid);
	//cout <<"matsize"<< matsize << endl;
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega + (sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].real + sqrt(omega) * acoo_ome_sqrt[i].real;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega + (sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].imag + sqrt(omega) * acoo_ome_sqrt[i].imag;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega + (sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].real + sqrt(omega) * m_ome_sqrt[i].real;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega + (sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].imag + sqrt(omega) * m_ome_sqrt[i].imag;
	}
	//for (int i = 0; i < matsize; i++) {
	//	acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
	//	acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
	//	m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
	//	m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	//}

	int nnz_temp= num_nzero_Pmatrix[myid];



	//if (myid == 0) {
	//	ofstream ofs0201("matrix_A11_2.csv", ios::app);
	//	for (int i = 0; i < nnz_temp; i++) {
	//		ofs0201 << rowind[i] << ',' << colind[i] << ','<< acoo[i].imag << ',' << acoo[i].real << endl;
	//	}
	//	ofstream ofs0202("matrix_A11_nnz_unknown.csv", ios::app);
	//	ofs0202 << nnz_temp << ',' << num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] << ',' << Edge_Boundary[myid] << endl;
	//}




	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	//cout << "myid is " << myid << endl;
	//cout << "matsize" << matsize << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	cout << "Edge_Boundary["<<myid<<"] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);

	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
	cout << " my id is " << myid << endl;
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
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}
	
	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}


	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	if (myid == 0) {
		//cout << "time is " << time << endl;
		ofstream ofs1119("data_DG_RTC_patch_method1.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< 1819 << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}

	long long size_face = 0;
	for (int j = 0; j < num_domain; j++) {
		size_face += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}

	//cout << "size_face = " << size_face << endl;








	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}



	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;



	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	time_t start, end;

	start = time(NULL);

	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;


		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;



		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}

		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		//clock_t tt1, tt2;
		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
		iparm[1] = 3; iparm[24] = 2;
		
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);

		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
		// calulate full solution Um

	}


	end = time(NULL);

	double time = (double)(end - start);

	//int num_opp_element_boundary = 0;
	//for (int j = 0; j < num_domain; j++) {
	//	num_opp_element_boundary += num_opp_element[j];
	//}
	//

	if (myid == 0) {
		cout << "time is " << time << endl;
		ofstream ofs1119("data_DG_RTC_awace_method1.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}
	

	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];

	
	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
		//ofstream ofs732("Solution_face.csv", ios::app);
		//for (int i = 0; i < num_unk_boundary; i++)
		//{
		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
		//}

	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;

//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary(Element* el, int myid) {

	Fbr_Setting(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].real + sqrt(omega) * acoo_ome_sqrt[i].real;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].imag + sqrt(omega) * acoo_ome_sqrt[i].imag;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].real + sqrt(omega) * m_ome_sqrt[i].real;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].imag + sqrt(omega) * m_ome_sqrt[i].imag;
	}

	int nnz_temp = num_nzero_Pmatrix[myid];



	time_t start, end;

	start = time(NULL);

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}

	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	

	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
	cout << " my id is " << myid << endl;
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
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}


	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	long long size_face = 0;
	for (int j = 0; j < num_domain; j++) {
		size_face += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}

	//cout << "size_face = " << size_face << endl;



	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}



	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;



	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;


		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}
		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		iparm[1] = 3; iparm[24] = 2;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
		// calulate full solution Um

	}


	end = time(NULL);

	

	//int num_opp_element_boundary = 0;
	//for (int j = 0; j < num_domain; j++) {
	//	num_opp_element_boundary += num_opp_element[j];
	//}
	//

	omp_set_num_threads(process_threads);
	MPI_Barrier(MPI_COMM_WORLD);
	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];
	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	end_out = time(NULL);

	double time_out = (double)(end_out - start_out);
	double time = (double)(end - start);
	if (myid == 0) {
		cout << "time is " << time << endl;
		int num_un_total = 0; 
		num_un_total = Accumulated_unknowns[num_domain - 1] + num_unknown_subdomain[num_domain - 1][0] + num_unknown_subdomain[num_domain - 1][1];
		ofstream ofs1119("data_DG_RTC_waveguide_method1.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_un_total << ',' << num_unk_boundary << ',' << time_out << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol0.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 1) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol1.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 2) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol2.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 3) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol3.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	free(fbr); delete[] unKnownUi;
	//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary2(Element* el, int myid) {

	Fbr_Setting2(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].real + sqrt(omega) * acoo_ome_sqrt[i].real;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * acoo_wave[i].imag + sqrt(omega) * acoo_ome_sqrt[i].imag;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].real + sqrt(omega) * m_ome_sqrt[i].real;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;// +(sqrt(omega * omega * epsilon0 * mur - (pi / w_wg) * (pi / w_wg)) / mur) * m_wave[i].imag + sqrt(omega) * m_ome_sqrt[i].imag;
	}

	int nnz_temp = num_nzero_Pmatrix[myid];



	time_t start, end;

	start = time(NULL);

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}

	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);


	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
	cout << " my id is " << myid << endl;
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
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}


	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	long long size_face = 0;
	for (int j = 0; j < num_domain; j++) {
		size_face += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}

	//cout << "size_face = " << size_face << endl;



	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}



	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;



	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;


		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}
		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		iparm[1] = 3; iparm[24] = 2;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
		// calulate full solution Um

	}


	end = time(NULL);



	//int num_opp_element_boundary = 0;
	//for (int j = 0; j < num_domain; j++) {
	//	num_opp_element_boundary += num_opp_element[j];
	//}
	//

	omp_set_num_threads(process_threads);
	MPI_Barrier(MPI_COMM_WORLD);
	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];
	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	end_out = time(NULL);

	double time_out = (double)(end_out - start_out);
	double time = (double)(end - start);
	if (myid == 0) {
		cout << "time is " << time << endl;
		int num_un_total = 0;
		num_un_total = Accumulated_unknowns[num_domain - 1] + num_unknown_subdomain[num_domain - 1][0] + num_unknown_subdomain[num_domain - 1][1];
		ofstream ofs1119("data_DG_RTC_waveguide_method1.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_un_total << ',' << num_unk_boundary << ',' << time_out << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol0.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 1) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol1.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 2) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol2.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	if (myid == 3) {
		if (Freq == 1.5e9) {
			ofstream ofs08100("Xsol3.csv");
			for (int j = 0; j < num_unknown_this_process; j++) {
				ofs08100 << Xsol[j].imag() << ',' << Xsol[j].real() << endl;
			}
		}
	}
	free(fbr); delete[] unKnownUi;
	//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary_E_T(Element* el, int myid) {

	Fbr_Setting(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}
	int nnz_temp = num_nzero_Pmatrix[myid];
	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure_E_T(myid);
	end_FETI = time(NULL);
	double time_FETI = (double)(end_FETI - start_FETI);
	//cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);

	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
	int* nnz_ZC_process = new int[num_domain];
	for (int j = 0; j < num_domain; j++) {
		nnz_ZC_process[j] = 0;
	}
	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	int face_mat_size = 0;
	for (int j = 0; j < num_domain; j++) {
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}
	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}
	long long size_face = 0;
	for (int j = 0; j < num_domain; j++) {
		size_face += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}

	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}
	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];
	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	time_t start, end;

	start = time(NULL);

	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);
		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;
		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}

		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		iparm[1] = 3; iparm[24] = 2;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;

	}
	end = time(NULL);
	double time = (double)(end - start);

	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary_iter(Element* el, int myid) {

	Fbr_Setting(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure_iter(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time_FETI is  " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);


	int matsize_of_C = 0;
	int* nnz_of_C = new int[num_domain];

	for (int i = 0; i < num_domain; i++) {
		nnz_of_C[i] = 0;
	}

	for (int i = 0; i < num_domain; i++) {
		nnz_of_C[i] = matsize_of_C;
		if (nnz_c[myid][i] != 0 && i != myid) {

			matsize_of_C += Edge_Boundary[myid] * Edge_Boundary[i];
		}
	}
	cout << "matsize_of_C = " << matsize_of_C << endl;
	MKL_Complex16* CC = new MKL_Complex16[matsize_of_C];
	for (int i = 0; i < matsize_of_C; i++)
	{
		CC[i].imag = 0;
		CC[i].real = 0;
	}
	Matrix_Generator_CY2(matsize_of_C, nnz_of_C, myid, CC);
	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	//cout << "myid is " << myid << "  num_unk_boundary = " << num_unk_boundary << endl;

	int size_fbb = Edge_Boundary[myid];


	//complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]]();
	//complex<double>* fbb_total = new complex<double>[num_unk_boundary]();
	complex<double>* fbb_temp = (complex<double>*)malloc(Edge_Boundary[myid] * sizeof(complex<double>));
	complex<double>* fbb_total = (complex<double>*)malloc(num_unk_boundary * sizeof(complex<double>));
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];
	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Bcast(fbb_total, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	unKnownUi = new MKL_Complex16[num_unk_boundary];


	time_t start, end;
	start = time(NULL);

	/* Jacobian iterative process  */
	ofstream ofs231("error_of_iteration.csv", ios::app);
	int num_cycle1 = 200;
	vector<float> v9(num_cycle1, 0.0);
	float bb = 1.0;
	int k = 0;
	while (k < num_cycle1) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = fbb_total[j].real();
			unKnownUi[j].imag = fbb_total[j].imag();
		}
		for (int j = 0; j < Edge_Boundary[myid]; j++) {
			fbb_temp[j].real(fbb[j].real);
			fbb_temp[j].imag(fbb[j].imag);
		}
		Matrix_Generator_JSIM6(nnz_of_C, myid, CC, fbb_temp, fbb_total);
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(fbb_total, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		float abs_uk = 0.0;
		for (int i = 0; i < num_unk_boundary; i++) {
			v9[k] += (unKnownUi[i].real - fbb_total[i].real()) * (unKnownUi[i].real - fbb_total[i].real()) + (unKnownUi[i].imag - fbb_total[i].imag()) * (unKnownUi[i].imag - fbb_total[i].imag());
			abs_uk += unKnownUi[i].real * unKnownUi[i].real + unKnownUi[i].imag * unKnownUi[i].imag;
		}
		v9[k] = v9[k] / abs_uk;
		bb = v9[k];
		k++;
		//cout << v9[k] << endl;

	}
	for (int j = 0; j < num_unk_boundary; j++) {
		unKnownUi[j].real = fbb_total[j].real();
		unKnownUi[j].imag = fbb_total[j].imag();
	}

	/*    */
	//ofstream ofs8026("fbb_802.csv", ios::app);
	//for (int i = 0; i < num_unk_boundary; i++) {
	//	ofs8026 << fbb[i].imag << ',' << fbb[i].real << endl;
	//}
	/*    */
	cout.precision(16);
	if (myid == 0) {
		for (int j = 0; j < k; j++)
		{
			ofs231 << j << ',' << v9[j] << endl;
		}
	}
	ofstream ofs241("error_of_iteration_last_one.csv", ios::app);
	if (myid == 0) {
		ofs241 << num_cycle1 << ',' << v9[num_cycle1 - 1] << endl;
	}



	/* Jacobian iterative process  */


	end = time(NULL);

	double time_face = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	if (myid == 0) {
		cout << "time_out is " << time_out << endl;
		cout << "time_face is " << time_face << endl;
		ofstream ofs1119("data_DG_RTC_Yagi_iteration.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_FETI << "s" << ',' \
			<< time_face << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << matsize_of_C << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
	}


	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	start = time(NULL);
	Solver_mpi(myid);
	end = time(NULL);

	time_face = (double)(end - start);

	if (myid == 0) {
		cout << "Solver_mpi time is " << time_face << endl;
	}


	free(fbr); delete[] unKnownUi;

	//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary_iter_long(Element* el, int myid) {

	Fbr_Setting(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	//ofstream ofs0122("Awace_A_id0.csv", ios::app);
	//ofstream ofs01221("Awace_A_id0_row.csv", ios::app);
	//ofstream ofs01222("Awace_A_id0_col.csv", ios::app);
	//if (myid == 0) {
	//	for (int i = 0; i < matsize; i++) {
	//		ofs0122 << acoo[i].imag << ',' << acoo[i].real << endl;
	//		ofs01221 << rowind[i] << endl;
	//		ofs01222 << colind[i] << endl;
	//	}
	//	cout << " Matrix output finished! " << endl;
	//}

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure_iter(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time_FETI is  " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);

	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);


	long long matsize_of_C = 0;
	long long* nnz_of_C = new long long[num_domain];

	for (int i = 0; i < num_domain; i++) {
		nnz_of_C[i] = 0;
	}

	for (int i = 0; i < num_domain; i++) {
		nnz_of_C[i] = matsize_of_C;
		if (nnz_c[myid][i] != 0 && i != myid) {
			matsize_of_C += (Edge_Boundary[myid] * Edge_Boundary[i]);
		}
	}
	cout << "long long matsize_of_C = " << matsize_of_C << endl;
	MKL_Complex16* CC = new MKL_Complex16[matsize_of_C];
	for (long long i = 0; i < matsize_of_C; i++)
	{
		CC[i].imag = 0;
		CC[i].real = 0;
	}

	//cout << "myid is " << myid << "  test1" << endl;

	Matrix_Generator_CY_long(matsize_of_C, nnz_of_C, myid, CC);

	//cout << "myid is " << myid << "  test2" << endl;
	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	//cout << "myid is " << myid << "  num_unk_boundary = " << num_unk_boundary << endl;

	int size_fbb = Edge_Boundary[myid];


	//complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]]();
	//complex<double>* fbb_total = new complex<double>[num_unk_boundary]();
	complex<double>* fbb_temp = (complex<double>*)malloc(Edge_Boundary[myid] * sizeof(complex<double>));
	complex<double>* fbb_total = (complex<double>*)malloc(num_unk_boundary * sizeof(complex<double>));
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];
	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	MPI_Bcast(fbb_total, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	unKnownUi = new MKL_Complex16[num_unk_boundary];


	time_t start, end;
	start = time(NULL);

	/* Jacobian iterative process  */
	ofstream ofs231("error_of_iteration_long.csv", ios::app);
	int num_cycle1 = 500;
	vector<float> v9(num_cycle1, 0.0);
	float bb = 1.0;
	int k = 0;
	while (k < num_cycle1&&bb>1e-4) {
		//cout << "myid is " << myid << "  test3" << endl;
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = fbb_total[j].real();
			unKnownUi[j].imag = fbb_total[j].imag();
		}
		for (int j = 0; j < Edge_Boundary[myid]; j++) {
			fbb_temp[j].real(fbb[j].real);
			fbb_temp[j].imag(fbb[j].imag);
		}
		Matrix_Generator_JSIM_long(nnz_of_C, myid, CC, fbb_temp, fbb_total);
		MPI_Barrier(MPI_COMM_WORLD);
		//cout << "myid is " << myid << "  test4" << endl;
		MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(fbb_total, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		float abs_uk = 0.0;
		for (int i = 0; i < num_unk_boundary; i++) {
			v9[k] += (unKnownUi[i].real - fbb_total[i].real()) * (unKnownUi[i].real - fbb_total[i].real()) + (unKnownUi[i].imag - fbb_total[i].imag()) * (unKnownUi[i].imag - fbb_total[i].imag());
			abs_uk += unKnownUi[i].real * unKnownUi[i].real + unKnownUi[i].imag * unKnownUi[i].imag;
		}
		v9[k] = v9[k] / abs_uk;
		bb = v9[k];
		k++;
		//cout << v9[k] << endl;

	}
	for (int j = 0; j < num_unk_boundary; j++) {
		unKnownUi[j].real = fbb_total[j].real();
		unKnownUi[j].imag = fbb_total[j].imag();
	}

	/*    */
	//ofstream ofs8026("fbb_802.csv", ios::app);
	//for (int i = 0; i < num_unk_boundary; i++) {
	//	ofs8026 << fbb[i].imag << ',' << fbb[i].real << endl;
	//}
	/*    */
	cout.precision(16);
	if (myid == 0) {
		for (int j = 0; j < k; j++)
		{
			ofs231 << j << ',' << v9[j] << endl;
		}
	}
	ofstream ofs241("error_of_iteration_last_one_long.csv", ios::app);
	if (myid == 0) {
		ofs241 << num_cycle1 << ',' << v9[num_cycle1 - 1] << endl;
	}



	/* Jacobian iterative process  */


	end = time(NULL);

	double time_face = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	if (myid == 0) {
		cout << "time_out is " << time_out << endl;
		cout << "time_face is " << time_face << endl;
		ofstream ofs1119("data_DG_RTC_awace_iteration.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time_face << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << matsize_of_C << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
	}


	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic; delete[]CC;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	start = time(NULL);
	Solver_mpi(myid);
	end = time(NULL);

	time_face = (double)(end - start);

	if (myid == 0) {
		cout << "Solver_mpi time is " << time_face << endl;
	}


	free(fbr); delete[] unKnownUi;

	//	free(fbr); delete[] unKnownUi;
	return  0;

}

//int Solve_E_H_boundary_Eigen(Element* el, int myid) {
//
//	Fbr_Setting(el, myid);
//	//cout <<"matsize"<< matsize << endl;
//#pragma omp  parallel for
//	for (int i = 0; i < matsize; i++) {
//		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
//		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
//		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
//		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
//	}
//
//	r_dm = new int[num_domain];
//	r_dm[0] = 0;
//	for (int i = 1; i < num_domain; ++i) {
//		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
//	}
//	//cout << "myid is " << myid << endl;
//	//cout << "matsize" << matsize << endl;
//	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
//	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;
//
//	time_t start_out, end_out;
//	start_out = time(NULL);
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	time_t start_FETI, end_FETI;
//	start_FETI = time(NULL);
//	FETI_like_procedure(myid);
//	end_FETI = time(NULL);
//
//	double time_FETI = (double)(end_FETI - start_FETI);
//	cout << "myid is " << myid << "time is " << time_FETI << endl;
//	MPI_Barrier(MPI_COMM_WORLD);
//	end_out = time(NULL);
//	double time_out = (double)(end_out - start_out);
//
//	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
//	complex<double>* m_final = new complex<double>[mat_size];
//	int* rowm = new int[mat_size], * colm = new int[mat_size];
//#	pragma omp parallel for 
//	for (int i = 0; i < mat_size; i++) {
//		m_final[i].real(TriList[i].value().real());
//		m_final[i].imag(TriList[i].value().imag());
//		rowm[i] = TriList[i].row() + 1;
//		colm[i] = TriList[i].col() + 1;
//	}
//	TriList.clear();
//	cout << " my id is " << myid << endl;
//	cout << " mat_size is " << mat_size << endl;
//	int* nnz_ZC_process = new int[num_domain];
//	for (int j = 0; j < num_domain; j++) {
//		nnz_ZC_process[j] = 0;
//	}
//	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);
//
//	//cout << "after,  my id is " << myid << endl;
//	int face_mat_size = 0;
//	for (int j = 0; j < num_domain; j++) {
//		//cout << nnz_ZC_process[j] << endl;
//		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
//	}
//	if (myid == 0) {
//		cout << "face_mat_size = " << face_mat_size << endl;
//	}
//
//	int* address_offset_mat = new int[num_domain];
//	address_offset_mat[0] = 0;
//	for (int j = 1; j < num_domain; j++) {
//		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
//	}
//
//	int num_unk_boundary = 0;
//	for (int j = 0; j < num_domain; j++) {
//		num_unk_boundary += Edge_Boundary[j];
//	}
//	int* rowm_total = nullptr;
//	int* colm_total = nullptr;
//	complex<double>* m_final_total = nullptr;
//	complex<double>* fbb_total = nullptr;
//	if (myid == 0) {
//		m_final_total = new complex<double>[face_mat_size];
//		rowm_total = new int[face_mat_size];
//		colm_total = new int[face_mat_size];
//		fbb_total = new complex<double>[num_unk_boundary];
//		for (int j = 0; j < num_unk_boundary; j++) {
//			fbb_total[j].imag(0.0);
//			fbb_total[j].real(0.0);
//		}
//
//
//
//	}
//	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
//	for (int j = 0; j < Edge_Boundary[myid]; j++) {
//		fbb_temp[j].real(fbb[j].real);
//		fbb_temp[j].imag(fbb[j].imag);
//	}
//	int* address_offset = new int[num_domain];
//
//	address_offset[0] = 0;
//	for (int j = 1; j < num_domain; j++) {
//		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
//	}
//
//	//cout << "test1111" << endl;
//	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//	//cout << "test2222" << endl;
//	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
//	//cout << "test3333" << endl;
//	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
//	//cout << "test4444" << endl;
//	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//	//cout << "test5555" << endl;
//
//
//
//	unKnownUi = new MKL_Complex16[num_unk_boundary];
//	MKL_Complex16* m_final_mkl;
//
//	complex<double>* xbb_total = nullptr;
//	xbb_total = new complex<double>[num_unk_boundary];
//
//
//	time_t start, end;
//
//	start = time(NULL);
//
//	omp_set_num_threads(64);
//
//	if (myid == 0) {
//
//		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
//			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
//			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
//			m_final_total[j].real(1.0);
//		}
//		vector<Eigen::Triplet<complex<double>>> Tri_tmp(face_mat_size);
//
//		for (int j = 0; j < face_mat_size; j++) {
//			Tri_tmp[j] = Triplet<complex<double>>(rowm_total[j], colm_total[j], m_final_total[j]);
//		}
//
//		Eigen::SparseMatrix <complex<double>> M_face(num_unk_boundary, num_unk_boundary);
//		M_face.setFromTriplets(Tri_tmp.begin(), Tri_tmp.end());
//
//		SparseLU<SparseMatrix<complex<double>>, COLAMDOrdering<int> >   solver;
//		// fill A and b;
//		// Compute the ordering permutation vector from the structural pattern of A
//		solver.analyzePattern(M_face);
//		// Compute the numerical factorization 
//		solver.factorize(M_face);
//		//Use the factors to solve the linear system 
//		xbb_total = solver.solve(fbb_total);
//
//		//Eigen::solver.compute(M_face);
//	/*	xbb_total = M_face.solve(fbb_total);
//		M_face.s*/
//
//
//
//
//		m_final_mkl = new MKL_Complex16[face_mat_size];
//		//#pragma omp  parallel for
//		for (int j = 0; j < face_mat_size; j++) {
//			m_final_mkl[j].real = m_final_total[j].real();
//			m_final_mkl[j].imag = m_final_total[j].imag();
//		}
//		delete[] m_final_total;
//		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
//		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
//		int info; int nnz_temp = face_mat_size;
//		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
//		//sparse_index_base_t indexM;
//		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);
//
//		delete[] m_final_mkl;
//		delete[] rowm_total;
//		delete[] colm_total;
//
//
//		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
//		//ofstream ofs11141("fbb_total.csv", ios::app);
//		for (int j = 0; j < num_unk_boundary; j++) {
//			fbb_total_mkl[j].real = fbb_total[j].real();
//			fbb_total_mkl[j].imag = fbb_total[j].imag();
//			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
//		}
//		delete[]  fbb_total;
//
//		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
//		MKL_INT iparm[64];
//		int pt[64];
//		MKL_INT* perm2 = nullptr;
//		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
//
//
//
//		for (int el = 0; el < 64; el++) {
//			pt[el] = 0; iparm[el] = 0;
//		}
//
//		iparm[0] = 1;
//		iparm[1] = 3;
//		iparm[2] = 16;
//		iparm[9] = 13;
//		iparm[10] = 1;
//		iparm[12] = 1;
//		iparm[17] = -1;
//		iparm[18] = -1;
//		MKL_INT  nrhs = 1;
//		//clock_t tt1, tt2;
//		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
//		iparm[1] = 3; iparm[24] = 2;
//
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
//		if (error != 0)cout << "error for boundary equation   " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
//
//		//tt2 = clock();
//		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
//		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
//		delete[]fbb_total_mkl;
//		// calulate full solution Um
//
//	}
//
//
//	end = time(NULL);
//
//	double time = (double)(end - start);
//
//	int num_opp_element_boundary = 0;
//	for (int j = 0; j < num_domain; j++) {
//		num_opp_element_boundary += num_opp_element[j];
//	}
//
//
//	if (myid == 0) {
//		cout << "time is " << time << endl;
//		ofstream ofs1119("data_DG_RTC_SIW_method1.csv", ios::app);
//		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
//			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
//			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
//
//	}
//
//
//	omp_set_num_threads(process_threads);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];
//
//
//	if (myid == 0) {
//		for (int j = 0; j < num_unk_boundary; j++) {
//			unKnownUi_mpic[j].real(unKnownUi[j].real);
//			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
//		}
//		cout << "num_unk_boundary = " << num_unk_boundary << endl;
//		//ofstream ofs732("Solution_face.csv", ios::app);
//		//for (int i = 0; i < num_unk_boundary; i++)
//		//{
//		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
//		//}
//
//	}
//	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//	if (myid != 0) {
//		for (int j = 0; j < num_unk_boundary; j++) {
//			unKnownUi[j].real = unKnownUi_mpic[j].real();
//			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
//		}
//	}
//	delete[]unKnownUi_mpic;
//	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
//	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
//	Solver_mpi(myid);
//	free(fbr); delete[] unKnownUi;
//
//	//	free(fbr); delete[] unKnownUi;
//	return  0;
//
//}

//int Solve_E_H_boundary_bicg(Element* el, int myid) {
//
//	Fbr_Setting(el, myid);
//	//cout <<"matsize"<< matsize << endl;
//#pragma omp  parallel for
//	for (int i = 0; i < matsize; i++) {
//		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
//		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
//		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
//		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
//	}
//
//	r_dm = new int[num_domain];
//	r_dm[0] = 0;
//	for (int i = 1; i < num_domain; ++i) {
//		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
//	}
//	//cout << "myid is " << myid << endl;
//	//cout << "matsize" << matsize << endl;
//	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
//	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;
//
//	time_t start_out, end_out;
//	start_out = time(NULL);
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	time_t start_FETI, end_FETI;
//	start_FETI = time(NULL);
//	FETI_like_procedure(myid);
//	end_FETI = time(NULL);
//
//	double time_FETI = (double)(end_FETI - start_FETI);
//	cout << "myid is " << myid << "time is " << time_FETI << endl;
//	MPI_Barrier(MPI_COMM_WORLD);
//	end_out = time(NULL);
//	double time_out = (double)(end_out - start_out);
//
//	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
//	int unknown_boundary = Edge_Boundary[myid];
//	complex<double>* m_final = new complex<double>[mat_size+ unknown_boundary];
//	int* rowm = new int[mat_size+ unknown_boundary], * colm = new int[mat_size+ unknown_boundary];
//#	pragma omp parallel for 
//	for (int i = 0; i < mat_size; i++) {
//		m_final[i].real(TriList[i].value().real());
//		m_final[i].imag(TriList[i].value().imag());
//		rowm[i] = TriList[i].row() + 1;
//		colm[i] = TriList[i].col() + 1;
//	}
//	TriList.clear();
//
//	for (int i = mat_size; i < mat_size+ unknown_boundary; i++) {
//		m_final[i].real(1.0);
//		m_final[i].imag(0.0);
//		rowm[i] = i- mat_size+1;
//		colm[i] = i - mat_size + 1+ r_dm[i];
//	}
//
//	int mat_size_all = mat_size + unknown_boundary;
//
//	int num_unk_boundary = 0;
//	for (int j = 0; j < num_domain; j++) {
//		num_unk_boundary += Edge_Boundary[j];
//	}
//	int* address_offset = new int[num_domain];
//
//	address_offset[0] = 0;
//	for (int j = 1; j < num_domain; j++) {
//		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
//	}
//	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
//	for (int j = 0; j < Edge_Boundary[myid]; j++) {
//		fbb_temp[j].real(fbb[j].real);
//		fbb_temp[j].imag(fbb[j].imag);
//	}
//	complex<double>* fbb_total = new complex<double>[num_unk_boundary];
//	for (int j = 0; j < num_unk_boundary; j++) {
//		fbb_total[j].imag(0.0);
//		fbb_total[j].real(0.0);
//	}
//	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//	MPI_Bcast(fbb_total, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//	double* rho = new double[num_domain];
//	double* alpha = new double[num_domain];
//	double* beta = new double[num_domain];
//	double rho_i = 0;
//	for (int i = 0; i < Edge_Boundary[myid]; i++) {
//		rho_i += fbb_temp[i].imag() * fbb_temp[i].imag() + fbb_temp[i].real() * fbb_temp[i].real();
//	}
//
//
//	MPI_Gather(&rho_i, 1, MPI_DOUBLE, rho, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	MPI_Bcast(rho, num_domain, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//	rho_i = 0;
//
//	for (int i = 0; i < num_domain; i++) {
//		rho_i += rho[i];
//	}
//
//	int num_bicg = 100;
//	double flag_error = 1e-5;
//
//	complex<double>* p_i = new complex<double>[num_unk_boundary];
//	for (int j = 0; j < num_unk_boundary; j++) {
//		p_i[j].imag(fbb_total[j].imag());
//		p_i[j].real(fbb_total[j].real());
//	}
//
//	for (int i = 1; i < num_bicg; i++) {
//
//		for (int j = 0; j < mat_size_all; j++) {
//			int row_temp = rowm[j];
//			int col_temp = colm[j];
//
//		}
//
//
//
//
//
//
//	}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	cout << " my id is " << myid << endl;
//	cout << " mat_size is " << mat_size << endl;
//	int* nnz_ZC_process = new int[num_domain];
//	for (int j = 0; j < num_domain; j++) {
//		nnz_ZC_process[j] = 0;
//	}
//	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);
//
//	//cout << "after,  my id is " << myid << endl;
//	int face_mat_size = 0;
//	for (int j = 0; j < num_domain; j++) {
//		//cout << nnz_ZC_process[j] << endl;
//		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
//	}
//	if (myid == 0) {
//		cout << "face_mat_size = " << face_mat_size << endl;
//	}
//
//	int* address_offset_mat = new int[num_domain];
//	address_offset_mat[0] = 0;
//	for (int j = 1; j < num_domain; j++) {
//		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
//	}
//
//	int num_unk_boundary = 0;
//	for (int j = 0; j < num_domain; j++) {
//		num_unk_boundary += Edge_Boundary[j];
//	}
//	int* rowm_total = nullptr;
//	int* colm_total = nullptr;
//	complex<double>* m_final_total = nullptr;
//	complex<double>* fbb_total = nullptr;
//	if (myid == 0) {
//		m_final_total = new complex<double>[face_mat_size];
//		rowm_total = new int[face_mat_size];
//		colm_total = new int[face_mat_size];
//		fbb_total = new complex<double>[num_unk_boundary];
//		for (int j = 0; j < num_unk_boundary; j++) {
//			fbb_total[j].imag(0.0);
//			fbb_total[j].real(0.0);
//		}
//
//
//
//	}
//	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
//	for (int j = 0; j < Edge_Boundary[myid]; j++) {
//		fbb_temp[j].real(fbb[j].real);
//		fbb_temp[j].imag(fbb[j].imag);
//	}
//	int* address_offset = new int[num_domain];
//
//	address_offset[0] = 0;
//	for (int j = 1; j < num_domain; j++) {
//		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
//	}
//
//	//cout << "test1111" << endl;
//	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//	//cout << "test2222" << endl;
//	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
//	//cout << "test3333" << endl;
//	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
//	//cout << "test4444" << endl;
//	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//	//cout << "test5555" << endl;
//
//
//
//	unKnownUi = new MKL_Complex16[num_unk_boundary];
//	MKL_Complex16* m_final_mkl;
//
//
//
//	time_t start, end;
//
//	start = time(NULL);
//
//	omp_set_num_threads(64);
//
//	if (myid == 0) {
//
//		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
//			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
//			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
//			m_final_total[j].real(1.0);
//		}
//		m_final_mkl = new MKL_Complex16[face_mat_size];
//		//#pragma omp  parallel for
//		for (int j = 0; j < face_mat_size; j++) {
//			m_final_mkl[j].real = m_final_total[j].real();
//			m_final_mkl[j].imag = m_final_total[j].imag();
//		}
//
//		delete[] m_final_total;
//		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
//		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
//		int info; int nnz_temp = face_mat_size;
//		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
//		//sparse_index_base_t indexM;
//		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);
//
//		delete[] m_final_mkl;
//		delete[] rowm_total;
//		delete[] colm_total;
//
//
//		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
//		//ofstream ofs11141("fbb_total.csv", ios::app);
//		for (int j = 0; j < num_unk_boundary; j++) {
//			fbb_total_mkl[j].real = fbb_total[j].real();
//			fbb_total_mkl[j].imag = fbb_total[j].imag();
//			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
//		}
//		delete[]  fbb_total;
//
//		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
//		MKL_INT iparm[64];
//		int pt[64];
//		MKL_INT* perm2 = nullptr;
//		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
//
//
//
//		for (int el = 0; el < 64; el++) {
//			pt[el] = 0; iparm[el] = 0;
//		}
//
//		iparm[0] = 1;
//		iparm[1] = 3;
//		iparm[2] = 16;
//		iparm[9] = 13;
//		iparm[10] = 1;
//		iparm[12] = 1;
//		iparm[17] = -1;
//		iparm[18] = -1;
//		MKL_INT  nrhs = 1;
//		//clock_t tt1, tt2;
//		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
//		iparm[1] = 3; iparm[24] = 2;
//
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
//		if (error != 0)cout << "error for boundary equation   " << error << endl;
//		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
//
//		//tt2 = clock();
//		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
//		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
//		delete[]fbb_total_mkl;
//		// calulate full solution Um
//
//	}
//
//
//	end = time(NULL);
//
//	double time = (double)(end - start);
//
//	int num_opp_element_boundary = 0;
//	for (int j = 0; j < num_domain; j++) {
//		num_opp_element_boundary += num_opp_element[j];
//	}
//
//
//	if (myid == 0) {
//		cout << "time is " << time << endl;
//		ofstream ofs1119("data_DG_RTC_SIW_method1.csv", ios::app);
//		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
//			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
//			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;
//
//	}
//
//
//	omp_set_num_threads(process_threads);
//
//	MPI_Barrier(MPI_COMM_WORLD);
//
//	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];
//
//
//	if (myid == 0) {
//		for (int j = 0; j < num_unk_boundary; j++) {
//			unKnownUi_mpic[j].real(unKnownUi[j].real);
//			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
//		}
//		cout << "num_unk_boundary = " << num_unk_boundary << endl;
//		//ofstream ofs732("Solution_face.csv", ios::app);
//		//for (int i = 0; i < num_unk_boundary; i++)
//		//{
//		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
//		//}
//
//	}
//	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
//
//	if (myid != 0) {
//		for (int j = 0; j < num_unk_boundary; j++) {
//			unKnownUi[j].real = unKnownUi_mpic[j].real();
//			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
//		}
//	}
//	delete[]unKnownUi_mpic;
//	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
//	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
//	Solver_mpi(myid);
//	free(fbr); delete[] unKnownUi;
//
//	//	free(fbr); delete[] unKnownUi;
//	return  0;
//
//}

int Solve_E_H_boundary_Iteration(Element* el, int myid) {

	Fbr_Setting(el, myid);
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}


	if (myid == 0) {
		ofstream ofs1119("data_DG_RTC_Awace_iteration.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary  << ',' \
			<< num_nzero_Pmatrix[myid] << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}


	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;


	if (myid == 0) {
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	//for (int i = 0; i < num_unk_boundary; i++) {
	//	unKnownUi[i].imag = unKnownUi_temp[i].imag; unKnownUi[i].real = unKnownUi_temp[i].real; //The initial value is the solution of the previous iteration
	//}
	//for (int i = 0; i < Edge_Boundary[myid]; i++) {
	//	unKnownUi2[i].imag = 0; unKnownUi2[i].real = 0;
	//}
	//for (int i = 0; i < num_unKnown_b; i++) {
	//	unKnownUi[i].imag = 0; unKnownUi[i].real = 1;
	//	unKnownUi2[i].imag = 0; unKnownUi2[i].real = 0;
	//}
	//cout << "1234" << endl;
	int unknown_dm = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	int nnz_dm = num_nzero_Pmatrix[myid];

	int job[8] = { 2, 1, 1, 0, nnz_dm, 0, 0, 0 }; int info;
	mkl_zcsrcoo(job, &unknown_dm, p_pardiso, jp, ip, &nnz_dm, acoo, rowind, colind, &info);
	int cnt_n1 = Edge_Boundary[myid];
	MKL_INT maxfct, mnum, mtype, phase1, phase2, nrhs, msglvl, error;
	MKL_INT iparm[64];
	int pt[64];
	MKL_INT* perm2 = nullptr;
	maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;
	for (int el = 0; el < 64; el++) {
		pt[el] = 0; iparm[el] = 0;
	}
	iparm[0] = 1;
	iparm[1] = 3;
	iparm[9] = 13;
	iparm[10] = 1;
	iparm[12] = 1;
	iparm[17] = -1;
	iparm[18] = -1;
	iparm[24] = 1;

	nrhs = 1;
	//cout << "test2" << endl;

	complex<double>* pro = (complex<double>*)malloc(unknown_dm * sizeof(complex<double>));
	fbb = new MKL_Complex16[Edge_Boundary[myid]];

	MKL_INT phase3, phase4;

	phase3 = 12; phase4 = 33;


	pardiso(pt, &maxfct, &mnum, &mtype, &phase3, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
	pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);


	//pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);// pro 为 A（-1）*fi
	//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &unknown_dm, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, fbr, pro, &error);
	for (int i = unknown_dm - cnt_n1; i < unknown_dm; i++) {
		fbb[i - unknown_dm + cnt_n1].real = pro[i].real();
		fbb[i - unknown_dm + cnt_n1].imag = pro[i].imag();// store boundy fi	
	}
	//free(yii);
	free(pro);


	//cout << "test1" << endl;

	complex<double>* ui = new complex<double>[num_unk_boundary];
	complex<double>* ui_temp = new complex<double>[num_unk_boundary];

	for (int i = 0; i < num_unk_boundary; i++) {
		ui[i].real(0.0);
		ui[i].imag(0.0);
		ui_temp[i].real(0.0);
		ui_temp[i].imag(0.0);
	}

	complex<double>* ui1 = new complex<double>[Edge_Boundary[myid]];

	//Matrix_Generator_CY2();
	//float error = 1e-1;
	//Matrix_ACA_Generator(error);

	fbb1 = new MKL_Complex16[Edge_Boundary[myid]];

	MKL_Complex16* temp_fi;
	MKL_Complex16* temp_fi1;

	time_t start, end;

	start = time(NULL);
	/* Jacobian iterative process  */

	int num_cycle1 = 50;
	vector<float> v9(num_cycle1, 0.0);
	float bb = 1.0;
	int k = 0;
	cout << "start" << endl;
	//for (int k = 0; k < num_cycle1; k++) {
	while (k < num_cycle1) {

		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			fbb1[i].real = fbb[i].real;
			fbb1[i].imag = fbb[i].imag;
		}


		for (int n1 = myid; n1 < myid + 1; n1++) {
			//int nn3 = n3_dm[n1];
			//int cnt_n1 = Edge_Boundary[n1];
			//int mm4 = 0;
			//int nn4 = 0;
#		pragma omp parallel for 
			for (int i = 0; i < Edge_Boundary[myid]; i++) {
				ui1[i].imag(0.0);
				ui1[i].real(0.0);
			}
			int end1 = num_unknown_subdomain[n1][0] + num_unknown_subdomain[n1][1];;
			int start1 = end1 - cnt_n1 + 1;
			int mm2 = 0; int nn2 = 0;
			complex<double>* ui1_temp = new complex<double>[cnt_n1];
			for (int n2 = 0; n2 < num_domain; n2++) {
				if (nnz_c[n1][n2] != 0) {
					mm2 = nn2 + nnz_c[n1][n2];
					if (n1 != n2) {
						int cnt_n2 = Edge_Boundary[n2];
						int end2 = num_unknown_subdomain[n2][0] + num_unknown_subdomain[n2][1];
						int start2 = end2 - cnt_n2 + 1;
						//int matsize = cnt_n1 * cnt_n2;
						long long cnt_n1_long = cnt_n1;
						long long cnt_n2_long = cnt_n2;
						long  long matsize = cnt_n1_long * cnt_n2_long;

						int start_ui = r_dm[n2];
						complex<double>* b_pro = new complex<double>[matsize];
						complex<double>* pro = new complex<double>[cnt_n2];

#		pragma omp parallel for 
						for (int i = 0; i < matsize; i++) {
							b_pro[i].real(0);
							b_pro[i].imag(0);

						}
						for (int i = 0; i < cnt_n2; i++) {
							pro[i].real(ui[start_ui + i].real());
							pro[i].imag(ui[start_ui + i].imag());
						}

						//#		pragma omp parallel for 
						//						for (int i = nn2; i < mm2; i++) {         //RmCmnRn
						//#	   pragma omp critical
						//							{
						//								if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
						//									ui1[mrow[i] - start1].real(ui1[mrow[i] - start1].real() + m[i].real * ui[start_ui + mcol[i] - start2].real() - m[i].imag * ui[start_ui + mcol[i] - start2].imag());
						//									ui1[mrow[i] - start1].imag(ui1[mrow[i] - start1].real() + m[i].real * ui[start_ui + mcol[i] - start2].imag() + m[i].imag * ui[start_ui + mcol[i] - start2].real());
						//								}
						//							}
						//						}


//#		pragma omp parallel for 
						for (int i = nn2; i < mm2; i++) {         //RmCmnRn
							if (mrow[i] >= start1 && mrow[i] <= end1 && mcol[i] >= start2 && mcol[i] <= end2) {
								long long long_temp1 = mrow[i] - start1;
								long long long_temp2 = mcol[i] - start2;
								b_pro[cnt_n2_long * (long_temp1)+long_temp2].real(m[i].real);
								b_pro[cnt_n2_long * (long_temp1)+long_temp2].imag(m[i].imag);
							}
						}

						MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
						cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, 1, cnt_n2, &alpa, b_pro, cnt_n2, pro, 1, &beta, ui1_temp, 1);  //Zm*RmCmnRn

						for (int i = 0; i < cnt_n1; i++) {
							ui1[i].real(ui1[i].real() + ui1_temp[i].real());
							ui1[i].imag(ui1[i].imag() + ui1_temp[i].imag());
						}
						//MKL_Complex16 alpa, beta; alpa.real = 1.0; alpa.imag = 0; beta.real = 0.0; beta.imag = 0;
						//cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cnt_n1, cnt_n2, cnt_n1, &alpa, mat_zm, cnt_n1, b_pro, cnt_n2, &beta, pro, cnt_n2);  //Zm*RmCmnRn
						delete[] b_pro;
						delete[] pro;

					}
					nn2 = mm2;
				}
			}
			delete[] ui1_temp;

		}

		//int cnt_n1 = Edge_Boundary[myid];
		int end1 = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
		int start1 = end1 - cnt_n1;
		temp_fi = new MKL_Complex16[end1];
		temp_fi1 = new MKL_Complex16[end1];
		for (int i = 0; i < end1; i++) {
			temp_fi[i].real = 0.0;
			temp_fi[i].imag = 0.0;
			temp_fi1[i].real = 0.0;
			temp_fi1[i].imag = 0.0;
		}

		
		for (int i = start1; i < end1; i++) {
			temp_fi[i].real = ui1[i - start1].real();
			temp_fi[i].imag = ui1[i - start1].imag();
		}




		nrhs = 1;



		pardiso(pt, &maxfct, &mnum, &mtype, &phase4, &end1, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, temp_fi, temp_fi1, &error);

		////cout << "start pardiso" << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &end1, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, temp_fi, temp_fi1, &error); // x 为 A（-1）
		//if (error != 0)cout << "Pardiso error " << error << endl;
		//pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &end1, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, temp_fi, temp_fi1, &error);   //释放内存（存储内部格式的矩阵）？

		for (int i = 0; i < Edge_Boundary[myid]; i++) {
			ui1[i].real(fbb1[i].real - temp_fi1[i + start1].real);
			ui1[i].imag(fbb1[i].imag - temp_fi1[i + start1].imag);
		}

		MPI_Gatherv(ui1, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, ui, Edge_Boundary, r_dm, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
		MPI_Bcast(ui, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

		//cout << "test7" << endl;
		//Matrix_JSIM_ACA();
		if (myid == 0) {
			float abs_uk = 0.0;
			for (int i = 0; i < num_unk_boundary; i++) {
				v9[k] += (ui[i].real() - ui_temp[i].real()) * (ui[i].real() - ui_temp[i].real()) + (ui[i].imag() - ui_temp[i].imag()) * (ui[i].imag() - ui_temp[i].imag());
				abs_uk += ui[i].real() * ui[i].real() + ui[i].imag() * ui[i].imag();
			}
			v9[k] = v9[k] / abs_uk;
			bb = v9[k];
		}

		MPI_Bcast(&bb, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		v9[k] = bb;
		if (myid == 0) {
			cout << k << endl;
		}
		k++;
		//cout << v9[k] << endl;
		if (k == num_cycle1) {
			pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &end1, p_pardiso, ip, jp, perm2, &nrhs, iparm, &msglvl, temp_fi, temp_fi1, &error);
		}
		for (int j = 0; j < num_unk_boundary; j++) {
			ui_temp[j].real(ui[j].real());
			ui_temp[j].imag(ui[j].imag());
		}

	}
	delete[]ui; delete[]ui1; delete[]temp_fi1; delete[]temp_fi;
	ui = nullptr; ui1 = nullptr; temp_fi1 = nullptr; temp_fi = nullptr;

	end = time(NULL);

	double time = (double)(end - start);

	if (myid == 0) {
		ofstream ofs11270("time_iteraton.csv", ios::app);
		ofs11270 << ',' << num_domain << ',' << num_cycle1 << ',' << time << endl;

	}

	unKnownUi = new MKL_Complex16[num_unk_boundary];
	for (int j = 0; j < num_unk_boundary; j++) {
		unKnownUi[j].real = ui_temp[j].real();
		unKnownUi[j].imag = ui_temp[j].imag();
	}

	ofstream ofs231("vectorfbb2314.csv", ios::app);
	if (myid == 0) {
		for (int j = 0; j < k; j++)
		{
			ofs231 << j << ',' << v9[j] << endl;
		}
	}
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;

	return  0;

}

int Solve_E_H_boundary_scalapack(Element* el, int myid) {

	Fbr_Setting(el, myid);
	//cout <<"matsize"<< matsize << endl;
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	//cout << "myid is " << myid << endl;
	//cout << "matsize" << matsize << endl;
	//cout << "num_unKnown_b is " << num_unKnown_b << endl;
	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure(myid);
	end_FETI = time(NULL);

	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);

	int mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n)
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
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
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}



	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;



	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	time_t start, end;

	start = time(NULL);

	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;


		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;



		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}

		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		//clock_t tt1, tt2;
		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
		iparm[1] = 3; iparm[24] = 2;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);


		//char trans = 'N'; MKL_INT* ipiv = nullptr; MKL_INT info;
		//MKL_INT  nrhs = 1; MKL_INT desca = 1;
		//int* im_B = new int[1 + 1]; im_B[0] = 0; im_B[num_unk_boundary] = 0;

		//int* jm_B = new int[num_unk_boundary];
		//for (int j = 0; j < num_unk_boundary; j++) {
		//	jm_B[j] = j;
		//}
		//pzgesv(&num_unk_boundary, &nrhs, m_pardiso, im, jm
		//	, &num_unk_boundary, ipiv, fbb_total_mkl, im_B, jm_B,
		//	&num_unk_boundary, &info);

		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
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
		ofstream ofs1119("data_DG_RTC_SIW.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}


	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
		//ofstream ofs732("Solution_face.csv", ios::app);
		//for (int i = 0; i < num_unk_boundary; i++)
		//{
		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
		//}

	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;

	//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary_method2(Element* el, int myid) {

	Fbr_Setting(el, myid);
	//cout <<"matsize"<< matsize << endl;
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	int unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];

	int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];

	int nnz_temp_1225 = num_nzero_Pmatrix[myid];

	int* num_nzero_Pmatrix_E_inter = new int[num_domain];


	MKL_Complex16* acoo_temp_1225 = new MKL_Complex16[nnz_temp_1225];



	int* rowind_temp_1225 = new int[nnz_temp_1225];
	int* colind_temp_1225 = new int[nnz_temp_1225];

	int count_nnz_E_inter = 0;
	for (int i = 0; i < nnz_temp_1225; i++) {
		if (rowind[i] <= unknown_E_inter && colind[i] <= unknown_E_inter) {
			acoo_temp_1225[count_nnz_E_inter] = acoo[i];
			rowind_temp_1225[count_nnz_E_inter] = rowind[i];
			colind_temp_1225[count_nnz_E_inter] = colind[i];
			count_nnz_E_inter++;
		}
	}
	num_nzero_Pmatrix_E_inter[myid] = count_nnz_E_inter;
	cout << "myid is " << myid << "   count_nnz_E_inter is " << count_nnz_E_inter << endl;

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure_method2(myid, acoo_temp_1225, num_nzero_Pmatrix_E_inter, rowind_temp_1225, colind_temp_1225);
	end_FETI = time(NULL);



	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);

	long long mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n) revise
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row() + 1;
		colm[i] = TriList[i].col() + 1;
	}
	TriList.clear();
	//cout << " my id is " << myid << endl;
	//cout << " mat_size is " << mat_size << endl;
	int* nnz_ZC_process = new int[num_domain];    //revise
	for (int j = 0; j < num_domain; j++) {
		nnz_ZC_process[j] = 0;
	}
	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	//cout << "after,  my id is " << myid << endl;
	long long face_mat_size = 0;  //revise
	for (int j = 0; j < num_domain; j++) {
		//cout << nnz_ZC_process[j] << endl;
		face_mat_size += nnz_ZC_process[j] + Edge_Boundary[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];  //revise
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}



	}
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;



	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;



	time_t start, end;

	start = time(NULL);

	omp_set_num_threads(64);

	if (myid == 0) {

		for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
			rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
			m_final_total[j].real(1.0);
		}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;


		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;



		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}

		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		//clock_t tt1, tt2;
		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
		iparm[1] = 3; iparm[24] = 2;

		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);

		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
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
		ofstream ofs1119("data_DG_RTC_Yagi_method2.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}



	//if (myid == 0) {
	//	cout << "time is " << time << endl;
	//	ofstream ofs1119("data_DG_RTC_SIW.csv", ios::app);
	//	ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
	//		<< time << "s" << ',' << matsize << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
	//		<< process_threads << endl;

	//}


	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
		//ofstream ofs732("Solution_face.csv", ios::app);
		//for (int i = 0; i < num_unk_boundary; i++)
		//{
		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
		//}

	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;

	//	free(fbr); delete[] unKnownUi;
	return  0;

}

int Solve_E_H_boundary_method3(Element* el, int myid) {

	Fbr_Setting(el, myid);
	//cout <<"matsize"<< matsize << endl;
#pragma omp  parallel for
	for (int i = 0; i < matsize; i++) {
		acoo[i].real = acoo1[i].real + acoo_ome[i].real * omega + acoo_ome_2[i].real * omega * omega;
		acoo[i].imag = acoo1[i].imag + acoo_ome[i].imag * omega + acoo_ome_2[i].imag * omega * omega;
		m[i].real = m1[i].real + m_ome[i].real * omega + m_ome_2[i].real * omega * omega;
		m[i].imag = m1[i].imag + m_ome[i].imag * omega + m_ome_2[i].imag * omega * omega;
	}

	int unknown_E_inter = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] - Edge_Boundary[myid];

	int unknown_E_boundary = Edge_Boundary[myid] - num_unknown_subdomain[myid][1];

	int nnz_temp_1225 = num_nzero_Pmatrix[myid];

	int* num_nzero_Pmatrix_E_inter = new int[num_domain];


	MKL_Complex16* acoo_temp_1225 = new MKL_Complex16[nnz_temp_1225];



	int* rowind_temp_1225 = new int[nnz_temp_1225];
	int* colind_temp_1225 = new int[nnz_temp_1225];

	int count_nnz_E_inter = 0;
	for (int i = 0; i < nnz_temp_1225; i++) {
		if (rowind[i] <= unknown_E_inter && colind[i] <= unknown_E_inter) {
			acoo_temp_1225[count_nnz_E_inter] = acoo[i];
			rowind_temp_1225[count_nnz_E_inter] = rowind[i];
			colind_temp_1225[count_nnz_E_inter] = colind[i];
			count_nnz_E_inter++;
		}
	}
	num_nzero_Pmatrix_E_inter[myid] = count_nnz_E_inter;
	cout << "myid is " << myid << "   count_nnz_E_inter is " << count_nnz_E_inter << endl;

	r_dm = new int[num_domain];
	r_dm[0] = 0;
	for (int i = 1; i < num_domain; ++i) {
		r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	}
	cout << "Edge_Boundary[" << myid << "] is " << Edge_Boundary[myid] << endl;

	time_t start_out, end_out;
	start_out = time(NULL);
	MPI_Barrier(MPI_COMM_WORLD);

	time_t start_FETI, end_FETI;
	start_FETI = time(NULL);
	FETI_like_procedure_method3(myid, acoo_temp_1225, num_nzero_Pmatrix_E_inter, rowind_temp_1225, colind_temp_1225);
	end_FETI = time(NULL);



	double time_FETI = (double)(end_FETI - start_FETI);
	cout << "myid is " << myid << "time is " << time_FETI << endl;
	MPI_Barrier(MPI_COMM_WORLD);
	end_out = time(NULL);
	double time_out = (double)(end_out - start_out);

	delete[]acoo_temp_1225; delete[]rowind_temp_1225; delete[]colind_temp_1225; delete[]num_nzero_Pmatrix_E_inter;




	int  mat_size = TriList.size();   //  nnz of Zi*Cij (j=1,2,3...,n) revise
	complex<double>* m_final = new complex<double>[mat_size];
	int* rowm = new int[mat_size], * colm = new int[mat_size];
#	pragma omp parallel for 
	for (int i = 0; i < mat_size; i++) {
		m_final[i].real(TriList[i].value().real());
		m_final[i].imag(TriList[i].value().imag());
		rowm[i] = TriList[i].row();
		colm[i] = TriList[i].col();
	}
	TriList.clear();
	cout << " my id is " << myid << endl;
	cout << " mat_size is " << mat_size << endl;
	int* nnz_ZC_process = new int[num_domain];    //revise
	for (int j = 0; j < num_domain; j++) {
		nnz_ZC_process[j] = 0;
	}
	MPI_Gather(&mat_size, 1, MPI_INT, nnz_ZC_process, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(nnz_ZC_process, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	//cout << "after,  my id is " << myid << endl;
	long long face_mat_size = 0;  //revise
	for (int j = 0; j < num_domain; j++) {
		//cout << nnz_ZC_process[j] << endl;
		face_mat_size += nnz_ZC_process[j];  //Edge_Boundary[j]  Number of non zero elements of identity matrix
	}
	if (myid == 0) {
		cout << "face_mat_size = " << face_mat_size << endl;
	}

	int* address_offset_mat = new int[num_domain];  //revise
	address_offset_mat[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset_mat[j] = address_offset_mat[j - 1] + nnz_ZC_process[j - 1];
	}

	int num_unk_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_unk_boundary += Edge_Boundary[j];
	}
	if (myid == 0) {
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
	}
	int* rowm_total = nullptr;
	int* colm_total = nullptr;
	complex<double>* m_final_total = nullptr;
	complex<double>* fbb_total = nullptr;
	if (myid == 0) {
		m_final_total = new complex<double>[face_mat_size];
		rowm_total = new int[face_mat_size];
		colm_total = new int[face_mat_size];
		fbb_total = new complex<double>[num_unk_boundary];
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total[j].imag(0.0);
			fbb_total[j].real(0.0);
		}
	}
	//cout << "test0000" << endl;
	complex<double>* fbb_temp = new complex<double>[Edge_Boundary[myid]];
	for (int j = 0; j < Edge_Boundary[myid]; j++) {
		fbb_temp[j].real(fbb[j].real);
		fbb_temp[j].imag(fbb[j].imag);
	}
	int* address_offset = new int[num_domain];

	address_offset[0] = 0;
	for (int j = 1; j < num_domain; j++) {
		address_offset[j] = address_offset[j - 1] + Edge_Boundary[j - 1];
	}

	//cout << "test1111" << endl;
	MPI_Gatherv(fbb_temp, Edge_Boundary[myid], MPI_DOUBLE_COMPLEX, fbb_total, Edge_Boundary, address_offset, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test2222" << endl;
	MPI_Gatherv(rowm, mat_size, MPI_INT, rowm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test3333" << endl;
	MPI_Gatherv(colm, mat_size, MPI_INT, colm_total, nnz_ZC_process, address_offset_mat, MPI_INT, 0, MPI_COMM_WORLD);
	//cout << "test4444" << endl;
	MPI_Gatherv(m_final, mat_size, MPI_DOUBLE_COMPLEX, m_final_total, nnz_ZC_process, address_offset_mat, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	//cout << "test5555" << endl;
	delete[]fbb_temp; delete[]rowm; delete[]colm; delete[]m_final;


	unKnownUi = new MKL_Complex16[num_unk_boundary];
	MKL_Complex16* m_final_mkl;


	//cout << "test1111" << endl;
	time_t start, end;

	start = time(NULL);

	omp_set_num_threads(64);

	if (myid == 0) {

		//for (int j = face_mat_size - num_unk_boundary; j < face_mat_size; j++) {
		//	rowm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
		//	colm_total[j] = j - (face_mat_size - num_unk_boundary) + 1;
		//	m_final_total[j].real(1.0);
		//}
		m_final_mkl = new MKL_Complex16[face_mat_size];
		//#pragma omp  parallel for
		for (int j = 0; j < face_mat_size; j++) {
			m_final_mkl[j].real = m_final_total[j].real();
			m_final_mkl[j].imag = m_final_total[j].imag();
		}

		delete[] m_final_total;
		MKL_Complex16* m_pardiso = new MKL_Complex16[face_mat_size];
		int* im = new int[num_unk_boundary + 1], * jm = new int[face_mat_size];
		int info; int nnz_temp = face_mat_size;
		int job[8] = { 2, 1, 1, 0, face_mat_size, 0, 0, 0 };
		//sparse_index_base_t indexM;
		mkl_zcsrcoo(job, &num_unk_boundary, m_pardiso, jm, im, &nnz_temp, m_final_mkl, rowm_total, colm_total, &info);

		delete[] m_final_mkl;
		delete[] rowm_total;
		delete[] colm_total;
		cout << "test7777" << endl;

		MKL_Complex16* fbb_total_mkl = new MKL_Complex16[num_unk_boundary];
		//ofstream ofs11141("fbb_total.csv", ios::app);
		for (int j = 0; j < num_unk_boundary; j++) {
			fbb_total_mkl[j].real = fbb_total[j].real();
			fbb_total_mkl[j].imag = fbb_total[j].imag();
			//ofs11141 << fbb_total[j].imag() << ',' << fbb_total[j].real() << endl;
		}
		delete[]  fbb_total;

		MKL_INT maxfct, mnum, mtype, phase1, phase2, msglvl, error;
		MKL_INT iparm[64];
		int pt[64];
		MKL_INT* perm2 = nullptr;
		maxfct = 1; mnum = 1; mtype = 13; phase1 = 13; phase2 = -1; msglvl = 0;



		for (int el = 0; el < 64; el++) {
			pt[el] = 0; iparm[el] = 0;
		}

		iparm[0] = 1;
		iparm[1] = 3;
		iparm[2] = 16;
		iparm[9] = 13;
		iparm[10] = 1;
		iparm[12] = 1;
		iparm[17] = -1;
		iparm[18] = -1;
		MKL_INT  nrhs = 1;
		//clock_t tt1, tt2;
		//tt1 = clock(); //iparm[59] = 2; msglvl = 1; 
		iparm[1] = 3; iparm[24] = 2;
		//cout << "test6666" << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase1, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);
		if (error != 0)cout << "error for boundary equation   " << error << endl;
		pardiso(pt, &maxfct, &mnum, &mtype, &phase2, &num_unk_boundary, m_pardiso, im, jm, perm2, &nrhs, iparm, &msglvl, fbb_total_mkl, unKnownUi, &error);

		//tt2 = clock();
		//cout << "task cost  " << (tt2 - tt1) / CLOCKS_PER_SEC << 's' << endl;
		delete[] fbb; delete[] m_pardiso; delete[] im; delete[] jm;
		delete[]fbb_total_mkl;
		// calulate full solution Um

	}
	//cout << "test2222" << endl;

	end = time(NULL);

	double time = (double)(end - start);

	int num_opp_element_boundary = 0;
	for (int j = 0; j < num_domain; j++) {
		num_opp_element_boundary += num_opp_element[j];
	}

	//cout << "test4444" << endl;

	if (myid == 0) {
		cout << "time is " << time << endl;

		ofstream ofs1119("data_DG_RTC_Yagi_method3.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
			<< time << "s" << ',' << num_nzero_Pmatrix[myid] << ',' << face_mat_size << ',' << Edge_Boundary[myid] - num_unknown_subdomain[0][1] << ','\
			<< process_threads << ',' << num_opp_element_boundary << ',' << matsize << endl;

	}



	//if (myid == 0) {
	//	cout << "time is " << time << endl;
	//	ofstream ofs1119("data_DG_RTC_SIW.csv", ios::app);
	//	ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_unk_boundary << ',' << time_out << "s" << ',' \
	//		<< time << "s" << ',' << matsize << ',' << face_mat_size << ',' << Edge_Boundary[myid] << ','\
	//		<< process_threads << endl;

	//}


	omp_set_num_threads(process_threads);

	MPI_Barrier(MPI_COMM_WORLD);

	complex<double>* unKnownUi_mpic = new complex<double>[num_unk_boundary];


	if (myid == 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi_mpic[j].real(unKnownUi[j].real);
			unKnownUi_mpic[j].imag(unKnownUi[j].imag);
		}
		cout << "num_unk_boundary = " << num_unk_boundary << endl;
		//ofstream ofs732("Solution_face.csv", ios::app);
		//for (int i = 0; i < num_unk_boundary; i++)
		//{
		//	ofs732 << unKnownUi[i].imag << ',' << unKnownUi[i].real << endl;
		//}

	}
	MPI_Bcast(unKnownUi_mpic, num_unk_boundary, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

	if (myid != 0) {
		for (int j = 0; j < num_unk_boundary; j++) {
			unKnownUi[j].real = unKnownUi_mpic[j].real();
			unKnownUi[j].imag = unKnownUi_mpic[j].imag();
		}
	}
	delete[]unKnownUi_mpic;
	int num_unknown_this_process = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1];
	Xsol = (complex<double>*)malloc(sizeof(complex<double>) * num_unknown_this_process);
	Solver_mpi(myid);
	free(fbr); delete[] unKnownUi;
	//cout << "test3333" << endl;
	//	free(fbr); delete[] unKnownUi;
	return  0;

}


int Eqn_Solver(Element* el, Element* oppo_element, int myid) {
	unordered_map<long long, complex<double>>mat_full;  //Frequency independent matrix terms
	unordered_map<long long, complex<double>>mat_full_ome;  //Matrix terms related to the power of frequency
	unordered_map<long long, complex<double>>mat_full_ome_2;  //Matrix terms related to the quadratic frequency
	unordered_map<long long, complex<double>>mat_full_ome_sqrt;  //Matrix terms related to the sqrt frequency
	unordered_map<long long, complex<double>> mat_full4;



	//prepare for code
	cout.precision(16);
	long long nnz = 0;  //Number of non zero elements in a sparse matrix 
	nnz_c = new size_t * [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_c[i] = new size_t[num_domain];

	int kk;

	int ndm, nn, mm, ii, jj, node_ii1, node_ii2, node_jj1, node_jj2;
	int nodeglb_ii1, nodeglb_ii2, nodeglb_jj1, nodeglb_jj2, op1, ofn1, edgeii_loc, edgejj_loc, edgeii_E, edgeii_H, edgejj_E, edgejj_H;
	double nx, ny, nz, Ni_dot_ncrosNj, ncrosNi_dot_ncrosNj, ncrosNix1, ncrosNix2, ncrosNiy1, ncrosNiy2, ncrosNiz1, ncrosNiz2;
	double ncrosNjx1, ncrosNjx2, ncrosNjy1, ncrosNjy2, ncrosNjz1, ncrosNjz2, signE_Ni, signE_Nj, signH_Ni, signH_Nj;
	double Li1[3], Li2[3], Lj1[3], Lj2[3], Ni11[3], Ni12[3], zb1[3], zb2[3], zb3[3], rs[3];
	double Nix1, Nix2, Niy1, Niy2, Niz1, Niz2, Njx1, Njx2, Njy1, Njy2, Njz1, Njz2;
	//complex <long double> nVec[3], Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	//complex<long double> EincX, EincY, EincZ, HincX, HincY, HincZ;

	double Ftemp[6][3];




	MKL_Complex16* Me_full = nullptr;
	int* iRow = nullptr;// store the rows of elements per matrix
	int* jCol = nullptr;// store the columns of elements per matrix

	long long temp = 0;

	//end prepare




	//  Matrix Fullfillment

	for (int i = 0; i < num_element_subdomain; i++) {
		for (ii = 0; ii < 6; ii++) {
			Ftemp[ii][0] = 2.0 * (el[i].c[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1] - el[i].d[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve); //The x component of the curl of the basis function
			Ftemp[ii][1] = 2.0 * (el[i].d[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1] - el[i].b[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);//The y component of the curl of the basis function
			Ftemp[ii][2] = 2.0 * (el[i].b[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1] - el[i].c[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);// The z component of the curl of the basis function
		}
		for (ii = 0; ii < 6; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			for (jj = 0; jj < 6; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					//temp = (long long)((edgeii_E - 1) * num_unknown + edgejj_E - 1);
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full[temp].real(mat_full[temp].real() + 1.0 / mur * signE_Ni * signE_Nj * el[i].length[ii] * el[i].length[jj] * el[i].Ve * \
						(Ftemp[ii][0] * Ftemp[jj][0] + Ftemp[ii][1] * Ftemp[jj][1] + Ftemp[ii][2] * Ftemp[jj][2]));
					mat_full_ome[temp] += 0.0;
					mat_full_ome_2[temp] += 0.0;
					mat_full4[temp] += 0.0;
					mat_full_ome_sqrt[temp] += 0.0;
				}
			}
		}
	}


	for (int i = 0; i < num_element_subdomain; i++) {
		if (el[i].domain - 1 != myid) {
			cout << "error" << endl;
		}
		for (ii = 0; ii < 12; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			node_ii1 = edge_node_local[ii][0]; node_ii2 = edge_node_local[ii][1];
			nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
			Li1[0] = sign_judge(ii, 5) * el[i].b[node_ii1 - 1] * el[i].length[ii];
			Li1[1] = sign_judge(ii, 5) * el[i].c[node_ii1 - 1] * el[i].length[ii];
			Li1[2] = sign_judge(ii, 5) * el[i].d[node_ii1 - 1] * el[i].length[ii];
			Li2[0] = el[i].b[node_ii2 - 1] * el[i].length[ii];
			Li2[1] = el[i].c[node_ii2 - 1] * el[i].length[ii];
			Li2[2] = el[i].d[node_ii2 - 1] * el[i].length[ii];
			for (jj = 0; jj < 12; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				node_jj1 = edge_node_local[jj][0]; node_jj2 = edge_node_local[jj][1];
				nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
				Lj1[0] = sign_judge(jj, 5) * el[i].b[node_jj1 - 1] * el[i].length[jj];
				Lj1[1] = sign_judge(jj, 5) * el[i].c[node_jj1 - 1] * el[i].length[jj];
				Lj1[2] = sign_judge(jj, 5) * el[i].d[node_jj1 - 1] * el[i].length[jj];
				Lj2[0] = el[i].b[node_jj2 - 1] * el[i].length[jj];
				Lj2[1] = el[i].c[node_jj2 - 1] * el[i].length[jj];
				Lj2[2] = el[i].d[node_jj2 - 1] * el[i].length[jj];

				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full_ome_2[temp].real(mat_full_ome_2[temp].real() - el[i].epsil * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);//
					//mat_full_ome[temp] += 0.0;
					mat_full_ome[temp].imag(mat_full_ome[temp].imag() + el[i].sigma * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);   //sigma!=0
					mat_full[temp] += 0.0;
					mat_full4[temp] += 0.0;
					mat_full_ome_sqrt[temp] += 0.0;

				}
			}
		}
	}



	//Absorbing boundary condition term

	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && (ofn1 == -7)) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							//cout << "error" << endl;
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
							mat_full4[temp] += 0.0;
							mat_full_ome_sqrt[temp] += 0.0;
						}
					}
				}
			}
		}
	}


	//IBC
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if (ofn1 == -55) {
				cout << "error" << endl;
				double sigma_face_temp = el[i].face[nn].sigma_face;
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);

					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;

						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / (Z0 * Wmicrostrip / Hsubstrate));
							//mat_full4[temp].imag(mat_full4[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
							//mat_full_ome_sqrt[temp].imag(mat_full_ome_sqrt[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj * sqrt(sigma_IBC / 2 / mur));
							//mat_full_ome_sqrt[temp].real(mat_full_ome_sqrt[temp].real() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj * sqrt(sigma_IBC / 2 / mur));
							mat_full_ome_sqrt[temp].imag(mat_full_ome_sqrt[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj * sqrt(sigma_face_temp / 2 / mur));
							mat_full_ome_sqrt[temp].real(mat_full_ome_sqrt[temp].real() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj * sqrt(sigma_face_temp / 2 / mur));
							mat_full4[temp] += 0.0;

							mat_full_ome[temp] += 0.0;
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
							//mat_full_ome_sqrt[temp] += 0.0;

						}
					}
				}
			}
		}
	}


	//Waveport1
	for (int i = 0; i < num_element_subdomain; ++i) {
		double zt = (el[i].node[0].zb[2] + el[i].node[1].zb[2] + el[i].node[2].zb[2] + el[i].node[3].zb[2]) / 4.0;
		//if (el[i].Material != wave_material) {
		if (zt< wave_zzz) {
			continue;
		}



		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && (-215<=ofn1&& ofn1 <= -200)) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							//cout << "error" << endl;
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0*signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
							mat_full4[temp] += 0.0;
							mat_full_ome_sqrt[temp] += 0.0;
						}
					}
				}
			}
		}
	}

	//Waveport2
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ( ofn1 == -66) {
				cout << "error" << endl;
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);

					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;

						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / (Z0 * Wmicrostrip / Hsubstrate));
							//mat_full4[temp].imag(mat_full4[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
							mat_full4[temp].imag(mat_full4[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
							mat_full_ome[temp] += 0.0;
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
							mat_full_ome_sqrt[temp] += 0.0;

						}
					}
				}
			}
		}
	}

	int size_temp = mat_full.size();
	//cout << "myid is " << myid << endl;

	//cout << "size_temp is " << size_temp << endl;

	// Fehii Feeii Fehij  Feeij
	// Fehii Feeii Fehij  Feeij

	//Total number of  elements in Stiffness matrix, mass matrix, and ABC matrix

	//contribution from numerical flux, numerical flux has two parts : outgoing fluxand incoming flux
	int count_oppo_element = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		if (mesh_at_RTC[i] != 0) {    //on the RTC boundary
			for (nn = 0; nn < 4; nn++) {
				op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
				if (op1 != 0) {
					//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
					if (el[i].face[nn].whether_boundary) {
						//if (el[op1 - 1].domain != el[i].domain) {
						if (oppo_element[count_oppo_element].Global_num == op1) {
							op1 = count_oppo_element + 1;
							count_oppo_element++;
							//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
						}
						//else {
						//	for (int temp = 0; temp < num_opp_element[myid]; temp++) {
						//		if (oppo_element[temp].Global_num == op1) {
						//			op1 = temp+1;
						//			cout << "wrong" << endl;
						//			count_oppo_element++;
						//			break;
						//		}
						//	}
						//}
						nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
							signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];//
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//outgoing flux
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in outgoing flux：\hat{ n } \times \mathbf{ H } ^ i
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//!second term in outgoing flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ i }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									//This term should be integraed
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//incoming Flux
								//if (count_oppo_element == 1) {
								//	cout << FactorH * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj << endl;
								//	//cout << "ttttttttttttttttttttt" << endl;
								//}
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];

								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];

								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;

								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//first term in incoming flux：\hat{ n } \times \mathbf{ H } ^ j
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//second term in incoming flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ j }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
							}
						}
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_H = el[i].Hedge_GBNO[edgeii_loc - 1];
							signH_Ni = 1.0 * el[i].Hedge_GBNO[edgeii_loc + 12 - 1];
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni;
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//RTC Terms in current subdomain
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in RTC：\hat{ n } \times \mathbf{ E }^ i
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
										/*  mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
										  mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;

									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//!second term in RTC : \hat{ n } \times  \mathbf{ H }^ i} \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;

									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//RTC terms in neighboring subdomain
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);

								//edgejj_E = el[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								//signE_Nj = 1.0 * el[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								//node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								//nodeglb_jj1 = el[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = el[op1 - 1].node[node_jj2 - 1].ver;
								////x, y, z components of n cross Nj
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//third term in RTC：\hat{n} \times \mathbf{E} ^ j
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;

									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//fourth term in RTC : \hat{n} \times  \mathbf{H}^ i } \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									/*mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									/*	mat_full_ome[temp] += 0;
										mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									mat_full_ome_2[temp] += 0;

									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
							}
						}
					}
				}
			}
		}
	}
	//cout << "myid" << " = " << myid << endl;
	//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
	//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;


	//eq H


	/*        
	for (int i = 0; i < num_element_subdomain; i++) {
		for (ii = 0; ii < 6; ii++) {
			Ftemp[ii][0] = 2.0 * (el[i].c[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1] - el[i].d[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve); //The x component of the curl of the basis function
			Ftemp[ii][1] = 2.0 * (el[i].d[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1] - el[i].b[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);//The y component of the curl of the basis function
			Ftemp[ii][2] = 2.0 * (el[i].b[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1] - el[i].c[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);// The z component of the curl of the basis function
		}
		for (ii = 0; ii < 6; ii++) {
			edgeii_H = el[i].Hedge_GBNO[ii];
			signH_Ni = 1.0 * el[i].Hedge_GBNO[ii + 12];

			for (jj = 0; jj < 6; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_H = el[i].Hedge_GBNO[jj];
				signH_Nj = 1.0 * el[i].Hedge_GBNO[jj + 12];
				if ((edgeii_H != 0) && (edgejj_H != 0)) {
					//temp = (long long)((edgeii_E - 1) * num_unknown + edgejj_E - 1);
					temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full[temp].real(mat_full[temp].real() + 1.0 / mur * signH_Ni * signH_Nj * el[i].length[ii] * el[i].length[jj] * el[i].Ve * \
						(Ftemp[ii][0] * Ftemp[jj][0] + Ftemp[ii][1] * Ftemp[jj][1] + Ftemp[ii][2] * Ftemp[jj][2]));
					mat_full_ome[temp] += 0.0;
					mat_full_ome_2[temp] += 0.0;
					mat_full4[temp] += 0.0;
					mat_full_ome_sqrt[temp] += 0.0;
				}
			}
		}
	}

	for (int i = 0; i < num_element_subdomain; i++) {
		if (el[i].domain - 1 != myid) {
			cout << "error" << endl;
		}
		for (ii = 0; ii < 12; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];

			edgeii_H = el[i].Hedge_GBNO[ii];
			signH_Ni = 1.0 * el[i].Hedge_GBNO[ii + 12];
			node_ii1 = edge_node_local[ii][0]; node_ii2 = edge_node_local[ii][1];
			nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
			Li1[0] = sign_judge(ii, 5) * el[i].b[node_ii1 - 1] * el[i].length[ii];
			Li1[1] = sign_judge(ii, 5) * el[i].c[node_ii1 - 1] * el[i].length[ii];
			Li1[2] = sign_judge(ii, 5) * el[i].d[node_ii1 - 1] * el[i].length[ii];
			Li2[0] = el[i].b[node_ii2 - 1] * el[i].length[ii];
			Li2[1] = el[i].c[node_ii2 - 1] * el[i].length[ii];
			Li2[2] = el[i].d[node_ii2 - 1] * el[i].length[ii];
			for (jj = 0; jj < 12; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				edgejj_H = el[i].Hedge_GBNO[jj];
				signH_Nj = 1.0 * el[i].Hedge_GBNO[jj + 12];
				node_jj1 = edge_node_local[jj][0]; node_jj2 = edge_node_local[jj][1];
				nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
				Lj1[0] = sign_judge(jj, 5) * el[i].b[node_jj1 - 1] * el[i].length[jj];
				Lj1[1] = sign_judge(jj, 5) * el[i].c[node_jj1 - 1] * el[i].length[jj];
				Lj1[2] = sign_judge(jj, 5) * el[i].d[node_jj1 - 1] * el[i].length[jj];
				Lj2[0] = el[i].b[node_jj2 - 1] * el[i].length[jj];
				Lj2[1] = el[i].c[node_jj2 - 1] * el[i].length[jj];
				Lj2[2] = el[i].d[node_jj2 - 1] * el[i].length[jj];

				if ((edgeii_H != 0) && (edgejj_H != 0)) {
					temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full_ome_2[temp].real(mat_full_ome_2[temp].real() - el[i].epsil * signH_Ni * signH_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);//
					//mat_full_ome[temp] += 0.0;
					mat_full_ome[temp] += 0.0;
					mat_full[temp] += 0.0;
					mat_full4[temp] += 0.0;
					mat_full_ome_sqrt[temp] += 0.0;

				}
			}
		}
	}

	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && (ofn1 == -7)) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];

					edgeii_H = el[i].Hedge_GBNO[ii];
					signH_Ni = 1.0 * el[i].Hedge_GBNO[ii + 12];

					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];

						edgejj_H = el[i].Hedge_GBNO[jj];
						signH_Nj = 1.0 * el[i].Hedge_GBNO[jj + 12];

						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_H != 0) && (edgejj_H != 0)) {
							//cout << "error" << endl;
							temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
							mat_full4[temp] += 0.0;
							mat_full_ome_sqrt[temp] += 0.0;
						}
					}
				}
			}
		}
	}

	count_oppo_element = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		if (mesh_at_RTC[i] != 0) {    //on the RTC boundary
			for (nn = 0; nn < 4; nn++) {
				op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
				if (op1 != 0) {
					//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
					if (el[i].face[nn].whether_boundary) {
						//if (el[op1 - 1].domain != el[i].domain) {
						if (oppo_element[count_oppo_element].Global_num == op1) {
							op1 = count_oppo_element + 1;
							count_oppo_element++;
							//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
						}
						//else {
						//	for (int temp = 0; temp < num_opp_element[myid]; temp++) {
						//		if (oppo_element[temp].Global_num == op1) {
						//			op1 = temp+1;
						//			cout << "wrong" << endl;
						//			count_oppo_element++;
						//			break;
						//		}
						//	}
						//}
						nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
							signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];//

							edgeii_H = el[i].Hedge_GBNO[edgeii_loc - 1];

							signH_Ni = 1.0 * el[i].Hedge_GBNO[edgeii_loc + 12 - 1];

							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//outgoing flux
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in outgoing flux：\hat{ n } \times \mathbf{ H } ^ i
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//!second term in outgoing flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ i }
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									//This term should be integraed
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() +FactorH * el[i].eta * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//incoming Flux
								//if (count_oppo_element == 1) {
								//	cout << FactorH * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj << endl;
								//	//cout << "ttttttttttttttttttttt" << endl;
								//}
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];

								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];

								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;

								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//first term in incoming flux：\hat{ n } \times \mathbf{ H } ^ j
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0 * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
								//second term in incoming flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ j }
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * el[i].eta * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
									mat_full4[temp] += 0.0;
									mat_full_ome_sqrt[temp] += 0.0;
								}
							}
						}
					}
				}
			}
		}
	}
	*/
	







	nnz = 0;
	matsize = mat_full.size();
	Me_full = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome_2 = new MKL_Complex16[matsize];
	iRow = new int[matsize]; jCol = new int[matsize];
#   pragma omp parallel for 
	for (int i = 0; i < matsize; ++i) {
		Me_full[i].real = (0.0); Me_full[i].imag = (0.0);
		Me_full_ome[i].real = (0.0); Me_full_ome[i].imag = (0.0);
		Me_full_ome_2[i].real = (0.0); Me_full_ome_2[i].imag = (0.0);
		iRow[i] = 0; jCol[i] = 0;
	}

	for (auto& m : mat_full) {
		Me_full[nnz].real = m.second.real();
		Me_full[nnz].imag = m.second.imag();
		iRow[nnz] = m.first / num_unknown + 1;
		jCol[nnz] = m.first % num_unknown + 1;
		++nnz;
	}
	mat_full.clear();
	//cout << "mat_full =  " << nnz << endl;


	//mat_full_ome
	nnz = 0;
	matsize = mat_full_ome.size();
	for (auto& m : mat_full_ome) {
		Me_full_ome[nnz].real = (m.second).real();//?
		Me_full_ome[nnz].imag = (m.second).imag();//?
		++nnz;
	}
	mat_full_ome.clear();


	nnz = 0;
	matsize = mat_full_ome_2.size();
	for (auto& m : mat_full_ome_2) {
		Me_full_ome_2[nnz].real = m.second.real();
		Me_full_ome_2[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full_ome_2.clear();
	//cout << "mat_full_ome_2 = " << nnz << endl;

	nnz = 0;
	MKL_Complex16* Me_full4 = new MKL_Complex16[matsize];
	matsize = mat_full4.size();
	for (auto& m : mat_full4) {
		Me_full4[nnz].real = m.second.real();
		Me_full4[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full4.clear();

	nnz = 0;
	MKL_Complex16* Me_full_ome_sqrt = new MKL_Complex16[matsize];
	matsize = mat_full_ome_sqrt.size();
	for (auto& m : mat_full_ome_sqrt) {
		Me_full_ome_sqrt[nnz].real = m.second.real();
		Me_full_ome_sqrt[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full_ome_sqrt.clear();


	int nn2, mm2;
	acoo1 = new MKL_Complex16[nnz];
	acoo = new MKL_Complex16[nnz];
	acoo_ome = new MKL_Complex16[nnz];
	acoo_ome_2 = new MKL_Complex16[nnz];
	acoo_wave = new MKL_Complex16[nnz];
	acoo_ome_sqrt = new MKL_Complex16[nnz];


	rowind = new int[nnz];
	colind = new int[nnz];

	num_nzero_Pmatrix = new int[num_domain];
	//for (int i = 0; i < num_domain; ++i) {
	//	num_nzero_Pmatrix[i] = 0;
	//}

	int nnz_tot = 0;

	m = new MKL_Complex16[nnz];
	m_ome = new MKL_Complex16[nnz];
	m_ome_2 = new MKL_Complex16[nnz];
	m_wave = new MKL_Complex16[nnz];
	m_ome_sqrt = new MKL_Complex16[nnz];
	m1 = new MKL_Complex16[nnz];

	mrow = new int[nnz];
	mcol = new int[nnz];


	size_t cc = 0;
	int nnz_dm;
	nn = 0;
	for (ndm = myid; ndm < myid + 1; ++ndm) {
		mm = nn + num_unknown_subdomain[ndm][0] + num_unknown_subdomain[ndm][1];
		nn2 = 0;
		for (size_t ndm2 = 0; ndm2 != num_domain; ++ndm2) {
			nnz_dm = 0;
			mm2 = nn2 + num_unknown_subdomain[ndm2][0] + num_unknown_subdomain[ndm2][1];
#      pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow[j] > nn && iRow[j] <= mm && jCol[j] > nn2 && jCol[j] <= mm2) {//?

#	   pragma omp critical
					{
						++nnz_dm;
						m1[cc] = Me_full[j];
						m_ome[cc] = Me_full_ome[j];
						m_ome_2[cc] = Me_full_ome_2[j];
						m_wave[cc] = Me_full4[j];
						m_ome_sqrt[cc] = Me_full_ome_sqrt[j];

						mrow[cc] = iRow[j] - nn;
						mcol[cc] = jCol[j] - nn2;
						++cc;
					}
					if (ndm == ndm2) {
#	 pragma omp critical
						{

							acoo1[nnz_tot] = Me_full[j];
							acoo_ome[nnz_tot] = Me_full_ome[j];
							acoo_ome_2[nnz_tot] = Me_full_ome_2[j];
							acoo_wave[nnz_tot] = Me_full4[j];
							acoo_ome_sqrt[nnz_tot] = Me_full_ome_sqrt[j];

							rowind[nnz_tot] = iRow[j] - nn;
							colind[nnz_tot] = jCol[j] - nn2;
							nnz_tot += 1;

						}
					}
				}
			}
			nn2 = mm2;
			nnz_c[ndm][ndm2] = nnz_dm;//The number of data to fill in each matrix
			if (ndm == ndm2) num_nzero_Pmatrix[ndm] = nnz_dm;
		}
		nn = mm;
	}
	delete[] Me_full; delete[] Me_full_ome; delete[] Me_full_ome_2; delete[] Me_full4; delete[] Me_full_ome_sqrt;
	delete[] iRow;
	delete[] jCol;
	cout << "cc = " << cc << endl;

	int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
	p_pardiso = new MKL_Complex16[nnz_tot];
	ip = new int[number_unknown_subdomain + 1];
	jp = new int[nnz_tot];

	cout << "myid is " << myid << endl;
	cout << "nnz_tot = " << nnz_tot << endl;
	//cout << "num_nzero_Pmatrix[myid] = " << num_nzero_Pmatrix[myid] << endl;

	/*int nzero1 = 0; int row1 = 0;*/

	//nnz_tot = 0;
	//nn2 = 0;

	row_dm = new int[num_domain];
	nzero_dm = new int[num_domain];
	nn2_dm = new int[num_domain];
	nn3_dm = new int[num_domain];

	n3_dm = new int[num_domain];

	//row_dm[0] = 0; nzero_dm[0] = 0; nn2_dm[0] = 0; nn3_dm[0] = 0; r_dm[0] = 0; n3_dm[0] = 0;
	//int unknown_dm;
	//for (int i = 1; i < num_domain; ++i) {
	//	unknown_dm = num_unknown_subdomain[i - 1][0] + num_unknown_subdomain[i - 1][1];
	//	nnz_dm = num_nzero_Pmatrix[i - 1];

	//	r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	//	n3_dm[i] = n3_dm[i - 1] + Edge_Boundary[i - 1];

	//	row_dm[i] = row_dm[i - 1] + unknown_dm + 1;//第i个domain的row坐标
	//	nzero_dm[i] = nzero_dm[i - 1] + nnz_dm;//第i个domain的坐标  统计的是P矩阵的nnz数目
	//	nn2_dm[i] = nn2_dm[i - 1] + unknown_dm;
	//	nn3_dm[i] = nn3_dm[i - 1];

	//	for (int j = 0; j < num_domain; j++) {
	//		nn3_dm[i] += nnz_c[i - 1][j];
	//	}
	//}

	return 0;
}

int Eqn_Solver_before(Element* el, Element* oppo_element,int myid) {
	unordered_map<long long, complex<double>>mat_full;  //Frequency independent matrix terms
	unordered_map<long long, complex<double>>mat_full_ome;  //Matrix terms related to the power of frequency
	unordered_map<long long, complex<double>>mat_full_ome_2;  //Matrix terms related to the quadratic frequency
	//prepare for code
	cout.precision(16);
	long long nnz = 0;  //Number of non zero elements in a sparse matrix 
	nnz_c = new size_t *[num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_c[i] = new size_t[num_domain];

	int kk;

	int ndm, nn, mm, ii, jj, node_ii1, node_ii2, node_jj1, node_jj2     ;
	int nodeglb_ii1, nodeglb_ii2, nodeglb_jj1, nodeglb_jj2, op1, ofn1, edgeii_loc, edgejj_loc, edgeii_E, edgeii_H, edgejj_E, edgejj_H;
	double nx, ny, nz, Ni_dot_ncrosNj, ncrosNi_dot_ncrosNj, ncrosNix1, ncrosNix2, ncrosNiy1, ncrosNiy2, ncrosNiz1, ncrosNiz2;
	double ncrosNjx1, ncrosNjx2, ncrosNjy1, ncrosNjy2, ncrosNjz1, ncrosNjz2, signE_Ni, signE_Nj, signH_Ni, signH_Nj;
	double Li1[3], Li2[3], Lj1[3], Lj2[3], Ni11[3], Ni12[3], zb1[3], zb2[3], zb3[3], rs[3];
	double Nix1, Nix2, Niy1, Niy2, Niz1, Niz2, Njx1, Njx2, Njy1, Njy2, Njz1, Njz2;
	//complex <long double> nVec[3], Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	//complex<long double> EincX, EincY, EincZ, HincX, HincY, HincZ;

	double Ftemp[6][3];




	MKL_Complex16* Me_full = nullptr;
	int* iRow = nullptr;// store the rows of elements per matrix
	int* jCol = nullptr;// store the columns of elements per matrix

	long long temp = 0;

	//end prepare




	//  Matrix Fullfillment

	for (int i = 0; i < num_element_subdomain; i++) {
		for (ii = 0; ii < 6; ii++) {
			Ftemp[ii][0] = 2.0 * (el[i].c[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1] - el[i].d[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve); //The x component of the curl of the basis function
			Ftemp[ii][1] = 2.0 * (el[i].d[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1] - el[i].b[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);//The y component of the curl of the basis function
			Ftemp[ii][2] = 2.0 * (el[i].b[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1] - el[i].c[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);// The z component of the curl of the basis function
		}
		for (ii = 0; ii < 6; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			for (jj = 0; jj < 6; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					//temp = (long long)((edgeii_E - 1) * num_unknown + edgejj_E - 1);
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full[temp].real(mat_full[temp].real() + 1.0 / mur * signE_Ni * signE_Nj * el[i].length[ii] * el[i].length[jj] * el[i].Ve *\
						(Ftemp[ii][0] * Ftemp[jj][0] + Ftemp[ii][1] * Ftemp[jj][1] + Ftemp[ii][2] * Ftemp[jj][2]));
					mat_full_ome[temp] += 0.0;
					mat_full_ome_2[temp] += 0.0;
				}
			}
		}
	}


	for (int i = 0; i < num_element_subdomain; i++) {
		if (el[i].domain - 1 != myid) {
			cout << "error" << endl;
		}
		for (ii = 0; ii < 12; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			node_ii1 = edge_node_local[ii][0]; node_ii2 = edge_node_local[ii][1];
			nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
			Li1[0] = sign_judge(ii, 5) * el[i].b[node_ii1 - 1] * el[i].length[ii];
			Li1[1] = sign_judge(ii, 5) * el[i].c[node_ii1 - 1] * el[i].length[ii];
			Li1[2] = sign_judge(ii, 5) * el[i].d[node_ii1 - 1] * el[i].length[ii];
			Li2[0] = el[i].b[node_ii2 - 1] * el[i].length[ii];
			Li2[1] = el[i].c[node_ii2 - 1] * el[i].length[ii];
			Li2[2] = el[i].d[node_ii2 - 1] * el[i].length[ii];
			for (jj = 0; jj < 12; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				node_jj1 = edge_node_local[jj][0]; node_jj2 = edge_node_local[jj][1];
				nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
				Lj1[0] = sign_judge(jj, 5) * el[i].b[node_jj1 - 1] * el[i].length[jj];
				Lj1[1] = sign_judge(jj, 5) * el[i].c[node_jj1 - 1] * el[i].length[jj];
				Lj1[2] = sign_judge(jj, 5) * el[i].d[node_jj1 - 1] * el[i].length[jj];
				Lj2[0] = el[i].b[node_jj2 - 1] * el[i].length[jj];
				Lj2[1] = el[i].c[node_jj2 - 1] * el[i].length[jj];
				Lj2[2] = el[i].d[node_jj2 - 1] * el[i].length[jj];

				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full_ome_2[temp].real(mat_full_ome_2[temp].real() - el[i].epsil * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);//
					//mat_full_ome[temp] += 0.0;
					mat_full_ome[temp].imag(mat_full_ome[temp].imag() + el[i].sigma * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);   //sigma!=0


					mat_full[temp] += 0.0;
				}
			}
		}
	}



	//Absorbing boundary condition term
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && ((ofn1 == -7) || (ofn1 == -8))) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);


					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
						}
					}
				}
			}
		}
	}

	//LP
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if (ofn1 == -6|| ofn1 == -66) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);

					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;

						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / (Z0 * Wmicrostrip / Hsubstrate));
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;

						}
					}
				}
			}
		}
	}

	int size_temp = mat_full.size();
	//cout << "myid is " << myid << endl;

	//cout << "size_temp is " << size_temp << endl;

	// Fehii Feeii Fehij  Feeij
	// Fehii Feeii Fehij  Feeij

	//Total number of  elements in Stiffness matrix, mass matrix, and ABC matrix

	//contribution from numerical flux, numerical flux has two parts : outgoing fluxand incoming flux
	int count_oppo_element = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		if (mesh_at_RTC[i] != 0) {    //on the RTC boundary
			for (nn = 0; nn < 4; nn++) {
				op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
				if (op1 != 0) {
					//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
					if(el[i].face[nn].whether_boundary){
					//if (el[op1 - 1].domain != el[i].domain) {
						if (oppo_element[count_oppo_element].Global_num == op1) {
							op1 = count_oppo_element+1;
							count_oppo_element++;
							//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
						}
						//else {
						//	for (int temp = 0; temp < num_opp_element[myid]; temp++) {
						//		if (oppo_element[temp].Global_num == op1) {
						//			op1 = temp+1;
						//			cout << "wrong" << endl;
						//			count_oppo_element++;
						//			break;
						//		}
						//	}
						//}
						nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
							signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];//
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//outgoing flux
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in outgoing flux：\hat{ n } \times \mathbf{ H } ^ i
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta)  * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//!second term in outgoing flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ i }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									//This term should be integraed
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//incoming Flux
								//if (count_oppo_element == 1) {
								//	cout << FactorH * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj << endl;
								//	//cout << "ttttttttttttttttttttt" << endl;
								//}
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];

								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];

								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;

								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//first term in incoming flux：\hat{ n } \times \mathbf{ H } ^ j
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta)  * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//second term in incoming flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ j }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) *  signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
							}
						}
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_H = el[i].Hedge_GBNO[edgeii_loc - 1];
							signH_Ni = 1.0 * el[i].Hedge_GBNO[edgeii_loc + 12 - 1];
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni;
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//RTC Terms in current subdomain
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in RTC：\hat{ n } \times \mathbf{ E }^ i
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
								/*	mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
								  /*  mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//!second term in RTC : \hat{ n } \times  \mathbf{ H }^ i} \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
								/*	mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
								/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//RTC terms in neighboring subdomain
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);

								//edgejj_E = el[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								//signE_Nj = 1.0 * el[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								//node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								//nodeglb_jj1 = el[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = el[op1 - 1].node[node_jj2 - 1].ver;
								////x, y, z components of n cross Nj
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//third term in RTC：\hat{n} \times \mathbf{E} ^ j
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
								/*	mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
								/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//fourth term in RTC : \hat{n} \times  \mathbf{H}^ i } \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									/*mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
								/*	mat_full_ome[temp] += 0;
									mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									mat_full_ome_2[temp] += 0;
								}
							}
						}
					}
				}
			}
		}
	}
	//cout << "myid" << " = " << myid << endl;
	//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
	//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;

	nnz = 0;
	matsize = mat_full.size();
	Me_full = new MKL_Complex16[matsize];
	MKL_Complex16 *Me_full_ome = new MKL_Complex16[matsize];
	MKL_Complex16 *Me_full_ome_2 = new MKL_Complex16[matsize];
	iRow = new int[matsize]; jCol = new int[matsize];
#   pragma omp parallel for 
	for (int i = 0; i < matsize; ++i) {
		Me_full[i].real = (0.0); Me_full[i].imag = (0.0);
		Me_full_ome[i].real = (0.0); Me_full_ome[i].imag = (0.0);
		Me_full_ome_2[i].real = (0.0); Me_full_ome_2[i].imag = (0.0);
		iRow[i] = 0; jCol[i] = 0;
	}

	for (auto& m : mat_full) {
		Me_full[nnz].real = m.second.real();
		Me_full[nnz].imag = m.second.imag();
		iRow[nnz] = m.first / num_unknown + 1;
		jCol[nnz] = m.first % num_unknown + 1;
		++nnz;
	}
	mat_full.clear();
	//cout << "mat_full =  " << nnz << endl;


	//mat_full_ome
	nnz = 0;
	matsize = mat_full_ome.size();
	for (auto& m : mat_full_ome) {
		Me_full_ome[nnz].real = (m.second).real();//?
		Me_full_ome[nnz].imag = (m.second).imag();//?
		++nnz;
	}
	mat_full_ome.clear();


	nnz = 0;
	matsize = mat_full_ome_2.size();
	for (auto& m : mat_full_ome_2) {
		Me_full_ome_2[nnz].real = m.second.real();
		Me_full_ome_2[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full_ome_2.clear();
	//cout << "mat_full_ome_2 = " << nnz << endl;

	int nn2, mm2;
	acoo1 = new MKL_Complex16[nnz];
	acoo = new MKL_Complex16[nnz];
	acoo_ome = new MKL_Complex16[nnz];
	acoo_ome_2 = new MKL_Complex16[nnz];

	rowind = new int[nnz];
	colind = new int[nnz];

	num_nzero_Pmatrix = new int[num_domain];
	//for (int i = 0; i < num_domain; ++i) {
	//	num_nzero_Pmatrix[i] = 0;
	//}
	
	int nnz_tot = 0;

	m = new MKL_Complex16[nnz];
	m_ome = new MKL_Complex16[nnz];
	m_ome_2 = new MKL_Complex16[nnz];
	m1 = new MKL_Complex16[nnz];

	mrow = new int[nnz];
	mcol = new int[nnz];


	size_t cc = 0;
	int nnz_dm;
	nn = 0;
	for (ndm = myid; ndm < myid+1; ++ndm) {
		mm = nn + num_unknown_subdomain[ndm][0] + num_unknown_subdomain[ndm][1];
		nn2 = 0;
		for (size_t ndm2 = 0; ndm2 != num_domain; ++ndm2) {
			nnz_dm = 0; 
			mm2 = nn2 + num_unknown_subdomain[ndm2][0] + num_unknown_subdomain[ndm2][1];
#      pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow[j] > nn && iRow[j] <= mm && jCol[j] > nn2 && jCol[j] <= mm2) {//?

#	   pragma omp critical
					{
						++nnz_dm;
						m1[cc] = Me_full[j];
						m_ome[cc] = Me_full_ome[j];
						m_ome_2[cc] = Me_full_ome_2[j];

						mrow[cc] = iRow[j] - nn;
						mcol[cc] = jCol[j] - nn2;
						++cc;
					}
					if (ndm == ndm2) {
#	 pragma omp critical
						{

							acoo1[nnz_tot] = Me_full[j];
							acoo_ome[nnz_tot] = Me_full_ome[j];
							acoo_ome_2[nnz_tot] = Me_full_ome_2[j];
							rowind[nnz_tot] = iRow[j] - nn;
							colind[nnz_tot] = jCol[j] - nn2;
							nnz_tot += 1;

						}
					}
				}
			}
			nn2 = mm2;
			nnz_c[ndm][ndm2] = nnz_dm;//The number of data to fill in each matrix
			if (ndm == ndm2) num_nzero_Pmatrix[ndm] = nnz_dm;
		}
		nn = mm;
	}
	delete[] Me_full; delete[] Me_full_ome; delete[] Me_full_ome_2;
	delete[] iRow;
	delete[] jCol;
	cout << "cc = " << cc << endl;

	int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
	p_pardiso = new MKL_Complex16[nnz_tot];
	ip = new int[number_unknown_subdomain + 1];
	jp = new int[nnz_tot];

	cout << "myid is " << myid << endl;
	cout << "nnz_tot = " << nnz_tot << endl;
	//cout << "num_nzero_Pmatrix[myid] = " << num_nzero_Pmatrix[myid] << endl;

	/*int nzero1 = 0; int row1 = 0;*/

	//nnz_tot = 0;
	//nn2 = 0;

	row_dm = new int[num_domain];
	nzero_dm = new int[num_domain];
	nn2_dm = new int[num_domain];
	nn3_dm = new int[num_domain];
	
	n3_dm = new int[num_domain];

	//row_dm[0] = 0; nzero_dm[0] = 0; nn2_dm[0] = 0; nn3_dm[0] = 0; r_dm[0] = 0; n3_dm[0] = 0;
	//int unknown_dm;
	//for (int i = 1; i < num_domain; ++i) {
	//	unknown_dm = num_unknown_subdomain[i - 1][0] + num_unknown_subdomain[i - 1][1];
	//	nnz_dm = num_nzero_Pmatrix[i - 1];

	//	r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	//	n3_dm[i] = n3_dm[i - 1] + Edge_Boundary[i - 1];

	//	row_dm[i] = row_dm[i - 1] + unknown_dm + 1;//第i个domain的row坐标
	//	nzero_dm[i] = nzero_dm[i - 1] + nnz_dm;//第i个domain的坐标  统计的是P矩阵的nnz数目
	//	nn2_dm[i] = nn2_dm[i - 1] + unknown_dm;
	//	nn3_dm[i] = nn3_dm[i - 1];

	//	for (int j = 0; j < num_domain; j++) {
	//		nn3_dm[i] += nnz_c[i - 1][j];
	//	}
	//}
	
	return 0;
}

int Eqn_Solver_T_E(Element* el, Element* oppo_element, int myid) {

	
	//delete[]nnz_c;


	unordered_map<long long, complex<double>>mat_full;  //Frequency independent matrix terms
	unordered_map<long long, complex<double>>mat_full_ome;  //Matrix terms related to the power of frequency
	unordered_map<long long, complex<double>>mat_full_ome_2;  //Matrix terms related to the quadratic frequency
	//prepare for code
	cout.precision(16);
	long long nnz = 0;  //Number of non zero elements in a sparse matrix 


	nnz_c = new size_t * [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_c[i] = new size_t[num_domain];

	int kk;

	int ndm, nn, mm, ii, jj, node_ii1, node_ii2, node_jj1, node_jj2;
	int nodeglb_ii1, nodeglb_ii2, nodeglb_jj1, nodeglb_jj2, op1, ofn1, edgeii_loc, edgejj_loc, edgeii_E, edgeii_H, edgejj_E, edgejj_H;
	double nx, ny, nz, Ni_dot_ncrosNj, ncrosNi_dot_ncrosNj, ncrosNix1, ncrosNix2, ncrosNiy1, ncrosNiy2, ncrosNiz1, ncrosNiz2;
	double ncrosNjx1, ncrosNjx2, ncrosNjy1, ncrosNjy2, ncrosNjz1, ncrosNjz2, signE_Ni, signE_Nj, signH_Ni, signH_Nj;
	double Li1[3], Li2[3], Lj1[3], Lj2[3], Ni11[3], Ni12[3], zb1[3], zb2[3], zb3[3], rs[3];
	double Nix1, Nix2, Niy1, Niy2, Niz1, Niz2, Njx1, Njx2, Njy1, Njy2, Njz1, Njz2;
	//complex <long double> nVec[3], Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	//complex<long double> EincX, EincY, EincZ, HincX, HincY, HincZ;

	double Ftemp[6][3];




	MKL_Complex16* Me_full = nullptr;
	int* iRow = nullptr;// store the rows of elements per matrix
	int* jCol = nullptr;// store the columns of elements per matrix

	long long temp = 0;

	//end prepare




	//  Matrix Fullfillment

	for (int i = 0; i < num_element_subdomain; i++) {
		for (ii = 0; ii < 6; ii++) {
			Ftemp[ii][0] = 2.0 * (el[i].c[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1] - el[i].d[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve); //The x component of the curl of the basis function
			Ftemp[ii][1] = 2.0 * (el[i].d[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1] - el[i].b[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);//The y component of the curl of the basis function
			Ftemp[ii][2] = 2.0 * (el[i].b[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1] - el[i].c[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);// The z component of the curl of the basis function
		}
		for (ii = 0; ii < 6; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			for (jj = 0; jj < 6; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					//temp = (long long)((edgeii_E - 1) * num_unknown + edgejj_E - 1);
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full[temp].real(mat_full[temp].real() + 1.0 / mur * signE_Ni * signE_Nj * el[i].length[ii] * el[i].length[jj] * el[i].Ve * \
						(Ftemp[ii][0] * Ftemp[jj][0] + Ftemp[ii][1] * Ftemp[jj][1] + Ftemp[ii][2] * Ftemp[jj][2]));
					mat_full_ome[temp] += 0.0;
					mat_full_ome_2[temp] += 0.0;
				}
			}
		}
	}


	for (int i = 0; i < num_element_subdomain; i++) {
		if (el[i].domain - 1 != myid) {
			cout << "error" << endl;
		}
		for (ii = 0; ii < 12; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			node_ii1 = edge_node_local[ii][0]; node_ii2 = edge_node_local[ii][1];
			nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
			Li1[0] = sign_judge(ii, 5) * el[i].b[node_ii1 - 1] * el[i].length[ii];
			Li1[1] = sign_judge(ii, 5) * el[i].c[node_ii1 - 1] * el[i].length[ii];
			Li1[2] = sign_judge(ii, 5) * el[i].d[node_ii1 - 1] * el[i].length[ii];
			Li2[0] = el[i].b[node_ii2 - 1] * el[i].length[ii];
			Li2[1] = el[i].c[node_ii2 - 1] * el[i].length[ii];
			Li2[2] = el[i].d[node_ii2 - 1] * el[i].length[ii];
			for (jj = 0; jj < 12; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				node_jj1 = edge_node_local[jj][0]; node_jj2 = edge_node_local[jj][1];
				nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
				Lj1[0] = sign_judge(jj, 5) * el[i].b[node_jj1 - 1] * el[i].length[jj];
				Lj1[1] = sign_judge(jj, 5) * el[i].c[node_jj1 - 1] * el[i].length[jj];
				Lj1[2] = sign_judge(jj, 5) * el[i].d[node_jj1 - 1] * el[i].length[jj];
				Lj2[0] = el[i].b[node_jj2 - 1] * el[i].length[jj];
				Lj2[1] = el[i].c[node_jj2 - 1] * el[i].length[jj];
				Lj2[2] = el[i].d[node_jj2 - 1] * el[i].length[jj];

				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full_ome_2[temp].real(mat_full_ome_2[temp].real() - el[i].epsil * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);//
					//mat_full_ome[temp] += 0.0;
					mat_full_ome[temp].imag(mat_full_ome[temp].imag() + el[i].sigma * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);   //sigma!=0


					mat_full[temp] += 0.0;
				}
			}
		}
	}



	//Absorbing boundary condition term
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && ((ofn1 == -7) || (ofn1 == -8))) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);


					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
						}
					}
				}
			}
		}
	}

	//LP
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if (ofn1 == -6 || ofn1 == -66) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);

					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;

						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / (Z0 * Wmicrostrip / Hsubstrate));
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;

						}
					}
				}
			}
		}
	}

	int size_temp = mat_full.size();
	//cout << "myid is " << myid << endl;

	//cout << "size_temp is " << size_temp << endl;

	// Fehii Feeii Fehij  Feeij
	// Fehii Feeii Fehij  Feeij

	//Total number of  elements in Stiffness matrix, mass matrix, and ABC matrix

	//contribution from numerical flux, numerical flux has two parts : outgoing fluxand incoming flux
	int count_oppo_element = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		if (mesh_at_RTC[i] != 0) {    //on the RTC boundary
			for (nn = 0; nn < 4; nn++) {
				op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
				if (op1 != 0) {
					//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
					if (el[i].face[nn].whether_boundary) {
						//if (el[op1 - 1].domain != el[i].domain) {
						if (oppo_element[count_oppo_element].Global_num == op1) {
							op1 = count_oppo_element + 1;
							count_oppo_element++;
							//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
						}
						//else {
						//	for (int temp = 0; temp < num_opp_element[myid]; temp++) {
						//		if (oppo_element[temp].Global_num == op1) {
						//			op1 = temp+1;
						//			cout << "wrong" << endl;
						//			count_oppo_element++;
						//			break;
						//		}
						//	}
						//}
						nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
							signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];//
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//outgoing flux
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in outgoing flux：\hat{ n } \times \mathbf{ H } ^ i
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//!second term in outgoing flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ i }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									//This term should be integraed
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//incoming Flux
								//if (count_oppo_element == 1) {
								//	cout << FactorH * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj << endl;
								//	//cout << "ttttttttttttttttttttt" << endl;
								//}
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];

								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];

								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;

								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//first term in incoming flux：\hat{ n } \times \mathbf{ H } ^ j
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//second term in incoming flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ j }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
							}
						}
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_H = el[i].Hedge_GBNO[edgeii_loc - 1];
							signH_Ni = 1.0 * el[i].Hedge_GBNO[edgeii_loc + 12 - 1];
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni;
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//RTC Terms in current subdomain
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in RTC：\hat{ n } \times \mathbf{ E }^ i
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
										/*  mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
										  mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//!second term in RTC : \hat{ n } \times  \mathbf{ H }^ i} \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//RTC terms in neighboring subdomain
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);

								//edgejj_E = el[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								//signE_Nj = 1.0 * el[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								//node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								//nodeglb_jj1 = el[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = el[op1 - 1].node[node_jj2 - 1].ver;
								////x, y, z components of n cross Nj
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//third term in RTC：\hat{n} \times \mathbf{E} ^ j
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									/*	mat_full[temp] += 0;
										mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);*/
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//fourth term in RTC : \hat{n} \times  \mathbf{H}^ i } \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									/*mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									/*	mat_full_ome[temp] += 0;
										mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									mat_full_ome_2[temp] += 0;
								}
							}
						}
					}
				}
			}
		}
	}
	cout << "myid" << " = " << myid << endl;
	//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
	//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;

	nnz = 0;
	matsize = mat_full.size();
	Me_full = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome_2 = new MKL_Complex16[matsize];
	iRow = new int[matsize]; jCol = new int[matsize];
#   pragma omp parallel for 
	for (int i = 0; i < matsize; ++i) {
		Me_full[i].real = (0.0); Me_full[i].imag = (0.0);
		Me_full_ome[i].real = (0.0); Me_full_ome[i].imag = (0.0);
		Me_full_ome_2[i].real = (0.0); Me_full_ome_2[i].imag = (0.0);
		iRow[i] = 0; jCol[i] = 0;
	}

	for (auto& m : mat_full) {
		Me_full[nnz].real = m.second.real();
		Me_full[nnz].imag = m.second.imag();
		iRow[nnz] = m.first / num_unknown + 1;
		jCol[nnz] = m.first % num_unknown + 1;
		++nnz;
	}
	mat_full.clear();
	//cout << "mat_full =  " << nnz << endl;


	//mat_full_ome
	nnz = 0;
	matsize = mat_full_ome.size();
	for (auto& m : mat_full_ome) {
		Me_full_ome[nnz].real = (m.second).real();//?
		Me_full_ome[nnz].imag = (m.second).imag();//?
		++nnz;
	}
	mat_full_ome.clear();


	nnz = 0;
	matsize = mat_full_ome_2.size();
	for (auto& m : mat_full_ome_2) {
		Me_full_ome_2[nnz].real = m.second.real();
		Me_full_ome_2[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full_ome_2.clear();
	//cout << "mat_full_ome_2 = " << nnz << endl;

	int nn2, mm2;

	//delete[]p_pardiso; delete[]ip; delete[]jp; delete[]acoo1; delete[]acoo; delete[]acoo_ome; delete[]acoo_ome_2; delete[]rowind; delete[]colind; delete[]num_nzero_Pmatrix;
	//delete[]m; delete[]m1; delete[]m_ome; delete[]m_ome_2;



	acoo1 = new MKL_Complex16[nnz];
	acoo = new MKL_Complex16[nnz];
	acoo_ome = new MKL_Complex16[nnz];
	acoo_ome_2 = new MKL_Complex16[nnz];

	rowind = new int[nnz];
	colind = new int[nnz];
	num_nzero_Pmatrix = new int[num_domain];


	int nnz_tot = 0;

	m = new MKL_Complex16[nnz];
	m_ome = new MKL_Complex16[nnz];
	m_ome_2 = new MKL_Complex16[nnz];
	m1 = new MKL_Complex16[nnz];

	mrow = new int[nnz];
	mcol = new int[nnz];


	size_t cc = 0;
	int nnz_dm;
	nn = 0;
	for (ndm = myid; ndm < myid + 1; ++ndm) {
		mm = nn + num_unknown_subdomain[ndm][0] + num_unknown_subdomain[ndm][1];
		nn2 = 0;
		for (size_t ndm2 = 0; ndm2 != num_domain; ++ndm2) {
			nnz_dm = 0;
			mm2 = nn2 + num_unknown_subdomain[ndm2][0] + num_unknown_subdomain[ndm2][1];
#      pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow[j] > nn && iRow[j] <= mm && jCol[j] > nn2 && jCol[j] <= mm2) {//?

#	   pragma omp critical
					{
						++nnz_dm;
						m1[cc] = Me_full[j];
						m_ome[cc] = Me_full_ome[j];
						m_ome_2[cc] = Me_full_ome_2[j];

						mrow[cc] = iRow[j] - nn;
						mcol[cc] = jCol[j] - nn2;
						++cc;
					}
					if (ndm == ndm2) {
#	 pragma omp critical
						{

							acoo1[nnz_tot] = Me_full[j];
							acoo_ome[nnz_tot] = Me_full_ome[j];
							acoo_ome_2[nnz_tot] = Me_full_ome_2[j];
							rowind[nnz_tot] = iRow[j] - nn;
							colind[nnz_tot] = jCol[j] - nn2;
							nnz_tot += 1;

						}
					}
				}
			}
			nn2 = mm2;
			nnz_c[ndm][ndm2] = nnz_dm;//The number of data to fill in each matrix
			if (ndm == ndm2) num_nzero_Pmatrix[ndm] = nnz_dm;
		}
		nn = mm;
	}
	delete[] Me_full; delete[] Me_full_ome; delete[] Me_full_ome_2;
	delete[] iRow;
	delete[] jCol;

	int number_unknown_subdomain = num_unknown_subdomain[myid][1] + num_unknown_subdomain[myid][0];
	p_pardiso = new MKL_Complex16[nnz_tot];
	ip = new int[number_unknown_subdomain + 1];
	jp = new int[nnz_tot];



	return 0;
}


int Eqn_Solver_Iteration(Element* el, Element* oppo_element, int myid) {
	unordered_map<long long, complex<double>>mat_full;  //Frequency independent matrix terms
	unordered_map<long long, complex<double>>mat_full_ome;  //Matrix terms related to the power of frequency
	unordered_map<long long, complex<double>>mat_full_ome_2;  //Matrix terms related to the quadratic frequency
	//prepare for code
	cout.precision(16);
	long long nnz = 0;  //Number of non zero elements in a sparse matrix 
	nnz_c = new size_t * [num_domain];//plus
	for (int i = 0; i < num_domain; ++i) nnz_c[i] = new size_t[num_domain];

	int kk;

	int ndm, nn, mm, ii, jj, node_ii1, node_ii2, node_jj1, node_jj2;
	int nodeglb_ii1, nodeglb_ii2, nodeglb_jj1, nodeglb_jj2, op1, ofn1, edgeii_loc, edgejj_loc, edgeii_E, edgeii_H, edgejj_E, edgejj_H;
	double nx, ny, nz, Ni_dot_ncrosNj, ncrosNi_dot_ncrosNj, ncrosNix1, ncrosNix2, ncrosNiy1, ncrosNiy2, ncrosNiz1, ncrosNiz2;
	double ncrosNjx1, ncrosNjx2, ncrosNjy1, ncrosNjy2, ncrosNjz1, ncrosNjz2, signE_Ni, signE_Nj, signH_Ni, signH_Nj;
	double Li1[3], Li2[3], Lj1[3], Lj2[3], Ni11[3], Ni12[3], zb1[3], zb2[3], zb3[3], rs[3];
	double Nix1, Nix2, Niy1, Niy2, Niz1, Niz2, Njx1, Njx2, Njy1, Njy2, Njz1, Njz2;
	//complex <long double> nVec[3], Einc[3], ncrosE[3], ncrosEcrosn[3], Hinc[3], ncrosH[3];
	//complex<long double> EincX, EincY, EincZ, HincX, HincY, HincZ;

	double Ftemp[6][3];




	MKL_Complex16* Me_full = nullptr;
	int* iRow = nullptr;// store the rows of elements per matrix
	int* jCol = nullptr;// store the columns of elements per matrix

	long long temp = 0;

	//end prepare




	//  Matrix Fullfillment

	for (int i = 0; i < num_element_subdomain; i++) {
		for (ii = 0; ii < 6; ii++) {
			Ftemp[ii][0] = 2.0 * (el[i].c[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1] - el[i].d[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve); //The x component of the curl of the basis function
			Ftemp[ii][1] = 2.0 * (el[i].d[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1] - el[i].b[edge_node_local[ii][0] - 1] * el[i].d[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);//The y component of the curl of the basis function
			Ftemp[ii][2] = 2.0 * (el[i].b[edge_node_local[ii][0] - 1] * el[i].c[edge_node_local[ii][1] - 1] - el[i].c[edge_node_local[ii][0] - 1] * el[i].b[edge_node_local[ii][1] - 1]) / (6.0 * el[i].Ve) / (6.0 * el[i].Ve);// The z component of the curl of the basis function
		}
		for (ii = 0; ii < 6; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			for (jj = 0; jj < 6; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					//temp = (long long)((edgeii_E - 1) * num_unknown + edgejj_E - 1);
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full[temp].real(mat_full[temp].real() + 1.0 / mur * signE_Ni * signE_Nj * el[i].length[ii] * el[i].length[jj] * el[i].Ve * \
						(Ftemp[ii][0] * Ftemp[jj][0] + Ftemp[ii][1] * Ftemp[jj][1] + Ftemp[ii][2] * Ftemp[jj][2]));
					mat_full_ome[temp] += 0.0;
					mat_full_ome_2[temp] += 0.0;
				}
			}
		}
	}


	for (int i = 0; i < num_element_subdomain; i++) {
		if (el[i].domain - 1 != myid) {
			cout << "error" << endl;
		}
		for (ii = 0; ii < 12; ii++) {
			edgeii_E = el[i].Eedge_GBNO[ii];
			signE_Ni = 1.0 * el[i].Eedge_GBNO[ii + 12];
			node_ii1 = edge_node_local[ii][0]; node_ii2 = edge_node_local[ii][1];
			nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
			Li1[0] = sign_judge(ii, 5) * el[i].b[node_ii1 - 1] * el[i].length[ii];
			Li1[1] = sign_judge(ii, 5) * el[i].c[node_ii1 - 1] * el[i].length[ii];
			Li1[2] = sign_judge(ii, 5) * el[i].d[node_ii1 - 1] * el[i].length[ii];
			Li2[0] = el[i].b[node_ii2 - 1] * el[i].length[ii];
			Li2[1] = el[i].c[node_ii2 - 1] * el[i].length[ii];
			Li2[2] = el[i].d[node_ii2 - 1] * el[i].length[ii];
			for (jj = 0; jj < 12; jj++) {
				//edgejj_E = el[i].Eedge_GBNO[jj] + Accumulated_unknowns[el[i].domain - 1];
				edgejj_E = el[i].Eedge_GBNO[jj];
				signE_Nj = 1.0 * el[i].Eedge_GBNO[jj + 12];
				node_jj1 = edge_node_local[jj][0]; node_jj2 = edge_node_local[jj][1];
				nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
				Lj1[0] = sign_judge(jj, 5) * el[i].b[node_jj1 - 1] * el[i].length[jj];
				Lj1[1] = sign_judge(jj, 5) * el[i].c[node_jj1 - 1] * el[i].length[jj];
				Lj1[2] = sign_judge(jj, 5) * el[i].d[node_jj1 - 1] * el[i].length[jj];
				Lj2[0] = el[i].b[node_jj2 - 1] * el[i].length[jj];
				Lj2[1] = el[i].c[node_jj2 - 1] * el[i].length[jj];
				Lj2[2] = el[i].d[node_jj2 - 1] * el[i].length[jj];

				if ((edgeii_E != 0) && (edgejj_E != 0)) {
					temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
					mat_full_ome_2[temp].real(mat_full_ome_2[temp].real() - el[i].epsil * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);//
					//mat_full_ome[temp] += 0.0;
					mat_full_ome[temp].imag(mat_full_ome[temp].imag() + el[i].sigma * signE_Ni * signE_Nj * 1.0 / (6.0 * el[i].Ve) / (6.0 * el[i].Ve) * \
						(Li2[0] * (Lj2[0] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[1] * (Lj2[1] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li2[2] * (Lj2[2] * Flabel(nodeglb_ii1, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii1, nodeglb_jj2)) + \
							Li1[0] * (Lj2[0] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[0] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[1] * (Lj2[1] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[1] * Flabel(nodeglb_ii2, nodeglb_jj2)) + \
							Li1[2] * (Lj2[2] * Flabel(nodeglb_ii2, nodeglb_jj1) + Lj1[2] * Flabel(nodeglb_ii2, nodeglb_jj2))) * el[i].Ve / 20.0);   //sigma!=0
					mat_full[temp] += 0.0;
				}
			}
		}
	}



	//Absorbing boundary condition term
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if ((op1 == 0) && ((ofn1 == -7) || (ofn1 == -8))) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);


					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / el[i].eta);
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;
						}
					}
				}
			}
		}
	}

	//Absorbing boundary condition term
	for (int i = 0; i < num_element_subdomain; ++i) {
		for (nn = 0; nn < 4; nn++) {
			op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
			if (ofn1 == -6 || ofn1 == -66) {
				nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
				for (ii = 0; ii < 6; ii++) {
					edgeii_loc = face_edge[nn][ii];
					edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
					signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];
					node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
					nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
					//n cross Ni
					ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
					ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
					ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);

					for (jj = 0; jj < 6; jj++) {
						edgejj_loc = face_edge[nn][jj];
						//edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1] + Accumulated_unknowns[el[i].domain - 1];
						edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1];
						signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1];
						node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
						nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
						//n cross Nj
						ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
						ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
						ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);

						ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
							ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
							ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;

						if ((edgeii_E != 0) && (edgejj_E != 0)) {
							temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
							mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj / (Z0 * Wmicrostrip / Hsubstrate));
							mat_full_ome_2[temp] += 0.0;
							mat_full[temp] += 0.0;

						}
					}
				}
			}
		}
	}

	int size_temp = mat_full.size();
	//cout << "myid is " << myid << endl;

	//cout << "size_temp is " << size_temp << endl;

	// Fehii Feeii Fehij  Feeij
	// Fehii Feeii Fehij  Feeij

	//Total number of  elements in Stiffness matrix, mass matrix, and ABC matrix

	//contribution from numerical flux, numerical flux has two parts : outgoing fluxand incoming flux
	int count_oppo_element = 0;
	for (int i = 0; i < num_element_subdomain; i++) {
		if (mesh_at_RTC[i] != 0) {    //on the RTC boundary
			for (nn = 0; nn < 4; nn++) {
				op1 = el[i].face[nn].opp[0]; ofn1 = el[i].face[nn].opp[1];
				if (op1 != 0) {
					//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
					if (el[i].face[nn].whether_boundary) {
						//if (el[op1 - 1].domain != el[i].domain) {
						if (oppo_element[count_oppo_element].Global_num == op1) {
							op1 = count_oppo_element + 1;
							count_oppo_element++;
							//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;
						}
						//else {
						//	for (int temp = 0; temp < num_opp_element[myid]; temp++) {
						//		if (oppo_element[temp].Global_num == op1) {
						//			op1 = temp+1;
						//			cout << "wrong" << endl;
						//			count_oppo_element++;
						//			break;
						//		}
						//	}
						//}
						nx = el[i].face[nn].N_Vector[0]; ny = el[i].face[nn].N_Vector[1]; nz = el[i].face[nn].N_Vector[2];
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_E = el[i].Eedge_GBNO[edgeii_loc - 1];
							signE_Ni = 1.0 * el[i].Eedge_GBNO[edgeii_loc + 12 - 1];//
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//outgoing flux
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in outgoing flux：\hat{ n } \times \mathbf{ H } ^ i
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//!second term in outgoing flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ i }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									//This term should be integraed
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//incoming Flux

								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];

								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];

								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;

								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//first term in incoming flux：\hat{ n } \times \mathbf{ H } ^ j
								if ((edgeii_E != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - FactorH * signE_Ni * signH_Nj * Ni_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
								//second term in incoming flux : \hat{ n } \times{ \hat{n} \times \mathbf{E} ^ j }
								if ((edgeii_E != 0) && (edgejj_E != 0)) {
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_E - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];

									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1.0 / (el[i].eta + oppo_element[op1 - 1].eta) * signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + signE_Ni * signE_Nj * ncrosNi_dot_ncrosNj);
									mat_full_ome_2[temp] += 0.0;
									mat_full[temp] += 0.0;
								}
							}
						}
						for (ii = 0; ii < 6; ii++) {
							edgeii_loc = face_edge[nn][ii];
							edgeii_H = el[i].Hedge_GBNO[edgeii_loc - 1];
							signH_Ni = 1.0 * el[i].Hedge_GBNO[edgeii_loc + 12 - 1];
							node_ii1 = edge_node_local[edgeii_loc - 1][0]; node_ii2 = edge_node_local[edgeii_loc - 1][1];
							nodeglb_ii1 = el[i].node[node_ii1 - 1].ver; nodeglb_ii2 = el[i].node[node_ii2 - 1].ver;
							//x, y, z components of Ni
							Nix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii2 - 1]; Nix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].b[node_ii1 - 1];
							Niy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii2 - 1]; Niy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].c[node_ii1 - 1];
							Niz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii2 - 1]; Niz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * el[i].d[node_ii1 - 1];
							//x, y, z components of n cross Ni;
							ncrosNix1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii2 - 1] - nz * el[i].c[node_ii2 - 1]); ncrosNix2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_ii1 - 1] - nz * el[i].c[node_ii1 - 1]);
							ncrosNiy1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii2 - 1] - nx * el[i].d[node_ii2 - 1]); ncrosNiy2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_ii1 - 1] - nx * el[i].d[node_ii1 - 1]);
							ncrosNiz1 = el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii2 - 1] - ny * el[i].b[node_ii2 - 1]); ncrosNiz2 = sign_judge(edgeii_loc, 6) * el[i].length[edgeii_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_ii1 - 1] - ny * el[i].b[node_ii1 - 1]);
							for (jj = 0; jj < 6; jj++) {
								//RTC Terms in current subdomain
								edgejj_loc = face_edge[nn][jj];
								edgejj_E = el[i].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[i].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * el[i].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[i].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = el[i].node[node_jj1 - 1].ver; nodeglb_jj2 = el[i].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj2 - 1] - nz * el[i].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (ny * el[i].d[node_jj1 - 1] - nz * el[i].c[node_jj1 - 1]);
								ncrosNjy1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj2 - 1] - nx * el[i].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nz * el[i].b[node_jj1 - 1] - nx * el[i].d[node_jj1 - 1]);
								ncrosNjz1 = el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj2 - 1] - ny * el[i].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[i].length[edgejj_loc - 1] / (6.0 * el[i].Ve) * (nx * el[i].c[node_jj1 - 1] - ny * el[i].b[node_jj1 - 1]);
								//first term in RTC：\hat{ n } \times \mathbf{ E }^ i
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[el[i].domain - 1];
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
										/*  mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
										  mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//!second term in RTC : \hat{ n } \times  \mathbf{ H }^ i} \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[el[i].domain - 1];
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * el[i].eta * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//RTC terms in neighboring subdomain
								edgejj_loc = face_edge[ofn1 - 1][jj];

								edgejj_E = oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								signE_Nj = 1.0 * oppo_element[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * oppo_element[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								nodeglb_jj1 = oppo_element[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = oppo_element[op1 - 1].node[node_jj2 - 1].ver;
								//x, y, z components of n cross Nj
								ncrosNjx1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj2 - 1] - nz * oppo_element[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (ny * oppo_element[op1 - 1].d[node_jj1 - 1] - nz * oppo_element[op1 - 1].c[node_jj1 - 1]);
								ncrosNjy1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj2 - 1] - nx * oppo_element[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nz * oppo_element[op1 - 1].b[node_jj1 - 1] - nx * oppo_element[op1 - 1].d[node_jj1 - 1]);
								ncrosNjz1 = oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj2 - 1] - ny * oppo_element[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * oppo_element[op1 - 1].length[edgejj_loc - 1] / (6.0 * oppo_element[op1 - 1].Ve) * (nx * oppo_element[op1 - 1].c[node_jj1 - 1] - ny * oppo_element[op1 - 1].b[node_jj1 - 1]);

								//edgejj_E = el[op1 - 1].Eedge_GBNO[edgejj_loc - 1]; edgejj_H = el[op1 - 1].Hedge_GBNO[edgejj_loc - 1];
								//signE_Nj = 1.0 * el[op1 - 1].Eedge_GBNO[edgejj_loc + 12 - 1]; signH_Nj = 1.0 * el[op1 - 1].Hedge_GBNO[edgejj_loc + 12 - 1];
								//node_jj1 = edge_node_local[edgejj_loc - 1][0]; node_jj2 = edge_node_local[edgejj_loc - 1][1];
								//nodeglb_jj1 = el[op1 - 1].node[node_jj1 - 1].ver; nodeglb_jj2 = el[op1 - 1].node[node_jj2 - 1].ver;
								////x, y, z components of n cross Nj
								//ncrosNjx1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj2 - 1] - nz * el[op1 - 1].c[node_jj2 - 1]); ncrosNjx2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (ny * el[op1 - 1].d[node_jj1 - 1] - nz * el[op1 - 1].c[node_jj1 - 1]);
								//ncrosNjy1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj2 - 1] - nx * el[op1 - 1].d[node_jj2 - 1]); ncrosNjy2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nz * el[op1 - 1].b[node_jj1 - 1] - nx * el[op1 - 1].d[node_jj1 - 1]);
								//ncrosNjz1 = el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj2 - 1] - ny * el[op1 - 1].b[node_jj2 - 1]); ncrosNjz2 = sign_judge(edgejj_loc, 6) * el[op1 - 1].length[edgejj_loc - 1] / (6.0 * el[op1 - 1].Ve) * (nx * el[op1 - 1].c[node_jj1 - 1] - ny * el[op1 - 1].b[node_jj1 - 1]);
								//third term in RTC：\hat{n} \times \mathbf{E} ^ j
								if ((edgeii_H != 0) && (edgejj_E != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_E - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * FactorH * signH_Ni * signE_Nj * Ni_dot_ncrosNj);
									mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() - 1 * el[i].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
										/*	mat_full[temp].real(mat_full[temp].real() - 1 * el[i].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * Ni_dot_ncrosNj);
											mat_full_ome[temp] += 0;*/
									mat_full_ome_2[temp] += 0;
								}
								//fourth term in RTC : \hat{n} \times  \mathbf{H}^ i } \times{ \hat{n}
								if ((edgeii_H != 0) && (edgejj_H != 0)) {
									Ni_dot_ncrosNj = (Nix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Nix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Nix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Nix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										Niz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + Niz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										Niz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + Niz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									ncrosNi_dot_ncrosNj = (ncrosNix1 * ncrosNjx1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNix1 * ncrosNjx2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNix2 * ncrosNjx1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNix2 * ncrosNjx2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiy1 * ncrosNjy1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiy1 * ncrosNjy2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiy2 * ncrosNjy1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiy2 * ncrosNjy2 * Flabel(nodeglb_ii2, nodeglb_jj2) + \
										ncrosNiz1 * ncrosNjz1 * Flabel(nodeglb_ii1, nodeglb_jj1) + ncrosNiz1 * ncrosNjz2 * Flabel(nodeglb_ii1, nodeglb_jj2) + \
										ncrosNiz2 * ncrosNjz1 * Flabel(nodeglb_ii2, nodeglb_jj1) + ncrosNiz2 * ncrosNjz2 * Flabel(nodeglb_ii2, nodeglb_jj2)) * el[i].face[nn].Area / 12.0;
									temp = (long long)(edgeii_H - 1) * num_unknown + edgejj_H - 1 + Accumulated_unknowns[oppo_element[op1 - 1].domain - 1];
									//mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);
									mat_full[temp] += 0;
									mat_full_ome[temp].imag(mat_full_ome[temp].imag() + 1 * FactorH * el[i].eta * oppo_element[op1 - 1].eta / (el[i].eta + oppo_element[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);

									/*	mat_full_ome[temp] += 0;
										mat_full[temp].real(mat_full[temp].real() + 1 * el[i].eta * el[op1 - 1].eta / (el[i].eta + el[op1 - 1].eta) * signH_Ni * signH_Nj * ncrosNi_dot_ncrosNj);*/
									mat_full_ome_2[temp] += 0;
								}
							}
						}
					}
				}
			}
		}
	}
	//cout << "myid" << " = " << myid << endl;
	//cout << "yes and count_oppo_element_temp" << " = " << count_oppo_element_temp << endl;
	//cout << "yes and count_oppo_element" << " = " << count_oppo_element << endl;

	nnz = 0;
	matsize = mat_full.size();
	Me_full = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome = new MKL_Complex16[matsize];
	MKL_Complex16* Me_full_ome_2 = new MKL_Complex16[matsize];
	iRow = new int[matsize]; jCol = new int[matsize];
#   pragma omp parallel for 
	for (int i = 0; i < matsize; ++i) {
		Me_full[i].real = (0.0); Me_full[i].imag = (0.0);
		Me_full_ome[i].real = (0.0); Me_full_ome[i].imag = (0.0);
		Me_full_ome_2[i].real = (0.0); Me_full_ome_2[i].imag = (0.0);
		iRow[i] = 0; jCol[i] = 0;
	}

	for (auto& m : mat_full) {
		Me_full[nnz].real = m.second.real();
		Me_full[nnz].imag = m.second.imag();
		iRow[nnz] = m.first / num_unknown + 1;
		jCol[nnz] = m.first % num_unknown + 1;
		++nnz;
	}
	mat_full.clear();
	//cout << "mat_full =  " << nnz << endl;


	//mat_full_ome
	nnz = 0;
	matsize = mat_full_ome.size();
	for (auto& m : mat_full_ome) {
		Me_full_ome[nnz].real = (m.second).real();//?
		Me_full_ome[nnz].imag = (m.second).imag();//?
		++nnz;
	}
	mat_full_ome.clear();


	nnz = 0;
	matsize = mat_full_ome_2.size();
	for (auto& m : mat_full_ome_2) {
		Me_full_ome_2[nnz].real = m.second.real();
		Me_full_ome_2[nnz].imag = m.second.imag();
		++nnz;
	}
	mat_full_ome_2.clear();
	//cout << "mat_full_ome_2 = " << nnz << endl;

	int nn2, mm2;
	acoo1 = new MKL_Complex16[nnz];
	acoo = new MKL_Complex16[nnz];
	acoo_ome = new MKL_Complex16[nnz];
	acoo_ome_2 = new MKL_Complex16[nnz];

	rowind = new int[nnz];
	colind = new int[nnz];

	num_nzero_Pmatrix = new int[num_domain];
	//for (int i = 0; i < num_domain; ++i) {
	//	num_nzero_Pmatrix[i] = 0;
	//}

	int nnz_tot = 0;

	m = new MKL_Complex16[nnz];
	m_ome = new MKL_Complex16[nnz];
	m_ome_2 = new MKL_Complex16[nnz];
	m1 = new MKL_Complex16[nnz];

	mrow = new int[nnz];
	mcol = new int[nnz];


	size_t cc = 0;
	int nnz_dm;
	nn = 0;
	for (ndm = myid; ndm < myid + 1; ++ndm) {
		mm = nn + num_unknown_subdomain[ndm][0] + num_unknown_subdomain[ndm][1];
		nn2 = 0;
		for (size_t ndm2 = 0; ndm2 != num_domain; ++ndm2) {
			nnz_dm = 0;
			mm2 = nn2 + num_unknown_subdomain[ndm2][0] + num_unknown_subdomain[ndm2][1];
#      pragma omp parallel for 
			for (int j = 0; j < nnz; j++) {
				if (iRow[j] > nn && iRow[j] <= mm && jCol[j] > nn2 && jCol[j] <= mm2) {//?

#	   pragma omp critical
					{
						++nnz_dm;
						m1[cc] = Me_full[j];
						m_ome[cc] = Me_full_ome[j];
						m_ome_2[cc] = Me_full_ome_2[j];

						mrow[cc] = iRow[j] - nn;
						mcol[cc] = jCol[j] - nn2;
						++cc;
					}
					if (ndm == ndm2) {
#	 pragma omp critical
						{

							acoo1[nnz_tot] = Me_full[j];
							acoo_ome[nnz_tot] = Me_full_ome[j];
							acoo_ome_2[nnz_tot] = Me_full_ome_2[j];
							rowind[nnz_tot] = iRow[j] - nn;
							colind[nnz_tot] = jCol[j] - nn2;
							nnz_tot += 1;

						}
					}
				}
			}
			nn2 = mm2;
			nnz_c[ndm][ndm2] = nnz_dm;//The number of data to fill in each matrix
			if (ndm == ndm2) num_nzero_Pmatrix[ndm] = nnz_dm;
		}
		nn = mm;
	}
	delete[] Me_full; delete[] Me_full_ome; delete[] Me_full_ome_2;
	delete[] iRow;
	delete[] jCol;
	cout << "cc = " << cc << endl;

	int number_unknown_subdomain = num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][0];
	p_pardiso = new MKL_Complex16[nnz_tot];
	ip = new int[number_unknown_subdomain + 1];
	jp = new int[nnz_tot];

	cout << "myid is " << myid << endl;
	cout << "nnz_tot = " << nnz_tot << endl;
	//cout << "num_nzero_Pmatrix[myid] = " << num_nzero_Pmatrix[myid] << endl;

	/*int nzero1 = 0; int row1 = 0;*/

	//nnz_tot = 0;
	//nn2 = 0;

	row_dm = new int[num_domain];
	nzero_dm = new int[num_domain];
	nn2_dm = new int[num_domain];
	nn3_dm = new int[num_domain];

	n3_dm = new int[num_domain];

	//row_dm[0] = 0; nzero_dm[0] = 0; nn2_dm[0] = 0; nn3_dm[0] = 0; r_dm[0] = 0; n3_dm[0] = 0;
	//int unknown_dm;
	//for (int i = 1; i < num_domain; ++i) {
	//	unknown_dm = num_unknown_subdomain[i - 1][0] + num_unknown_subdomain[i - 1][1];
	//	nnz_dm = num_nzero_Pmatrix[i - 1];

	//	r_dm[i] = r_dm[i - 1] + Edge_Boundary[i - 1];
	//	n3_dm[i] = n3_dm[i - 1] + Edge_Boundary[i - 1];

	//	row_dm[i] = row_dm[i - 1] + unknown_dm + 1;//第i个domain的row坐标
	//	nzero_dm[i] = nzero_dm[i - 1] + nnz_dm;//第i个domain的坐标  统计的是P矩阵的nnz数目
	//	nn2_dm[i] = nn2_dm[i - 1] + unknown_dm;
	//	nn3_dm[i] = nn3_dm[i - 1];

	//	for (int j = 0; j < num_domain; j++) {
	//		nn3_dm[i] += nnz_c[i - 1][j];
	//	}
	//}

	return 0;
}

	double Flabel(int vnode1, int vnode2) {
		double Flable0;
		if (vnode1 == vnode2) {
			Flable0 = 2.0;
		}
		else {
			Flable0 = 1.0;
		}
		return Flable0;
	}
	double sign_judge(int edgeNo1, int edgeNo2) {
		double sign_judge0;
		if (edgeNo1 <= edgeNo2) {
			sign_judge0 = -1.0;
		}
		else {
			sign_judge0 = 1.0;
		}
		return sign_judge0;
	}
