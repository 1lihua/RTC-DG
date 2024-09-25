#include<iostream>
#include<vector>
#include"mkl.h"
#include <fstream>
#include"ReadFile.h"
#include"Opp_Sct.h"
#include"Face_Vector.h"
#include"General_Parameter.h"
#include"Epsil_Mu.h"
#include"Global_Data.h"
#include"LookUp_EdgeIndex.h"
#include"Eqn_Solver.h"
#include"OutField.h"
#include<algorithm>
#include"time.h"
#include "metis1.h"
#include<mpi.h>
#include<omp.h>
#include"Solve_involve_Freq.h"
#include "Lookup_unknown_T_q.h"
#include "Material_Parm_T_q.h"
#include "Matrix_Generator_T_q.h"
#include "Matrix_Generator_stru.h"
#include"read_input.h"
using namespace std;

int main(int argc, char** argv) {

	time_t start, end;

	start = time(NULL);

	int myid, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	//const int myid = myid1;
	if (numprocs == 1)
	{
		cout << "Can't be one process" << endl;
		return 0;
	}
	else {
		num_domain = numprocs;
	}
	
#pragma omp parallel
	{
		process_threads = omp_get_num_threads();
	}
	if (myid == 0) {
		cout << "process_threads = " << process_threads << endl;
	}
	char input_file1[512] ="inputfile.txt";

	char* input_file = input_file1;

	//string input_file="inputfile.txt";
	//FILE* fp = fopen(input_file, "r");
	Read_inputfile(input_file);

	if (myid == 0) {
		cout << "  /***********************************/            " << endl;

		cout << "SI_meter = " << SI_meter << endl;
		cout << "fre_start = " << fre_start << endl;
		cout << "fre_step = " << fre_step << endl;
		cout << "fre_end = " << fre_end << endl;

		cout << "wave_material = " << wave_material << endl;
		cout << "num_wave = " << num_wave << endl;
		cout << "rb_wave = " << rb_wave << endl;
		cout << "ra_wave = " << ra_wave << endl;

		for (int i = 0; i < num_wave; i++) {
			cout << "x_wave[" << i << "] = " << x_wave[i] << endl;
		}

		for (int i = 0; i < num_wave; i++) {
			cout << "Amplitude_wave[" << i << "] = " << Amplitude_wave[i] << endl;
		}
		for (int i = 0; i < 2; i++) {
			cout << "DIEeps[" << i << "] = " << DIEeps[i] << endl;
		}
		for (int i = 0; i < 2; i++) {
			cout << "DIEalpha[" << i << "] = " << DIEalpha[i] << endl;
		}
		cout << "  /***********************************/            " << endl;
	}





	Element* el;
	if (myid == 0) {
		metis11();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	BoundaryRead();
	el = (Element*)ReadFile();
	if (myid == 0) {
		cout << "FileEnd" << endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);
	num_element_subdomain = 0;

	for (int i = 0; i < num_element; i++) {
		if (el[i].domain == myid + 1) {
			num_element_subdomain++;
		}
	}
	Element* el_subdomain = (Element*)malloc(num_element_subdomain * sizeof(Element));
	int num_count = 0;
	for (int i = 0; i < num_element; i++) {
		if (el[i].domain == myid + 1) {
			el_subdomain[num_count] = el[i];
			num_count++;
		}
	}
	num_opp_element = new int[num_domain];
	for (int i = 0; i < num_domain; i++) {
		num_opp_element[i] = 0;
	}
	Opp_Sct(el_subdomain, num_element_subdomain,myid);
	int temp_num_opp_element = num_opp_element[myid];
	MPI_Gather(&temp_num_opp_element, 1, MPI_INT, num_opp_element, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(num_opp_element, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	int num_element_boundary_total = 0;
	int* address_offset_face= new int[num_domain];
	address_offset_face[0] = 0;
	for (int i = 0; i < num_domain; i++) {
		num_element_boundary_total += num_opp_element[i];
	}
	if (myid == 0) {
		cout << "num_element_boundary_total = " << num_element_boundary_total << endl;
	}
	for (int i = 1; i < num_domain; i++) {
		address_offset_face[i] = address_offset_face[i-1]+num_opp_element[i-1];
	}

	int *Global_num_element_boundary_total = new int[num_element_boundary_total];
	int *Local_num_face_boundary_total = new int[num_element_boundary_total];


	//for (int i = 0; i < num_element_boundary_total; i++) {
	//	Global_num_element_boundary_total[i]=0;
	//	Local_num_face_boundary_total[i] = 0;
	//}


	MPI_Gatherv(Global_num_element_boundary, num_opp_element[myid], MPI_INT, Global_num_element_boundary_total, num_opp_element, address_offset_face, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(Local_num_face_boundary, num_opp_element[myid], MPI_INT, Local_num_face_boundary_total, num_opp_element, address_offset_face, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(Global_num_element_boundary_total, num_element_boundary_total, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Local_num_face_boundary_total, num_element_boundary_total, MPI_INT, 0, MPI_COMM_WORLD);


	MPI_Barrier(MPI_COMM_WORLD);
	Opp_Sct1(el, num_element_boundary_total, myid, Global_num_element_boundary_total, Local_num_face_boundary_total, address_offset_face);


	if (myid == 0) {
		cout << "num_node_boundary = " << num_node_boundary << endl;
		cout << "num_element_boundary=" << num_element_boundary << endl;
		cout << "num_node=" << num_node << endl;
		cout << "num_element=" << num_element << endl;
	}
	cout << "myid is " << myid << "num oppo is " << num_opp_element[myid] << endl;

	
#   pragma omp parallel for
	for (int ith = 0; ith < num_element; ++ith) {
		for (int j = 0; j < 4; ++j) {
			el[ith].node[j].zb[zbX] *= SI_meter;
			el[ith].node[j].zb[zbY] *= SI_meter;
			el[ith].node[j].zb[zbZ] *= SI_meter;
		}
	}
#   pragma omp parallel for
	for (int ith = 0; ith < num_element_subdomain; ++ith) {
		for (int j = 0; j < 4; ++j) {
			el_subdomain[ith].node[j].zb[zbX] *= SI_meter;
			el_subdomain[ith].node[j].zb[zbY] *= SI_meter;
			el_subdomain[ith].node[j].zb[zbZ] *= SI_meter;
		}
	}



	/*   */
	MPI_Barrier(MPI_COMM_WORLD);
	end = time(NULL);
	double time1 = (double)(end - start);
	if (myid == 0) {
		cout << "time1 is " << time1 << endl;
	}

	/*          */

	for (int i = 0; i < num_element_subdomain; i++) {
		for (int nn = 0; nn < 4; nn++) {
			int op0 = el_subdomain[i].face[nn].opp[0]; int ofn = el_subdomain[i].face[nn].opp[1];
			if (op0 == 0 && ofn == 0) {
				int Global_num = el_subdomain[i].Global_num-1;
				int op00 = el[Global_num].face[nn].opp[0];
				int op11 = el[Global_num].face[nn].opp[1];
				el_subdomain[i].face[nn].opp[0] = op00;
				el_subdomain[i].face[nn].opp[1] = op11;
			}
		}
	}

	/*     */
	int* oppo_domain = new int[num_opp_element[myid]];
	Element* oppo_element =  (Element*)malloc(num_opp_element[myid] * sizeof(Element));

	int count=0;
	for (int i = 0; i < num_element; i++) {
		if (el[i].domain - 1 == myid) {
			for (int nn = 0; nn < 4; nn++) {
				int op0 = el[i].face[nn].opp[0]; int ofn = el[i].face[nn].opp[1];
				if (op0 != 0) {
					if (el[op0 - 1].domain - 1 != myid) {//Interface between two sub domains
						oppo_element[count] = el[op0 - 1];
						oppo_domain[count] = el[op0 - 1].domain;
						count++;
					}
				}
			}
		}
	}
	/*     */

	MPI_Barrier(MPI_COMM_WORLD);
	Face_Vector(el_subdomain,num_element_subdomain);
	General_Parameter(el_subdomain, num_element_subdomain);
	Epsil_Mu(el_subdomain,num_element_subdomain);
	//Material_Parm_T_q(el_subdomain, num_element_subdomain);


	int ttt = num_opp_element[myid];
	Face_Vector(oppo_element, ttt);
	General_Parameter(oppo_element, ttt);
	Epsil_Mu(oppo_element, ttt);
	//Material_Parm_T_q(oppo_element, ttt);



	LookUp_EdgeIndex(el_subdomain);


	Vertice_Boundary_T_q = new int[num_domain];
	num_unknown_subdomain_T_q = new int* [num_domain];
	for (int i = 0; i < num_domain; i++) {
		num_unknown_subdomain_T_q[i] = new int[2];
	}

	for (int i = 0; i < num_domain; i++) {
		num_unknown_subdomain_T_q[i][0] = 0;
		num_unknown_subdomain_T_q[i][1] = 0;
		Vertice_Boundary_T_q[i] = 0;
	}

	unknownIndex_Thermal(el_subdomain, myid);



	/*      thermal stress          */
	Node_Boundary_stru = new int[num_domain];
	num_unKnown_subdomain_stru = new int* [num_domain];
	for (int i = 0; i < num_domain; i++) {
		num_unKnown_subdomain_stru[i] = new int[2];
	}
	for (int i = 0; i < num_domain; i++) {
		num_unKnown_subdomain_stru[i][0] = 0;
		num_unKnown_subdomain_stru[i][1] = 0;
		Node_Boundary_stru[i] = 0;
	}
	unknownIndex_stru(el_subdomain, myid);
	/*         thermal stress     */










	MPI_Barrier(MPI_COMM_WORLD);
	/*   gain the information of opposite element T and q */
	MPI_Status status;
	int tag = 0;
	int tag1 = 0;
	int tag2 = 0;
	int tag3 = 0;
	int tag4 = 0;
	int tag5 = 0;
	int tag6 = 0;

	int myid_comm = 1;
	int Global_number;
	int temp_E[24];
	int temp_H[24];
	int temp_T[4];
	int temp_q[4];
	int temp_T2[10];
	int temp_q2[10];

	int temp_uvw[12];
	int temp_abc[12];


	for (int i = 0; i < num_domain; i++) {
		tag1 += num_opp_element[i];
	}
	tag2 = tag1 * 2;
	tag3 = tag1 * 3;
	tag4 = tag1 * 4;
	tag5 = tag1 * 5;
	tag6 = tag1 * 6;
	for (int i = 0; i < num_domain; i++) {
		
		for (int j = 0; j < num_opp_element[i]; j++) {
			myid_comm = 10000;
			MPI_Barrier(MPI_COMM_WORLD);
			tag3++;
			tag4++;
			tag++;
			tag1++;
			tag2++;
			tag5++;
			tag6++;
			if (myid == i) {
				myid_comm = oppo_domain[j] - 1;  // the id of opposite element
				Global_number = oppo_element[j].Global_num;// the global number of opposite element
			}
			MPI_Bcast(&myid_comm, 1, MPI_INT, i, MPI_COMM_WORLD);
			MPI_Bcast(&Global_number, 1, MPI_INT, i, MPI_COMM_WORLD);


			if (i == myid) {
				MPI_Recv(temp_T, 4, MPI_INT, oppo_domain[j] - 1, tag, MPI_COMM_WORLD, &status);
				MPI_Recv(temp_q, 4, MPI_INT, oppo_domain[j] - 1, tag1, MPI_COMM_WORLD, &status);
				for (int ii = 0; ii < 4; ii++) {
					oppo_element[j].node[ii].unknown_T = temp_T[ii];
					oppo_element[j].node[ii].unknown_q = temp_q[ii];

				}

				//MPI_Recv(temp_T2, 10, MPI_INT, oppo_domain[j] - 1, tag, MPI_COMM_WORLD, &status);
				//MPI_Recv(temp_q2, 10, MPI_INT, oppo_domain[j] - 1, tag1, MPI_COMM_WORLD, &status);
				//for (int ii = 0; ii < 10; ii++) {
				//	oppo_element[j].Eedge_T[ii] = temp_T2[ii];
				//	oppo_element[j].Eedge_q[ii] = temp_q2[ii];

				//}

				MPI_Recv(temp_E, 24, MPI_INT, oppo_domain[j] - 1, tag3, MPI_COMM_WORLD, &status);
				MPI_Recv(temp_H, 24, MPI_INT, oppo_domain[j] - 1, tag4, MPI_COMM_WORLD, &status);
				for (int ii = 0; ii < 24; ii++) {
					oppo_element[j].Eedge_GBNO[ii] = temp_E[ii];
					oppo_element[j].Hedge_GBNO[ii] = temp_H[ii];
				}

				MPI_Recv(temp_uvw, 12, MPI_INT, oppo_domain[j] - 1, tag5, MPI_COMM_WORLD, &status);
				MPI_Recv(temp_abc, 12, MPI_INT, oppo_domain[j] - 1, tag6, MPI_COMM_WORLD, &status);
				for (int ii = 0; ii < 4; ii++) {
					for (int jj = 0; jj < 3; jj++) {
						oppo_element[j].node[ii].unknown_uvw[jj] = temp_uvw[ii * 3 + jj];
						oppo_element[j].node[ii].unknown_abc[jj] = temp_abc[ii * 3 + jj];
					}
				}



			}
			if (myid_comm == myid) {
				for (int k = 0; k < num_element_subdomain; k++) {
					if (el_subdomain[k].Global_num == Global_number) {
						for (int ii = 0; ii < 4; ii++) {
							temp_T[ii] = el_subdomain[k].node[ii].unknown_T;
							temp_q[ii] = el_subdomain[k].node[ii].unknown_q;
						}
						MPI_Send(temp_T, 4, MPI_INT, i, tag, MPI_COMM_WORLD);
						MPI_Send(temp_q, 4, MPI_INT, i, tag1, MPI_COMM_WORLD);

						//for (int ii = 0; ii < 10; ii++) {
						//	temp_T2[ii] = el_subdomain[k].Eedge_T[ii];
						//	temp_q2[ii] = el_subdomain[k].Eedge_q[ii];
						//}
						//MPI_Send(temp_T2, 10, MPI_INT, i, tag, MPI_COMM_WORLD);
						//MPI_Send(temp_q2, 10, MPI_INT, i, tag1, MPI_COMM_WORLD);


						for (int ii = 0; ii < 24; ii++) {
							temp_E[ii] = el_subdomain[k].Eedge_GBNO[ii];
							temp_H[ii] = el_subdomain[k].Hedge_GBNO[ii];
						}
						MPI_Send(temp_E, 24, MPI_INT, i, tag3, MPI_COMM_WORLD);
						MPI_Send(temp_H, 24, MPI_INT, i, tag4, MPI_COMM_WORLD);


						for (int ii = 0; ii < 4; ii++) {
							for (int jj = 0; jj < 3; jj++) {
								temp_uvw[ii * 3 + jj] = el_subdomain[k].node[ii].unknown_uvw[jj];
								temp_abc[ii * 3 + jj] = el_subdomain[k].node[ii].unknown_abc[jj];
							}
						}
						MPI_Send(temp_uvw, 12, MPI_INT, i, tag5, MPI_COMM_WORLD);
						MPI_Send(temp_abc, 12, MPI_INT, i, tag6, MPI_COMM_WORLD);

					}
				}

			}
		}
	}
	/*  gain the information of opposite element T and q */


	int* num_unknown_temp_T_q = new int[num_domain * 2];
	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unknown_temp_T_q[i * 2 + k] = num_unknown_subdomain_T_q[i][k];

		}
	}
	int arr_T_q[2]; arr_T_q[0] = 0; arr_T_q[1] = 0;
	for (int i = 0; i < num_domain; i++) {
		if (i == myid) {
			arr_T_q[0] = num_unknown_subdomain_T_q[i][0];
			arr_T_q[1] = num_unknown_subdomain_T_q[i][1];

		}
	}

	MPI_Gather(arr_T_q, 2, MPI_INT, num_unknown_temp_T_q, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(num_unknown_temp_T_q, num_domain * 2, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unknown_subdomain_T_q[i][k] = num_unknown_temp_T_q[i * 2 + k];

		}
	}
	num_unknown_T_q = 0; num_unKnown_Thermal = 0;
	for (int i = 0; i < num_domain; i++) {
		num_unknown_T_q += num_unknown_subdomain_T_q[i][0] + num_unknown_subdomain_T_q[i][1];
		num_unKnown_Thermal += num_unknown_subdomain_T_q[i][0] + num_unknown_subdomain_T_q[i][1];
	}

	Accumulated_unknowns_T_q = new int[num_domain];
	Accumulated_unknowns_T_q[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		Accumulated_unknowns_T_q[i] = Accumulated_unknowns_T_q[i - 1] + num_unknown_subdomain_T_q[i - 1][0] + num_unknown_subdomain_T_q[i - 1][1];
	}
	//cout << "myid is " << myid << "num un is " << num_unknown_subdomain_T_q[myid][0] + num_unknown_subdomain_T_q[myid][1] << endl;

	if (myid == 0) {
		int num_un_total_T_q = 0;
		num_un_total_T_q = Accumulated_unknowns_T_q[num_domain - 1] + num_unknown_subdomain_T_q[num_domain - 1][0] + num_unknown_subdomain_T_q[num_domain - 1][1];
		cout << "num_un_total_T_q   " << num_un_total_T_q << endl;
	}

	int temp_boundary_T_q = Vertice_Boundary_T_q[myid];
	MPI_Gather(&temp_boundary_T_q, 1, MPI_INT, Vertice_Boundary_T_q, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Vertice_Boundary_T_q, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		int num_un_total = 0; int num_unk_boundary = 0;
		num_un_total = Accumulated_unknowns_T_q[num_domain - 1] + num_unknown_subdomain_T_q[num_domain - 1][0] + num_unknown_subdomain_T_q[num_domain - 1][1];
		for (int i = 0; i < num_domain; i++) {
			num_unk_boundary += Vertice_Boundary_T_q[i];
		}
		int num_opp_element_boundary = 0;
		for (int j = 0; j < num_domain; j++) {
			num_opp_element_boundary += num_opp_element[j];
		}
		ofstream ofs1119("face_subdomain_unknown_T_q.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain_T_q[0][0] + num_unknown_subdomain_T_q[0][1] << ',' << num_un_total << ',' \
			<< num_unk_boundary << ',' << Vertice_Boundary_T_q[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << endl;
	}

	//Matrix_Generator_T_q(el_subdomain, oppo_element, myid);
	//Thermal_Solver_T_q(el_subdomain, oppo_element, myid);





	/*      thermal stress          */
	int* num_unknown_temp_stru = new int[num_domain * 2];
	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unknown_temp_stru[i * 2 + k] = num_unKnown_subdomain_stru[i][k];

		}
	}
	int arr_stru[2]; arr_stru[0] = 0; arr_stru[1] = 0;
	for (int i = 0; i < num_domain; i++) {
		if (i == myid) {
			arr_stru[0] = num_unKnown_subdomain_stru[i][0];
			arr_stru[1] = num_unKnown_subdomain_stru[i][1];

		}
	}

	MPI_Gather(arr_stru, 2, MPI_INT, num_unknown_temp_stru, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(num_unknown_temp_stru, num_domain * 2, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unKnown_subdomain_stru[i][k] = num_unknown_temp_stru[i * 2 + k];

		}
	}
	num_unKnown_stru = 0; //num_unKnown_Thermal = 0;
	for (int i = 0; i < num_domain; i++) {
		num_unKnown_stru += num_unKnown_subdomain_stru[i][0] + num_unKnown_subdomain_stru[i][1];
		//num_unKnown_Thermal += num_unKnown_subdomain_stru[i][0] + num_unKnown_subdomain_stru[i][1];
	}

	Accumulated_unknowns_stru = new int[num_domain];
	Accumulated_unknowns_stru[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		Accumulated_unknowns_stru[i] = Accumulated_unknowns_stru[i - 1] + num_unKnown_subdomain_stru[i - 1][0] + num_unKnown_subdomain_stru[i - 1][1];
	}
	//cout << "myid is " << myid << "num un is " << num_unKnown_subdomain_stru[myid][0] + num_unKnown_subdomain_stru[myid][1] << endl;

	if (myid == 0) {
		int num_un_total_stru = 0;
		num_un_total_stru = Accumulated_unknowns_stru[num_domain - 1] + num_unKnown_subdomain_stru[num_domain - 1][0] + num_unKnown_subdomain_stru[num_domain - 1][1];
		cout << "num_un_total_stru   " << num_un_total_stru << endl;
	}

	int temp_boundary_stru = Node_Boundary_stru[myid];
	MPI_Gather(&temp_boundary_stru, 1, MPI_INT, Node_Boundary_stru, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Node_Boundary_stru, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		int num_un_total = 0; int num_unk_boundary = 0;
		num_un_total = Accumulated_unknowns_stru[num_domain - 1] + num_unKnown_subdomain_stru[num_domain - 1][0] + num_unKnown_subdomain_stru[num_domain - 1][1];
		for (int i = 0; i < num_domain; i++) {
			num_unk_boundary += Node_Boundary_stru[i];
		}
		int num_opp_element_boundary = 0;
		for (int j = 0; j < num_domain; j++) {
			num_opp_element_boundary += num_opp_element[j];
		}
		ofstream ofs1119("face_subdomain_unknown_stru.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unKnown_subdomain_stru[0][0] + num_unKnown_subdomain_stru[0][1] << ',' << num_un_total << ',' \
			<< num_unk_boundary << ',' << Node_Boundary_stru[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << endl;
	}

	/*      thermal stress          */















	/*      E H          */



	int* num_unknown1 = (int*)malloc(num_domain * 2 * sizeof(int));
	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unknown1[i * 2 + k] = num_unknown_subdomain[i][k];

		}
	}
	int arr[2];
	for (int i = 0; i < num_domain; i++) {
		if (i == myid) {
			arr[0] = num_unknown_subdomain[i][0];
			arr[1] = num_unknown_subdomain[i][1];
		
		}
	}


	MPI_Gather(arr, 2, MPI_INT, num_unknown1, 2, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(num_unknown1, num_domain*2, MPI_INT, 0, MPI_COMM_WORLD);



	for (int i = 0; i < num_domain; i++) {
		for (int k = 0; k < 2; k++) {
			num_unknown_subdomain[i][k] = num_unknown1[i * 2 + k];

		}
	}
	num_unknown = 0;
	for (int i = 0; i < num_domain; i++) {
		num_unknown += num_unknown_subdomain[i][0] + num_unknown_subdomain[i][1];
	}

	MPI_Barrier(MPI_COMM_WORLD);

	Accumulated_unknowns = new int[num_domain];
	Accumulated_unknowns[0] = 0;
	for (int i = 1; i < num_domain; i++) {
		Accumulated_unknowns[i] = Accumulated_unknowns[i - 1] + num_unknown_subdomain[i - 1][0] + num_unknown_subdomain[i - 1][1];
	}
	cout << "myid is " << myid << "num un is " << num_unknown_subdomain[myid][0] + num_unknown_subdomain[myid][1] << endl;

	if (myid == 0) {
		int num_un_total = 0;
		num_un_total = Accumulated_unknowns[num_domain-1] + num_unknown_subdomain[num_domain-1][0] + num_unknown_subdomain[num_domain-1][1];
		cout << "num_un_total   " << num_un_total << endl;
	}

	int temp_boundary = Edge_Boundary[myid];
	MPI_Gather(&temp_boundary, 1, MPI_INT, Edge_Boundary, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Edge_Boundary, num_domain, MPI_INT, 0, MPI_COMM_WORLD);

	if (myid == 0) {
		int num_un_total = 0; int num_unk_boundary = 0;
		num_un_total = Accumulated_unknowns[num_domain - 1] + num_unknown_subdomain[num_domain - 1][0] + num_unknown_subdomain[num_domain - 1][1];
		for (int i = 0; i < num_domain; i++) {
			num_unk_boundary += Edge_Boundary[i];
		}
		int num_opp_element_boundary = 0;
		for (int j = 0; j < num_domain; j++) {
			num_opp_element_boundary += num_opp_element[j];
		}
		ofstream ofs1119("face_subdomain_unknown.csv", ios::app);
		ofs1119 << num_domain << ',' << num_unknown_subdomain[0][0] + num_unknown_subdomain[0][1] << ',' << num_un_total << ',' \
			<< num_unk_boundary << ',' << Edge_Boundary[myid] << ','\
			<< process_threads << ',' << num_opp_element_boundary << endl;

	}

	/*       E H            */

	//if (method < 0) {
	//	Eqn_Solver_Iteration(el_subdomain, oppo_element, myid);
	//}
	//else {
	//	Eqn_Solver(el_subdomain, oppo_element, myid);
	//}
	if (EM_solve) {
		Eqn_Solver(el_subdomain, oppo_element, myid);
		Solve_involve_Freq(el_subdomain, myid);
	}





	//delete[]p_pardiso; delete[]ip; delete[]jp; delete[]acoo1; delete[]acoo; delete[]acoo_ome; delete[]acoo_ome_2; delete[]rowind; delete[]colind; delete[]num_nzero_Pmatrix;
	//delete[]m; delete[]m1; delete[]m_ome; delete[]m_ome_2;    delete[]nnz_c; delete[]mrow; delete[]mcol; delete[]r_dm; delete[]fbb;
	//delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm; delete[]r_dm;


	if (EM_solve&&E_T) {
		Matrix_Generator_T_q(el_subdomain, oppo_element, myid);
		Thermal_Solver_T_q(el_subdomain, oppo_element, myid);
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
		for (int j = 0; j < num_element_subdomain; j++) {
			for (int ii = 0; ii < 4; ii++) {
				if (el_subdomain[j].node[ii].unknown_T != 0) {
					node_global_temp[el_subdomain[j].node[ii].unknown_T - 1] = el_subdomain[j].node[ii].ver;
				}
			}
		}
		for (int j = 0; j < num_unknown_subdomain_T_q[myid][0]; j++) {
			X_h_temp[j] = Xh[j];
		}
		MPI_Gatherv(X_h_temp, num_unknown_subdomain_T_q[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Gatherv(node_global_temp, num_unknown_subdomain_T_q[myid][0], MPI_INT, V_T, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
		//delete[]Accumulated_T_temp; delete[]unknown_T_temp;  delete[]node_global_temp; delete[]X_h_temp; //delete[]Xh;
		double* XX1 = new double[num_node];
		double* XX2 = new double[num_node];
		for (int i = 0; i < num_node; i++) {
			XX1[i] = T0;
			XX2[i] = T0;
		}
		if (myid == 0) {
			for (int i = 0; i < un_T_temp; i++) {
				int temp_ver = V_T[i];
				XX1[temp_ver - 1] = X_T[i];
			}
		}

		double error_sigma = 1;
		double error_T = 1;
		int iteration_number_E_T = 1;
		double sigma0 = 1 / (1.75e-8);
		double sigma2 = 0.004;
		double epsil0 = 9.0;
		double sigma1;
		//double alepa_T = 0.0039;
		//double alepa_T = 0.0025;
		int num_sigma = 0;
		for (int el = 0; el < num_element_subdomain; el++) {
			int material = el_subdomain[el].Material;
			if (material == 2 || material == 3) {
				//sigma0 = el_subdomain[el].sigma;
				num_sigma++;
			}
		}
		double* sigma_v1 = new double[num_sigma];
		double* sigma_v2 = new double[num_sigma];
		for (int i = 0; i < num_sigma; i++) {
			sigma_v1[i] = sigma0;
			sigma_v2[i] = sigma0;
		}
		double* error_sigma_v = new double[num_domain];
		vector<double>error_sigma_iteration;
		vector<double>error_T_iteration;
		cout << "TTTTTTTT" << endl;
		Freq = 0;


		for (int nfreq = 0; Freq < fre_end; ++nfreq) {
			Freq = ((double)(nfreq)) * fre_step + fre_start;
			omega = 2.0 * pi * Freq;
			//Solve_E_H_boundary(el_subdomain, myid);
			error_T = 1;
			iteration_number_E_T = 0;
			while (error_T > 1e-10 && iteration_number_E_T < 20) {
				iteration_number_E_T++;
				int count_sigma = 0;
				for (int el = 0; el < num_element_subdomain; el++) {
					for (int nn = 0; nn < 4; nn++) {
						int op = el_subdomain[el].face[nn].opp[0]; int ofn = el_subdomain[el].face[nn].opp[1];
					}
					double sigma__temp = el_subdomain[el].sigma;
					double Material__temp = el_subdomain[el].Material;
					if (sigma__temp != 0) { //Material__temp==3
						double T_sum = 0;
						for (int node_ii = 0; node_ii < 4; node_ii++) {
							int unknown_TT = el_subdomain[el].node[node_ii].unknown_T;
							if (unknown_TT == 0) {
								T_sum += T0;
							}
							else {
								T_sum += Xh[unknown_TT - 1];
							}
						}
						double T_average = T_sum / 4;
						sigma2 = el_subdomain[el].sigma_temp;
						sigma1 = sigma2 / (1 + alepa_T * (T_average - T0));
						el_subdomain[el].sigma = sigma1;
					}
				}
				Eqn_Solver(el_subdomain, oppo_element, myid);
				Solve_E_H_boundary(el_subdomain, myid);
				Matrix_Generator_T_q(el_subdomain, oppo_element, myid);
				Thermal_Solver_T_q(el_subdomain, oppo_element, myid);
				for (int j = 0; j < num_unknown_subdomain_T_q[myid][0]; j++) {
					X_h_temp[j] = Xh[j];
				}
				MPI_Gatherv(X_h_temp, num_unknown_subdomain_T_q[myid][0], MPI_DOUBLE, X_T, unknown_T_temp, Accumulated_T_temp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				//MPI_Gatherv(node_global_temp, num_unknown_subdomain_T_q[myid][0], MPI_INT, V_T, unknown_T_temp, Accumulated_T_temp, MPI_INT, 0, MPI_COMM_WORLD);
				if (myid == 0) {
					for (int i = 0; i < un_T_temp; i++) {
						int temp_ver = V_T[i];
						XX2[temp_ver - 1] = X_T[i];
					}
					double error_T2_fz = 0;
					double error_T2_fm = 0;
					for (int i = 0; i < num_node; i++) {
						error_T2_fz += (XX2[i] - XX1[i]) * (XX2[i] - XX1[i]);
						error_T2_fm += XX1[i] * XX1[i];
						XX1[i] = XX2[i];
					}
					error_T = error_T2_fz / error_T2_fm;
				}
				MPI_Bcast(&error_T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				error_T_iteration.push_back(error_T);
				if (myid == 0) {
					cout << "error_T = " << error_T << endl;
					cout << "iteration_number_E_T = " << iteration_number_E_T << endl;
				}
			}
			Out_field_subdomain_new(el_subdomain, nfreq, myid, 1);

			cout << "myid is " << myid << " iteration_number_E_T = " << iteration_number_E_T << endl;
			if (myid == 0) {
				ofstream ofs11190("error_T_iteration.csv", ios::app);
				for (int i = 0; i < error_T_iteration.size(); i++) {
					ofs11190 << error_T_iteration[i] << endl;
				}
				ofstream ofs11191("error_T_vf.csv", ios::app);
				ofs11191 << (double)(Freq / 1e9) << ',' << error_T << endl;
				ofstream ofs11192("error_T_iteration_number.csv", ios::app);
				ofs11192 << (double)(Freq / 1e9) << ',' << iteration_number_E_T << endl;
			}


		}
		Matrix_Generator_stru(el_subdomain, oppo_element, myid);
	}
	else if(!TF_only){
		cout << "dnhbdhbcfjdnn" << endl;
		Freq = fre_E_T_F;
		omega = 2.0 * pi * Freq;
		Solve_E_H_boundary2(el_subdomain, myid);
		Matrix_Generator_T_q(el_subdomain, oppo_element, myid);
		Thermal_Solver_T_q(el_subdomain, oppo_element, myid);
		Matrix_Generator_stru(el_subdomain, oppo_element, myid);
	}

	if (TF_only) {
		Matrix_Generator_T_q(el_subdomain, oppo_element, myid);
		Thermal_Solver_T_q2(el_subdomain, oppo_element, myid);
		Matrix_Generator_stru(el_subdomain, oppo_element, myid);
	}



	end = time(NULL);

	double time11 = (double)(end - start);
	if (myid == 0) {
		cout << "time is " << time11 << endl;
	}
	
	MPI_Finalize();
	return 0;
}






