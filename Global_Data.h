#pragma once

#include"mkl.h"
#include<cmath>
#include<complex>
#include <vector>
#include<Eigen>
#define el_adj 0
#define face_suerposition 1
using namespace std;

extern double SI_meter;

extern double fre_start;
extern double fre_step;
extern double fre_end;


extern double powerful;


extern bool EM_solve;
extern bool TF_only;

extern double fre_pattern;

extern double fre_E_T_F;

extern bool E_T;

extern double alepa_T;


extern bool LUMP_O;


extern double* H_lump;
extern double* W_lump;


extern bool WAVE_O;

extern int wave_material;

extern double wave_zzz;

extern int num_wave;


extern double rb_wave;
extern double ra_wave;
extern double* x_wave;
extern double* y_wave;
extern double* z_wave;

extern double* phase_wave;
extern double* Amplitude_wave;

extern double* DIEeps;
extern double* DIEmu;
extern double* DIEsigma;
extern double* DIEk1;
extern double* DIEk2;
extern double* DIEk3;
extern double* DIERho;
extern double* DIESHC;
extern double* DIEQ;
extern double* DIEE;
extern double* DIEnu;
extern double* DIEalpha;


extern double T0;
extern double Td;
extern double Hcvt;

extern double Ta;










extern const double sigma_IBC;
extern double Waveguide_P;

extern double weight0_2D[6];
extern double alpha0_2D[6];
extern double beta0_2D[6];
extern double gamma0_2D[6];

extern double weight03_3D[5];
extern double alpha0_3D[5];
extern double beta0_3D[5];
extern double gamma0_3D[5];
extern double delta0_3D[5];

extern int num_domain;
extern int num_node;
extern int num_element;
extern int num_unKnown_b;
extern int num_edge;
extern int* Edge_Boundary;
extern int** num_unknown_subdomain;
extern int* mesh_at_RTC;
extern int num_face_RTC;
extern int num_unknown;
extern int* num_nzero_Pmatrix;
extern int* row_dm;
extern int* nzero_dm;
extern int* nn2_dm;
extern int* nn3_dm;
extern int* nn4_dm;
extern MKL_Complex16* unKnownUi;
extern MKL_Complex16* unKnownUi_temp;
extern MKL_Complex16* unKnownUi2;//added for face cycle
//extern size_t** nnz_c;
extern size_t** nnz_c;
extern size_t** nnz_c1;
extern int* mrow;
extern int* mcol;
extern  MKL_Complex16* m;
extern MKL_Complex16* fbr;

extern MKL_Complex16* acoo;
extern  MKL_Complex16* C;
extern int* rowind;
extern int* colind;

extern int* rowindA;
extern int* colindA;

extern int* rowindC;
extern int* colindC;

extern int* row_C;
extern int* col_C;

//
//extern MKL_Complex16** p_pardiso;
//extern int** ip ;
//extern int** jp ;
//extern std::vector<MKL_Complex16*>p_pardiso;
//extern std::vector<int*>ip;
//extern std::vector<int*>jp;
extern std::complex<double>* Xsol;
extern MKL_Complex16* p_pardiso;
extern int* ip, *jp;
extern MKL_Complex16* m_ome;
extern MKL_Complex16* m_ome_2;
extern MKL_Complex16* m_wave;
extern MKL_Complex16* m_ome_sqrt;

extern MKL_Complex16* acoo_ome;
extern MKL_Complex16* acoo_ome_2;
extern MKL_Complex16* acoo_wave;
extern MKL_Complex16* acoo_ome_sqrt;

extern MKL_Complex16* m1 ;
extern MKL_Complex16* acoo1;
//! free space light speed and pi
extern const double pi;
extern const double VC;

//!free space permittivity, permeabilityand intrinsic wave impedance
extern const double epsilon0;
extern const double mur;
extern const double eta0;

extern double Hsubstrate;
extern double Wmicrostrip;
//extern double Hsubstrate1;
//extern double Wmicrostrip1;

//!local node number in each face
extern const MKL_INT face_node[4][3];

//! edge number in each face
extern const MKL_INT face_edge[4][6];

//!local node number for each edge
extern const MKL_INT  edge_node_local[12][2];

//!EiHi_index including the faces for edge ii and the edge number in that face
extern const MKL_INT edge_by_face[12][2];




extern  double *epsil_r;
extern  double Freq;
extern  double omega;

//extern const double Kwave[3];

//!Factors to alter the condition number of the matrix system
extern const double FactorH;
extern const double FactorE;
extern double FactorOmega;

extern int kpmc1_GB ;
extern int kpmc2_GB ;
extern int port_num;
extern std::vector<std::complex<double>> V11, V21, V12, V22;



//BoundaryInformation

extern int** Bondver;
extern Eigen::MatrixXd zb_boundary;
extern int num_element_boundary;
extern int num_node_boundary;
extern int** Bondver;
extern int* boundaryInfo;


/***parameters for the structure***/
//extern const double ra_wave , rb_wave , height;
//extern const double xrod1, yrod1;
extern double xrod1, yrod1;
extern double factor;
extern double Z0;

extern Eigen::MatrixXd zb;
extern Eigen::MatrixXi ver;
extern complex<double> S11;
extern complex<double> S21;



extern int num_element_subdomain;

extern int* num_opp_element;  //number of opposite element
extern int* Accumulated_unknowns;

//const double ra_wave = 0.6e-3;
//const double rb_wave = 1.5e-3;

extern int process_threads;

extern int num_face_boundary;



const double h_wg = 9.52500e-3;
const double w_wg = 19.050e-3;
//const double y_o = -w_wg / 2;
const double y_o = -w_wg / 2;
const double x_o = -2.5e-2;



extern int* Global_num_element_boundary;
extern int* Local_num_face_boundary;



extern Eigen::VectorXd fbr_T_q;
extern Eigen::VectorXd fbt;








//extern double* Me_full;
extern int* iRow, * jCol;
//extern double* M_full;
//extern int* mrowind, * mcolind;
extern double* m_pardiso;
extern int* im, * jm;
extern double* fbb_T_q_G;

extern double* acoo_T_q;
//extern int* rowind, * colind;
extern double* acoo1_T_q;
extern int* rowind1, * colind1;

extern double* acoo_T_q;
extern double* acoo_1;//useless
extern double* acoo_2;//useless
extern double* m_T_q;
extern double* m_1;//useless
extern double* m_2;//useless
extern int* mrow_T_q, * mcol_T_q;
extern Eigen::MatrixXi nnz_C;

extern double* a_pardiso;
extern int* ia, * ja;
extern double* at_pardiso;
extern int* iat, * jat;

extern double* p_pardiso_T_q;
//extern int* ip, * jp;
extern double* pt_pardiso;
extern int* ipt, * jpt;

extern int* num_nzero_Pmatrix;

extern double* Xsol_T_q;
extern double* Xs0;
extern double* Xsol_H;
extern double* Xs0_H;


extern int* nn3_dm_T_q;


extern int nnzzz;







extern double* Xe, * Xh, * Xhold;
extern double* xx;
extern double* Xphi, * Xtemp, * Xtemp1;
extern double* Zj, * Vj;
extern Eigen::MatrixXd Hn, Qn, Qzn;
extern int num_convergence;



extern double time_step;

extern double time_total;



extern int num_unKnown_eb;


extern double* fbr_T_q1;


extern int* Vertice_Boundary_T_q;

extern int** num_unknown_subdomain_T_q;

extern int num_unknown_T_q;
extern int num_unKnown_Thermal;

extern int* Accumulated_unknowns_T_q;


extern int* ip_T_q, * jp_T_q;

extern int* num_nzero_Pmatrix_T_q;

extern int* r_dm_T_q;

extern double weight0[5];
extern double alpha0[5] ;
extern double beta0[5] ;
extern double gamma0[5];
extern double delta0[5];





/*  stru   */


extern int num_unKnown_stb;
extern int num_unKnown_stru;
extern int* Node_Boundary_stru;

extern int** num_unKnown_subdomain_stru;

//extern int num_unknown_stru;


extern int* Vertice_Boundary_stru ;


extern int* Accumulated_unknowns_stru;

extern int* ip_stru, * jp_stru;
extern double* p_pardiso_stru;

extern int* num_nzero_Pmatrix_stru;

extern int* r_dm_stru;

extern double* fbr_stru;

extern double* nu;
extern double* mu;
extern double* lambda;
extern double* alpha;
extern int* iRow_stru, * jCol_stru;

extern double* acoo_stru;
extern int* rowind_stru, * colind_stru;
extern double* acoo1_stru;
extern int* rowind1_stru, * colind1_stru;
extern double* acoo_1_stru;
extern double* acoo_2_stru;
extern double* m_stru;
extern double* m_1_stru;
extern double* m_2_stru;
extern int* mrow_stru, * mcol_stru;

extern double* fbb_stru;
extern int** nnz_C_stru;

extern double* xx_stru;

extern double* Xstru;