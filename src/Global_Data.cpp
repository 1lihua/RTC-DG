#include"mkl.h"
#include<cmath>
#include<complex>
#include<vector>
#include<Eigen>
using namespace std;
//General Parameter
extern const double sigma_IBC = 1 / (1.75e-8);
//double Waveguide_P = 1.1 / 1.4468 * 2.85 / 3.19 * 2.85 / 3.19 / 10 / 36.85 * 14.85*1.2*10.4*10.4/7.15/7.15/10*4*7.86*7.86/9/9*5.84*5.84/5.94/5.94/9.57/9.57*10.2*10.2*4.51/4.58*4.51/4.58;//10w     3.174164462785

double SI_meter = 1e-3;//4w
bool EM_solve = true;
bool TF_only = false;
double powerful=1.0;
double fre_start=1e9;
double fre_step = 0.01e9;

double fre_end=1e9;

double fre_pattern=1e15;

double fre_E_T_F = 3e9;

bool LUMP_O = false;
bool E_T = false;
double alepa_T = 0.0;



double* H_lump = nullptr;
double* W_lump = nullptr;

bool WAVE_O = false;;
int wave_material = 1;

double wave_zzz=0.0;

int num_wave=0;
double rb_wave;
double ra_wave;
double* x_wave = nullptr;
double* y_wave = nullptr;
double* z_wave = nullptr;

double* phase_wave = nullptr;
double* Amplitude_wave = nullptr;


double* DIEeps = nullptr;
double* DIEmu = nullptr;
double* DIEsigma = nullptr;
double* DIEk1 = nullptr;
double* DIEk2 = nullptr;
double* DIEk3 = nullptr;
double* DIERho = nullptr;
double* DIESHC = nullptr;
double* DIEQ = nullptr;
double* DIEE = nullptr;
double* DIEnu = nullptr;
double* DIEalpha = nullptr;


double T0 = 293.15;
double Td = 293.15;
double Hcvt = 10.0;
double Ta = 293.15;







double Waveguide_P = 0.020169206198702/6.35/5.645;//4w


double weight0_2D[6] = { {0.223381589678011},{0.223381589678011},{0.223381589678011},{0.109951743655322},{0.109951743655322},{0.109951743655322} };
double alpha0_2D[6] = { {0.108103018168070},{0.445948490915965},{0.445948490915965},{0.816847572980459},{0.091576213509771},{0.091576213509771} };
double beta0_2D[6] = { {0.445948490915965},{0.108103018168070},{0.445948490915965},{0.091576213509771},{0.816847572980459},{0.091576213509771} };
double gamma0_2D[6] = { {0.445948490915965},{0.445948490915965},{0.108103018168070},{0.091576213509771},{0.091576213509771},{0.816847572980459} };



double weight03_3D[5] = { {-0.8},{0.45},{0.45},{0.45},{0.45} };
double alpha0_3D[5] = { {0.25},{0.50000000000000},{0.166666666666667},{0.166666666666667},{0.166666666666667} };
double beta0_3D[5] = { {0.25},{0.166666666666667},{0.50000000000000},{0.166666666666667},{0.166666666666667} };
double gamma0_3D[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.50000000000000},{0.166666666666667} };
double delta0_3D[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.166666666666667},{0.50000000000000} };







int num_domain = 0;
int num_node = 0;
int num_element = 0;
int num_unKnown_b = 0;
int num_unknown = 0;
int num_edge = 0;
int* Edge_Boundary = nullptr;
int** num_unknown_subdomain = nullptr;
int *mesh_at_RTC = nullptr;
int num_face_RTC = 0;
int* num_nzero_Pmatrix = 0;
//Solution Module
int* row_dm = nullptr;
int* nzero_dm = nullptr;
int* nn2_dm = nullptr;
int* nn3_dm = nullptr;
int* nn4_dm = nullptr;
MKL_Complex16 *unKnownUi = nullptr;
MKL_Complex16* unKnownUi_temp = nullptr;
MKL_Complex16* unKnownUi2 = nullptr;//added for face cycle
//size_t** nnz_c = nullptr;
size_t** nnz_c = nullptr;
size_t** nnz_c1 = nullptr;//
int* mrow = nullptr;
int* mcol = nullptr;
MKL_Complex16* m = nullptr;
MKL_Complex16* fbr = nullptr;


//MKL_Complex16* uk = nullptr;
//
//MKL_Complex16* uk1 = nullptr;

MKL_Complex16* acoo = nullptr;
int* rowind = nullptr;
int* colind = nullptr;



int* rowindC = nullptr;
int* colindC = nullptr;


int* rowindA = nullptr;
int* colindA = nullptr;

int* row_C = nullptr;//
int* col_C = nullptr;//

MKL_Complex16* C = nullptr;//added for face cycle
MKL_Complex16* m1 = nullptr;
MKL_Complex16* acoo1 = nullptr;
MKL_Complex16* m_ome;
MKL_Complex16* m_ome_2;
MKL_Complex16* m_wave;
MKL_Complex16* m_ome_sqrt;


MKL_Complex16* acoo_ome;
MKL_Complex16* acoo_ome_2;
MKL_Complex16* acoo_wave;
MKL_Complex16* acoo_ome_sqrt;

//MKL_Complex16** p_pardiso = nullptr;
//int** ip = nullptr;
//int** jp = nullptr;
//vector<MKL_Complex16*>p_pardiso;
//vector<int*>ip;
//vector<int*>jp;
MKL_Complex16* p_pardiso;
int* ip, *jp;
std::complex<double>* Xsol  = nullptr;

//! free space light speed and pi
extern const double pi = 3.14159265358979;
extern const double VC = 299792458;

//!free space permittivity, permeabilityand intrinsic wave impedance
extern const double epsilon0 = 8.85418781762039e-12;
extern const double mur = 1.25663706143592e-6;
extern const double eta0 = sqrt(mur / epsilon0);

//!local node number in each face
extern const MKL_INT face_node[4][3] = { {1,2,3},{ 1, 3, 4 },{ 1, 2, 4 },{ 2, 3, 4 } };

//! edge number in each face
extern const MKL_INT face_edge[4][6] = { {1,2,4,7,8,10}, {2,3,6,8,9,12},{1,3,5,7,9,11},{4,5,6,10,11,12}};

//!local node number for each edge
extern const MKL_INT  edge_node_local[12][2] = { {1, 2} , {1, 3} , {1, 4} , {2, 3}, {4, 2} , {3, 4} , {1, 2}, {1, 3}, {1, 4} , {2, 3} , {4, 2}, {3, 4} };

//!EiHi_index including the faces for edge ii and the edge number in that face
extern const MKL_INT edge_by_face[12][2] = { {1, 3} , {1, 2} , {2, 3} , {1, 4} , {3, 4} , {2, 4} , {1, 3} , {1, 2} , {2, 3} , {1, 4} , {3, 4} , {2, 4} };



// Have changed parameters as list: 
double *epsil_r = nullptr;
double Freq, omega;
//double Hsubstrate = 3.5e-3;
//double Wmicrostrip = 20e-3;
//double Hsubstrate1 = 3.5e-3;
//double Wmicrostrip1 = 20e-3;

//double Hsubstrate = 1.524E-3;
//double Wmicrostrip = 3.2E-3;

double Hsubstrate = 1.524e-3;
double Wmicrostrip = 3.2e-3;

//double Hsubstrate = 1.6e-2;
//double Wmicrostrip = 2.5e-2;

//double Hsubstrate = 1e-3;
//double Wmicrostrip = 2.76e-3;

double factor = 10;
double Z0 = 50;

//extern const double Freq = 10.0e9;
//extern const double omega = 2.0 * pi * Freq;
//extern const double Kwave[3] = { 0.0, omega * sqrt(epsilon0 * mur), 0.0 };

//!Factors to alter the condition number of the matrix system
extern const double FactorH = 1.0 / eta0;
extern const double FactorE = eta0;
double FactorOmega;

int kpmc1_GB = 0;
int kpmc2_GB = 0;

int port_num = 0;


vector<complex<double>>V11, V21, V12, V22;



//BoundaryInformation
int* boundaryInfo = nullptr;
int** Bondver = nullptr;
Eigen::MatrixXd zb_boundary;
int num_element_boundary = 0;
int num_node_boundary = 0;


Eigen::MatrixXd zb;
Eigen::MatrixXi ver;
complex<double> S11;
complex<double> S21;




int num_element_subdomain = 0;


extern int* num_opp_element = nullptr;

int* Accumulated_unknowns = nullptr;

double xrod1;
double yrod1;
//const double ra_wave, rb_wave, height;


int process_threads;
int num_face_boundary;





int* Global_num_element_boundary = nullptr;
int* Local_num_face_boundary = nullptr;





Eigen::VectorXd fbr_T_q;
Eigen::VectorXd fbt;






//double* Me_full;
int* iRow, * jCol;
//double* M_full;
int* mRow, * mCol;
double* fbb_T_q_G;

double* acoo_T_q;
//int* rowind, * colind;
//double* acoo1;
int* rowind1, * colind1;
double* acoo1_T_q;
double* acoo_1;
double* acoo_2;
double* m_T_q;
//double* m;
double* m_1;
double* m_2;
int* mrow_T_q, * mcol_T_q;
Eigen::MatrixXi nnz_C;

//int* mrowind, * mcolind;

double* m_pardiso;
int* im, * jm;

double* a_pardiso;
int* ia, * ja;
double* at_pardiso;
int* iat, * jat;

double* p_pardiso_T_q;
double* pt_pardiso;
int* ipt, * jpt;

//int* num_nzero_Pmatrix;

//double* Xsol;
double* Xs0;
double* Xsol_H;
double* Xs0_H;


int* nn3_dm_T_q;


int nnzzz;


double* Xe, * Xh, * Xhold;
double* xx;
double* Xphi, * Xtemp, * Xtemp1;
double* Zj, * Vj;
Eigen::MatrixXd Hn, Qn, Qzn;

int num_convergence;



double time_step;
double time_total;



int num_unKnown_eb = 0;

double* fbr_T_q1 = nullptr;


int num_unknown_T_q;
int num_unKnown_Thermal;

int* Vertice_Boundary_T_q = nullptr;
int** num_unknown_subdomain_T_q = nullptr;

int* Accumulated_unknowns_T_q;

int* ip_T_q, * jp_T_q;

int* num_nzero_Pmatrix_T_q;

int* r_dm_T_q;


//
//double weight0[5] = { {-0.8},{0.45},{0.45},{0.45},{0.45} };
//double alpha0[5] = { {0.25},{0.50000000000000},{0.166666666666667},{0.166666666666667},{0.166666666666667} };
//double beta0[5] = { {0.25},{0.166666666666667},{0.50000000000000},{0.166666666666667},{0.166666666666667} };
//double gamma0[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.50000000000000},{0.166666666666667} };
//double delta0[5] = { {0.25},{0.166666666666667},{0.166666666666667},{0.166666666666667},{0.50000000000000} };






/*  stru   */

int num_unKnown_stb;
int num_unKnown_stru;
int* Node_Boundary_stru;

int** num_unKnown_subdomain_stru = nullptr;

//int num_unknown_stru;


int* Vertice_Boundary_stru = nullptr;


int* Accumulated_unknowns_stru;

int* ip_stru, * jp_stru;
double* p_pardiso_stru;

int* num_nzero_Pmatrix_stru;

int* r_dm_stru;

double* fbr_stru = nullptr;

double* nu;
double* mu;
double* lambda;
double* alpha;
int* iRow_stru, * jCol_stru;


double* acoo_stru;
int* rowind_stru, * colind_stru;
double* acoo1_stru;
int* rowind1_stru, * colind1_stru;
double* acoo_1_stru;
double* acoo_2_stru;
double* m_stru;
double* m_1_stru;
double* m_2_stru;
int* mrow_stru, * mcol_stru;

double* fbb_stru;

int** nnz_C_stru;
double* xx_stru;
double* Xstru;