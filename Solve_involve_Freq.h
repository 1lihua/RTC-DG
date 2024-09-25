#pragma once
#include"Eqn_Solver.h"
#include"OutField.h"
#include"Global_Data.h"
#include<Eigen>
#include <fstream>
#include<mpi.h>
#include<omp.h>
#define method  1   //1 or 2
void   Solve_involve_Freq(Element* el_subdomain, int myid);