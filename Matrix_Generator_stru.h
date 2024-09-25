#pragma once
#include <iostream>
#include <Eigen>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include "mkl.h"
#include "Global_Data.h"
#include "Matrix_Module.h"
#include "Element.h"
#include<mpi.h>
#include<omp.h>
int Matrix_Generator_stru(Element* element, Element* oppo_element, int myid);
//int Thermal_Solver_T_q(Element* element, Element* oppo_element, int myid);

