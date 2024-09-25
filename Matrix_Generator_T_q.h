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
int Matrix_Generator_T_q(Element* element, Element* oppo_element, int myid);
int Thermal_Solver_T_q(Element* element, Element* oppo_element, int myid);
int Thermal_Solver_T_q2(Element* element, Element* oppo_element, int myid);

int Matrix_Generator_T_q_2(Element* element, Element* oppo_element, int myid);
int Thermal_Solver_T_q_2(Element* element, Element* oppo_element, int myid);