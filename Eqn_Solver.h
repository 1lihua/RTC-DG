#pragma once
#include"Element.h"
int Eqn_Solver(Element* el, Element* oppo_element,int myid);
int Eqn_Solver_Iteration(Element* el, Element* oppo_element, int myid);

int Solve_E_H_boundary_Iteration(Element* el, int myid);
int Solve_E_H_boundary(Element* el, int myid);
int Solve_E_H_boundary2(Element* el, int myid);
int Solve_E_H_boundary_method2(Element* el, int myid);
//void Matrix_Fullfill(Element* el);

int Solve_E_H_boundary_scalapack(Element* el, int myid);


int Solve_E_H_boundary_method3(Element* el, int myid);

int Solve_E_H_boundary_iter(Element* el, int myid);
int Solve_E_H_boundary_iter_long(Element* el, int myid);



int Solve_E_H_boundary_E_T(Element* el, int myid);
int Eqn_Solver_T_E(Element* el, Element* oppo_element, int myid);