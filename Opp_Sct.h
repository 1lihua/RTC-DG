#pragma once
#include "Element.h"
int Opp_Sct( Element*& el,int num_element_subdomain,int myid);
int Opp_Sct1(Element*& el, int num_element_boundary_total, int myid, int* Global_num_element_boundary_total, int* Local_num_face_boundary_total, int* address_offset_face);