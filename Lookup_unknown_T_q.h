#pragma once
#include <iostream>
#include <Eigen>
#include"Element.h"
#include"Global_Data.h"

int unknownIndex_Thermal(Element* el, int myid);


int unknownIndex_Thermal_2(Element* element, int myid);

int unknownIndex_stru(Element* element, int myid);