#pragma once
#include"Element.h"
int Out_field(Element* el, const int& nfreq);
int Out_field_subdomain(Element* el, const int& nfreq);
int Out_field_sub(Element* el, const int& nfreq);
int Out_field_subdomain_double(Element* el, const int& nfreq);

int Out_field_subdomain_new(Element* el, const int& nfreq, int myid, int ET);