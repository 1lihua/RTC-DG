#include"General_Parameter.h"



double Det_Cal(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33) {

    double Det_Cal_o = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a13 * a22 * a31 - a11 * a23 * a32 - a12 * a21 * a33;
    return Det_Cal_o;
}
int General_Parameter(Element* el,int num_element_subdomain) {   //revise

    
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, s;
    double a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44;
    double vtemp[3];
    double vol_temp[4][4];
    double temp_a[3][3];
    double temp1[2][2], temp2[2][2], temp3[2][2];
  

    for (size_t i = 0; i != num_element_subdomain; ++i) {
        x1 = el[i].node[0].zb[zbX]; y1 = el[i].node[0].zb[zbY]; z1 = el[i].node[0].zb[zbZ];
        x2 = el[i].node[1].zb[zbX]; y2 = el[i].node[1].zb[zbY]; z2 = el[i].node[1].zb[zbZ];
        x3 = el[i].node[2].zb[zbX]; y3 = el[i].node[2].zb[zbY]; z3 = el[i].node[2].zb[zbZ];
        x4 = el[i].node[3].zb[zbX]; y4 = el[i].node[3].zb[zbY]; z4 = el[i].node[3].zb[zbZ];

        a11 = 1; a12 = x1; a13 = y1; a14 = z1;
        a21 = 1; a22 = x2; a23 = y2; a24 = z2;
        a31 = 1; a32 = x3; a33 = y3; a34 = z3;
        a41 = 1; a42 = x4; a43 = y4; a44 = z4;
        //Calculate the volume of tetrahedron
        el[i].Ve = (double)abs((Det_Cal(a22, a23, a24, a32, a33, a34, a42, a43, a44) - Det_Cal(a12, a13, a14, a32, a33, a34, a42, a43, a44) + Det_Cal(a12, a13, a14, a22, a23, a24, a42, a43, a44) - Det_Cal(a12, a13, a14, a22, a23, a24, a32, a33, a34))) / 6;

        // calculate ai
        a11 = x2; a12 = y2; a13 = z2; a21 = x3; a22 = y3; a23 = z3; a31 = x4; a32 = y4; a33 = z4;
        el[i].a[0] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = z1; a21 = x3; a22 = y3; a23 = z3; a31 = x4; a32 = y4; a33 = z4;
        el[i].a[1] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = z1; a21 = x2; a22 = y2; a23 = z2; a31 = x4; a32 = y4; a33 = z4;
        el[i].a[2] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = z1; a21 = x2; a22 = y2; a23 = z2; a31 = x3; a32 = y3; a33 = z3;
        el[i].a[3] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        // calculate bi
        a11 = 1; a12 = y2; a13 = z2; a21 = 1; a22 = y3; a23 = z3; a31 = 1; a32 = y4; a33 = z4;
        el[i].b[0] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = 1; a12 = y1; a13 = z1; a21 = 1; a22 = y3; a23 = z3; a31 = 1; a32 = y4; a33 = z4;
        el[i].b[1] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = 1; a12 = y1; a13 = z1; a21 = 1; a22 = y2; a23 = z2; a31 = 1; a32 = y4; a33 = z4;
        el[i].b[2] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = 1; a12 = y1; a13 = z1; a21 = 1; a22 = y2; a23 = z2; a31 = 1; a32 = y3; a33 = z3;
        el[i].b[3] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        // calculate ci
        a11 = x2; a12 = 1; a13 = z2; a21 = x3; a22 = 1; a23 = z3; a31 = x4; a32 = 1; a33 = z4;
        el[i].c[0] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = 1; a13 = z1; a21 = x3; a22 = 1; a23 = z3; a31 = x4; a32 = 1; a33 = z4;
        el[i].c[1] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = 1; a13 = z1; a21 = x2; a22 = 1; a23 = z2; a31 = x4; a32 = 1; a33 = z4;
        el[i].c[2] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = 1; a13 = z1; a21 = x2; a22 = 1; a23 = z2; a31 = x3; a32 = 1; a33 = z3;
        el[i].c[3] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        // calculate di
        a11 = x2; a12 = y2; a13 = 1.0; a21 = x3; a22 = y3; a23 = 1.0; a31 = x4; a32 = y4; a33 = 1.0;
        el[i].d[0] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = 1.0; a21 = x3; a22 = y3; a23 = 1.0; a31 = x4; a32 = y4; a33 = 1.0;
        el[i].d[1] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = 1.0; a21 = x2; a22 = y2; a23 = 1.0; a31 = x4; a32 = y4; a33 = 1.0;
        el[i].d[2] = -Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        a11 = x1; a12 = y1; a13 = 1.0; a21 = x2; a22 = y2; a23 = 1.0; a31 = x3; a32 = y3; a33 = 1.0;
        el[i].d[3] = Det_Cal(a11, a12, a13, a21, a22, a23, a31, a32, a33);

        // calculate length of six edges


        for (int ii = 0; ii < 12; ii++) {
            for (int j = 0; j < 3; j++) {
                vtemp[j] = el[i].node[edge_node_local[ii][0] - 1].zb[j] - el[i].node[edge_node_local[ii][1] - 1].zb[j];

            }
            el[i].length[ii] = sqrt(vtemp[0] * vtemp[0] + vtemp[1] * vtemp[1] + vtemp[2] * vtemp[2]);
        }

        // calculate area of four faces

        s = (double)(el[i].length[0] + el[i].length[1] + el[i].length[3]) / 2.0;
        el[i].face[0].Area = (double)sqrt(s * (s - el[i].length[0]) * (s - el[i].length[1]) * (s - el[i].length[3]));

        s = (double)(el[i].length[1] + el[i].length[2] + el[i].length[5]) / 2.0;
        el[i].face[1].Area = (double)sqrt(s * (s - el[i].length[1]) * (s - el[i].length[2]) * (s - el[i].length[5]));

        s = (double)(el[i].length[0] + el[i].length[2] + el[i].length[4]) / 2.0;
        el[i].face[2].Area = (double)sqrt(s * (s - el[i].length[0]) * (s - el[i].length[2]) * (s - el[i].length[4]));

        s = (double)(el[i].length[3] + el[i].length[4] + el[i].length[5]) / 2.0;
        el[i].face[3].Area = (double)sqrt(s * (s - el[i].length[3]) * (s - el[i].length[4]) * (s - el[i].length[5]));

    }
    return 0;
}












