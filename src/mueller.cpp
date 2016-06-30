/********* mueller.cpp ***********
 * Author: Cody Bezik
 * Last modified: 6-27-2016
 * Project: swarm.cpp
 * 
 * Description:
 * This file contains the definitions of functions and parameters necessary to describe the Mueller potential.
 * If one is careful, future potentials could be designed in the same way so that they can be seamlessly integrated
 * into the main swarm files.  
 *
 * ***************************/

#include <vector>
#include <iostream>
#include <cmath>
#include "structs.h"
#include "mueller.h"

using namespace std;

mueller::mueller()
{//Mueller potential constructor
   gradient.resize(2);
   gradient[0] = 0.0;
   gradient[1] = 0.0;
   
   int i; //Loop indexing
   A.resize(4);
   b.resize(4);
   x0.resize(4);
   a.resize(4);
   c.resize(4);
   y0.resize(4);

   for(i = 0; i <= 3; i++)
   {//Initialize the constants
       if(i == 0)
       {
           A[i] = -200;
           b[i] = 0;
           x0[i] = 1;
           a[i] = -1;
           c[i] = -10;
           y0[i] = 0;
       }
       if(i == 1)
       {
           A[i] = -100;
           b[i] = 0;
           x0[i] = 0;
           a[i] = -1;
           c[i] = -10;
           y0[i] = 0.5;
       }
       if(i == 2)
       {
           A[i] = -170;
           b[i] = 11;
           x0[i] = -0.5;
           a[i] = -6.5;
           c[i] = -6.5;
           y0[i] = 1.5;
       }
       if(i == 3)
       {
           A[i] = 15;
           b[i] = 0.6;
           x0[i] = -1;
           a[i] = 0.7;
           c[i] = 0.7;
           y0[i] = 1;
       }
   }
}

double mueller::get_potential(double &x_point, double &y_point)
{
    double x = x_point;
    double y = y_point;
    double V_xy = 0.0; //Potential
    int i; //Loop indexing
    double V_temp;

    for(i = 0; i <= 3; i++)
    {
        //Potential is a sum from k=1 to k=4
        V_temp = A[i]*exp(a[i]*pow((x-x0[i]),2) + b[i]*(x-x0[i])*(y-y0[i]) + c[i]*pow((y-y0[i]),2));
        V_xy += V_temp;
    }

    return V_xy; 
}

void mueller::get_gradient(double &x_point, double &y_point)
{//Get the value of the gradient of the potential (the force)
    double x = x_point;
    double y = y_point;
    int i; //Loop indexing
    double dVx = 0.0, dVy = 0.0, dVx_temp = 0.0, dVy_temp = 0.0; //Gradient

    for(i = 0; i <= 3; i++)
    {
        //Gradient is a sum from k=1 to k=4
        dVx_temp = (A[i]*exp(a[i]*pow((x-x0[i]),2) + b[i]*(x-x0[i])*(y-y0[i]) + c[i]*pow((y-y0[i]),2)))*(2*a[i]*(x-x0[i]) + b[i]*(y-y0[i]));
        dVy_temp = (A[i]*exp(a[i]*pow((x-x0[i]),2) + b[i]*(x-x0[i])*(y-y0[i]) + c[i]*pow((y-y0[i]),2)))*(b[i]*(x-x0[i]) + 2*c[i]*(y-y0[i]));
        dVx += dVx_temp;
        dVy += dVy_temp;
    }
    gradient[0] = dVx;
    gradient[1] = dVy;
}


