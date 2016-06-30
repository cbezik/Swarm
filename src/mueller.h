#ifndef __MUEL_H__
#define __MUEL_H__

/********* mueller.h ***********
 * Author: Cody Bezik
 * Last modified: 6-27-2016
 * Project: swarm.cpp
 * 
 * Description:
 * This file contains the declaration of functions and parameters necessary to describe the Mueller potential.
 * If one is careful, future potentials could be designed in the same way so that they can be seamlessly integrated
 * into the main swarm files.  
 *
 * ***************************/

#include <vector>
#include "structs.h"

using namespace std;

class mueller
{//Describes the Mueller potential

    public:

        //Objects
        vector<double> gradient; //Contains x and y component of the gradient

        //Constants for the potential
        vector<double> A;// = vector<double>(4); 
        vector<double> b;// = vector<double>(4);
        vector<double> x0;// = vector<double>(4);
        vector<double> a;// = vector<double>(4);
        vector<double> c;// = vector<double>(4);
        vector<double> y0;// = vector<double>(4);


        //Functions
        mueller(); //Constructor
        double get_potential(double &x_point, double &y_point); //Returns the value of the potential
        void get_gradient(double &x_point, double &y_point); //Returns the x and y component of the gradient
};

#endif
