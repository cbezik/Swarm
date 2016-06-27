#ifndef __STRUCT_H__
#define __STRUCT_H__

/************* structs.h ************
 * Author: Cody Bezik
 * Last modified: 6-26-2016
 * Project: swarm.cpp
 *
 * Description:
 * This simple header file has a few structures necessary that are key to the swarm method.
 *
 *
 * ********************************/

using namespace std;

struct parameters
{//Holds simulation parameters
    int images; //Number of images for the string
    int load_string_flag; //Whether to load a full string or not
    int iterations; //Max number of iterations to simulate
    int n_trajectories; //Number and length of trajectories
    int l_trajectories;
    double friction; //Friction coefficient for langevin dynamics
    double time_step; //Time step for Langevin dynamics
    double mass; //Mass for LD
    double kT; //Temperature for LD
    int potential_flag; //Potential flag, 1 = Mueller, 
};

struct image
{//The images along the string
    double x, y; //For Mueller potential the collective variables are just the x and y coordinatesi
    double driftx, drifty; //Estimates of the drift of the collective variables
};


#endif
