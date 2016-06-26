#ifndef __SWARM_H__
#define __SWARM_H__

/**************  swarm.h **********
 *
 * Author: Cody Bezik
 * Last modified: 06-26-2016
 * Project: swarm.cpp
 *
 * Description:
 * This header file contains declarations of the critical functions to run a swarm simulation.
 * The definitions of these functions is included in swarm.cpp.
 *
 * ********************************/

#include <ctime>
#include <vector>
#include "structs.h"

using namespace std;

class Swarm_method
{
    public:

    //Vectors and important objects
    clock_t start, end;
    time_t sim_time;
    parameters sim; 
    vector<image> swarm_string; //Vector containing all the images on the string

    //Function declarations

    Swarm_method(); //Constructor
    void run_swarm(); //Run the simulation
    void read_input(); //Read input
    void initialize(); //Initialize the system
    void load_string(); //Load a full initial string
};

#endif
