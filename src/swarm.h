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
#include <random>
#include "structs.h"
#include "mueller.h"
#include <mpi.h>

using namespace std;

class Swarm_method
{
    public:
        //Vectors and important objects
        clock_t start, end;
        time_t sim_time;
        parameters sim; 
        vector<image> swarm_string; //Vector containing all the images on the string
        mueller potential; 
        int iter_counter; //Check number of iterationsi
        int image_number; //Number of the image for parallel code
        
        //Function declarations

        Swarm_method(); //Constructor
        void run_swarm(); //Run the simulation
        void read_input(); //Read input
        void initialize(); //Initialize the system
        void load_string(); //Load a full initial string
        void write_string(); //Write the full string
        void restrained_sampling(); //Do restrained sampling for equilibrium states
        void generate_swarms(); //Generate swarms of unbiased trajectories
        void evolve_cv(); //Evolves the CVs using an estimate of the drift
        void reparametrize(); //Reparametrize the string
        void write_log(); //Write a simulation log
        void clear_output(); //Clear the output files
};

#endif
