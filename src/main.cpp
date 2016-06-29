/**************** main.cpp **************
 * Author: Cody Bezik
 * Last modified: 6-26-2016
 * Project: Swarm.cpp
 *
 * Description:
 *
 * Runs the swarm of trajectories method on the specified system.  Currently this testing version only runs
 * the Mueller potential.  Users can input via Input.swarm the simulation parameters and specify via the command
 * line the number of nodes and images (which must be equal) to run the string method on.
 *
 * *************************************/
#include "swarm.h"
#include "structs.h"
#include <mpi.h>
#include "mueller.h"
#include <iostream>

using namespace std;

//Runs the swarm of trajectories method

int main(int argc, char **argv)
{
    int image, num_tasks; //For MPI
    Swarm_method *myswarm; 
    myswarm = new Swarm_method(); 
    
    //Begin parallel code

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &image);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks); 

    if(num_tasks != myswarm->sim.images)
    {//Break if number of tasks isn't equal to the number of nodes
        cout << "Number of images not equal to number of processors" << endl;
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(EXIT_FAILURE);
    }

    myswarm.image_number = image; //Identify each image 
    myswarm->run_swarm();

    MPI_Finalize();

    //End parallel code

    cout << "Simulation done" << endl;
    delete myswarm;
    return 0;
}
