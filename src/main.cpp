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


using namespace std;

//Runs the swarm of trajectories method

int main(int argc, char **argv)
{
    Swarm_method *myswarm;
    myswarm = new Swarm_method();
    myswarm->run_swarm();

    delete myswarm;
    return 0;
}
