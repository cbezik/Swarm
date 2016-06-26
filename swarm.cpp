#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include "swarm.h"
#include "structs.h"

using namespace std;

Swarm_method::Swarm_method()
{//Swarm method constructor
    sim_time = time(NULL); 
    start = clock();
    read_input(); //Scans in input and sets variables
    initialize(); //Initialize vectors and images
}

void Swarm_method::read_input()
{//Scans in input and sets variables
    string trash; //Holds unnecessary input
    ifstream input ("Input.swarm");
    if(input.is_open())
    {
        //Read input
        input >> sim.images; getline(input, trash); //Number of images
        input >> sim.load_string_flag; getline(input, trash); //Whether to load a string or not
    }
    else
    {
        cout << "No input file! Exiting..." << endl;
        exit(EXIT_FAILURE);
    }
}

void Swarm_method::initialize()
{//Initialize the vectors
    int i; //Loop indexing
    double x_increment, y_increment; //For interpolation

    swarm_string.resize(sim.images); //Resize the vector so that it contains the right number of images
    //Initialize the collective variables
    if(sim.load_string_flag == 1)
    {
        load_string(); //Load a full initial string
    }
    else
    {
        //Linearly interpolate string from initial and final states
        ifstream initial ("initial.state");
        ifstream end ("final.state");

        //Read initial state in
        initial >> swarm_string[0].x >> swarm_string[0].y;
        end >> swarm_string[sim.images-1].x >> swarm_string[sim.images-1].y;
        x_increment = (swarm_string[sim.images-1].x - swarm_string[0].x) / (sim.images-1); 
        y_increment = (swarm_string[sim.images-1].y - swarm_string[0].y) / (sim.images-1);
                
        //Linearly interpolate the CVs
        cout << "0 " << swarm_string[0].x << " " << swarm_string[0].y << endl; //Debugging
        for(i = 1; i < sim.images - 1; i++)
        {
            swarm_string[i].x = swarm_string[i-1].x + x_increment;
            swarm_string[i].y = swarm_string[i-1].y + y_increment;
            cout << i << " " << swarm_string[i].x << " " << swarm_string[i].y << endl; //Debugging
        }
        cout << sim.images-1 << " " << swarm_string[sim.images-1].x << " " << swarm_string[sim.images-1].y << endl; //Debugging
    }    
}

void Swarm_method::load_string()
{
    //Load a full string
}

void Swarm_method::run_swarm()
{//Run the simulation

    cout << "Running swarm..." << endl;
    return;
}


