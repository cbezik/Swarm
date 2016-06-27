#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
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
        getline(input, trash); getline(input, trash);
        
        input >> sim.images; getline(input, trash);
        input >> sim.load_string_flag; getline(input, trash); 
        input >> sim.iterations; getline(input, trash); 
        input >> sim.n_trajectories; getline(input, trash);
        input >> sim.l_trajectories; getline(input, trash);

        getline(input, trash); getline(input, trash);

        input >> sim.friction; getline(input, trash);
        input >> sim.time_step; getline(input, trash);
        input >> sim.mass; getline(input, trash);
        input >> sim.kT; getline(input, trash);
        
        getline(input, trash); getline(input, trash);

        input >> sim.potential_flag; getline(input, trash);
        //cout << sim.potential_flag << endl; //Debugging
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
    ifstream string_in ("string.out"); 

    swarm_string.resize(sim.images); //Resize the vector so that it contains the right number of images
    //Initialize the collective variables
    if(sim.load_string_flag == 1 && string_in.is_open())
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
        //cout << "0 " << swarm_string[0].x << " " << swarm_string[0].y << endl; //Debugging
        for(i = 1; i < sim.images - 1; i++)
        {
            swarm_string[i].x = swarm_string[i-1].x + x_increment;
            swarm_string[i].y = swarm_string[i-1].y + y_increment;
            //cout << i << " " << swarm_string[i].x << " " << swarm_string[i].y << endl; //Debugging
        }
        //cout << sim.images-1 << " " << swarm_string[sim.images-1].x << " " << swarm_string[sim.images-1].y << endl; //Debugging
    }   

    string_in.close();
}

void Swarm_method::load_string()
{
    //Load a full string
    ifstream string_in ("string.out"); 
    string trash; //Holds unnecessary input
    int i; //Loop indexing

    getline(string_in, trash); getline(string_in, trash);

    for(i = 0; i < sim.images; i++)
    {
        string_in >> swarm_string[i].x >> swarm_string[i].y; getline(string_in, trash);
        //cout << swarm_string[i].x << " " << swarm_string[i].y << endl; //Debugging
    }
    string_in.close(); 
}

void Swarm_method::write_string()
{
    ofstream string_out ("string.out");
    int i; //Loop indexing

    string_out << "X\t\t\t\tY" << endl; //Write header
    string_out << setprecision(6) << fixed << endl; //Set precision
    for(i = 0; i < sim.images; i++)
    {
        string_out << swarm_string[i].x << "\t\t" << swarm_string[i].y << endl;
    }

    string_out.close();
}

void Swarm_method::run_swarm()
{//Run the simulation
    int i; //Loop indexing
    int iter_counter; //Counting value

    iter_counter = 1; 

    cout << "Running swarm..." << endl;
    write_string(); //Writing the current string

    for(i = 1; i <= sim.iterations; i++)
    {
        //Run simulation
        restrained_sampling(); //Do restrained sampling for equilibrium states 
        generate_swarms(); //Generate swarms of unrestrained trajectories
        evolve_cv(); //Estimate drift and evolve CVs
        reparametrize(); //Reparametrize string

        if(iter_counter % 10 == 0)
        {
            cout << "Iteration: " << iter_counter << endl; 
        }
        
        iter_counter++; 
    }

    write_string(); 
}

void Swarm_method::restrained_sampling()
{
    int i;
    if(sim.potential_flag == 1)
    {
        //Do nothing - for Mueller potential the dynamics will act directly on the Cartesian coordinates
        //Simulate initialization error
        for(i = 0; i < sim.images; i++)
        {
            default_random_engine generator;
            normal_distribution<double> distributionx(swarm_string[i].x, 0.005);
            swarm_string[i].x = distributionx(generator); //Randomize initial coordinates, mean z_i, var 0.005
            normal_distribution<double> distributiony(swarm_string[i].y, 0.005); 
            swarm_string[i].y = distributiony(generator); //Randomize y coordinate
        }
    }
}

void Swarm_method::generate_swarms()
{
    int i, j, k; //Loop indexing
    double x_n, y_n, vx_n, vy_n; //Trajectory at time n
    double x_n1, y_n1, vx_n1, vy_n1; //Trajectory at time n+1
    double A_xn, A_yn; //Random term for Langevin dynamics
    for(i = 0; i < sim.images; i++)
    {
        for(j = 1; j <= sim.n_trajectories; j++)
        {
            //Run a trajectory, calculate the drift of the cv, and add it to the running total
            for(k = 1; k <= sim.l_trajectories; k++)
            {
                //cout << i << " " << j << " " << k << endl; //Debugging
                if(k==1)
                {
                    //Initialize the trajectory
                    x_n = swarm_string[i].x; //Start at the current string position
                    y_n = swarm_string[i].y;
                    vx_n = 0; //Start velocities at zero - random motion will initiate them
                    vy_n = 0; 
                }

                A_xn = 0.5*sim.time_step*sim.time_step*(1/sim.mass)*gradx_n
                x_n1 = x_n + sim.time_step*v_n + A_xn; 
            }

        }
    }
}

void Swarm_method::evolve_cv()
{
    return;
}

void Swarm_method::reparametrize()
{
    return;
}
