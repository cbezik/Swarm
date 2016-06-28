/********* swarm.cpp ***********
 * Author: Cody Bezik
 * Last modified: 6-27-2016
 * Project: swarm.cpp
 * 
 * Description:
 * This file contains the bulk of the function definitions required to run the swarm of trajectories method
 * on a 2D system evolving on the Mueller potential according to Langevin dynamics.
 *
 * ***************************/


#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include "swarm.h"
#include "structs.h"
#include "mueller.h"
#include "spline.h" i


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
        input >> sim.tolerance; getline(input, trash); 

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
    double dx, dy, dr2;

    string_out << "X\t\t\t\tY\t\t\t\tLinear Distance" << endl; //Write header
    string_out << setprecision(6) << fixed << endl; //Set precision
    for(i = 0; i < sim.images; i++)
    {
        if(i == 0)
        {
            string_out << swarm_string[i].x << "\t\t" << swarm_string[i].y << "\t\t" << "0" << endl;
        }
        else
        {
            dx = swarm_string[i].x - swarm_string[i-1].x;
            dy = swarm_string[i].y - swarm_string[i-1].y;
            dr2 = dx*dx + dy*dy;
            string_out << swarm_string[i].x << "\t\t" << swarm_string[i].y << "\t\t" << sqrtl(dr2) << endl;
        }
    }

    string_out.close();
}

void Swarm_method::run_swarm()
{//Run the simulation
    int i, j; //Loop indexing
    vector<double> displacements_x = vector<double>(sim.images); //For checking displacement over the course of an interation
    vector<double> displacements_y = vector<double>(sim.images);

    iter_counter = 1; 

    cout << "Starting swarm..." << endl;
    write_string(); //Writing the current string

    for(i = 1; i <= sim.iterations; i++)
    {
        //Run simulation
        for(j = 0; j <= sim.images; j++)
        {
            displacements_x[j] = swarm_string[j].x;
            displacements_y[j] = swarm_string[j].y; 
        }
        restrained_sampling(); //Do restrained sampling for equilibrium states 
        generate_swarms(); //Generate swarms of unrestrained trajectories
        evolve_cv(); //Estimate drift and evolve CVs
        reparametrize(); //Reparametrize string

        if(iter_counter % 10 == 0)
        {
            cout << "Iteration: " << iter_counter << endl; 
            write_string();
            end = clock();
            cout << "Time elapsed: " << (double)(end - start) / CLOCKS_PER_SEC << " Seconds" << endl;
        }

        //Check convergence
        for(j = 0; j <= sim.images; j++)
        {
            displacements_x[j] *= -1;
            displacements_y[j] *= -1;
            displacements_x[j] += swarm_string[j].x; //Displacement equal x(n+1) - x(n), etc
            displacements_y[j] += swarm_string[j].y; 
        }
        if(abs(*max_element(displacements_x.begin(), displacements_x.end())) < sim.tolerance && abs(*max_element(displacements_y.begin(), displacements_y.end())) < sim.tolerance)
        {
            cout << "Converged at iteration " << iter_counter << endl;
            break; //Exit loop if converged
        }
        
        iter_counter++;
    }

    write_string();

    end = clock(); 
    write_log(); 
}

void Swarm_method::restrained_sampling()
{//Do restrained sampling, first to find a system with atomic coordinates near the cvs, then to prepare initial configurations for swarms
    int i;
    default_random_engine generator;

    if(sim.potential_flag == 1)
    {
        //Do nothing - for Mueller potential the dynamics will act directly on the Cartesian coordinates
        //Simulate initialization error
        /*for(i = 0; i < sim.images; i++)
        {
            normal_distribution<double> distributionx(swarm_string[i].x, 0.005);
            swarm_string[i].x = distributionx(generator); //Randomize initial coordinates, mean z_i, var 0.005
            normal_distribution<double> distributiony(swarm_string[i].y, 0.005); 
            swarm_string[i].y = distributiony(generator); //Randomize y coordinate
        }*/
    }
}

void Swarm_method::generate_swarms()
{//Generate swarms of trajectories and calculate the drift of the CVs 
    int i, j, k; //Loop indexing
    double x_n, y_n, vx_n, vy_n; //Trajectory at time n
    double x_n1, y_n1, vx_n1, vy_n1; //Trajectory at time n+1
    double grad_xn, grad_yn; //X and Y component of the gradient
    double grad_xn1, grad_yn1; //Components of the gradient at next time step
    double A_xn, A_yn; //Random term for Langevin dynamics
    double sigma = sqrtl(2*sim.kT*sim.friction*(1/sim.mass)); //Parameter for langevin dynamics
    default_random_engine generator;
    normal_distribution<double> normal_dist(0, 1); //Random number generator, mean 0, variance 1
    double zeta_xn, eta_xn, zeta_yn, eta_yn; //Random numbers
    string debug; //Debugging

    //Zero drift before starting the calculation
    for(i = 0; i < sim.images; i++)
    {
        swarm_string[i].driftx = 0.0;
        swarm_string[i].drifty = 0.0; 
    }

    for(i = 0; i < sim.images; i++)
    {
        //For each image...
        for(j = 1; j <= sim.n_trajectories; j++)
        {
            //...run n_trajectories number of trajectories...
            for(k = 1; k <= sim.l_trajectories; k++)
            {
                //...each of length l_trajectories
                
                //cout << i << " " << j << " " << k << endl; //Debugging
                if(k==1)
                {
                    normal_distribution<double> distributionx(swarm_string[i].x, 0.005);
                    normal_distribution<double> distributiony(swarm_string[i].y, 0.005);
                    //Initialize the trajectory
                    //Normally would pluck from trajectories generated during restrained dynamics, but for model potential, instead we simulate the initialization error
                    x_n = distributionx(generator); //Start at the current string position
                    y_n = distributiony(generator);
                    vx_n = 0; //Start velocities at zero - random motion will initiate them
                    vy_n = 0; 
                }

                //Generate random numbers
                zeta_xn = normal_dist(generator);
                zeta_yn = normal_dist(generator);
                eta_xn = normal_dist(generator);
                eta_yn = normal_dist(generator);

                //cout << zeta_xn << " " << zeta_yn << " " << endl; //Debugging
               
                /*if(i == 23)
                {
                    cout << x_n << " " << y_n << endl;
                }*/
                potential.get_gradient(x_n, y_n); //Calculate and set gradient
                grad_xn = -potential.gradient[0];
                grad_yn = -potential.gradient[1];

                /*if( i == 23)
                {
                    cout << "Image " << i << " Trajectory  " << j << " Time step " << k << endl;
                    cout << "Gradient x " << grad_xn << " Gradient y " << grad_yn << endl; 
                }*/

                //cout << grad_xn << " " << grad_yn << " " << endl; //Debugging
                //cin >> debug; //Debugging

                //Calculate A constant
                A_xn = 0.5*sim.time_step*sim.time_step*((1/sim.mass)*grad_xn-sim.friction*vx_n)+sigma*pow(sim.time_step, 1.5)*(0.5*zeta_xn+(1/(2*sqrtl(3)))*eta_xn);
                A_yn = 0.5*sim.time_step*sim.time_step*((1/sim.mass)*grad_yn-sim.friction*vy_n)+sigma*pow(sim.time_step, 1.5)*(0.5*zeta_yn+(1/(2*sqrtl(3)))*eta_yn);

                //Progress forward the positions
                x_n1 = x_n + sim.time_step*vx_n + A_xn; 
                y_n1 = y_n + sim.time_step*vy_n + A_yn; 

                //Calculate the potential at the new positions
                potential.get_gradient(x_n1, y_n1);
                grad_xn1 = -potential.gradient[0];
                grad_yn1 = -potential.gradient[1];

                //Calculate the new velocity
                vx_n1 = vx_n + 0.5*sim.time_step*((grad_xn1+grad_xn)/sim.mass)-sim.time_step*sim.friction*vx_n+sqrt(sim.time_step)*sigma*zeta_xn-sim.friction*A_xn;
                vy_n1 = vy_n + 0.5*sim.time_step*((grad_yn1+grad_yn)/sim.mass)-sim.time_step*sim.friction*vy_n+sqrt(sim.time_step)*sigma*zeta_yn-sim.friction*A_yn;

                //Prepare for next iteration - move time step forward
                x_n = x_n1;
                y_n = y_n1;
                vx_n = vx_n1;
                vy_n = vy_n1;

                if(k == sim.l_trajectories)
                {//After the last trajectory, record the amount of drift in the CV
                    swarm_string[i].driftx += (x_n-swarm_string[i].x);
                    swarm_string[i].drifty += (y_n-swarm_string[i].y); 
                    //cout << swarm_string[i].driftx << " " << swarm_string[i].drifty << endl; //Debugging 
                }
            }
        }
    }
}

void Swarm_method::evolve_cv()
{//Average the drift and update the CVs
    int i; //Loop indexing
    for(i = 0; i < sim.images; i++)
    {
        swarm_string[i].driftx /= sim.n_trajectories; //Average the drift over the number of trajectories
        swarm_string[i].drifty /= sim.n_trajectories;
        swarm_string[i].x += swarm_string[i].driftx; //Update the CV's
        swarm_string[i].y += swarm_string[i].drifty;
    }
}

void Swarm_method::reparametrize()
{//Enforce a cubic spline interpolation onto a regular mesh

    vector<double> arc_lengths = vector<double>(sim.images); //Holds the arc length between images
    vector<double> alpha_star = vector<double>(sim.images); //Irregular mesh, normalized arc lengths
    vector<double> alpha = vector<double>(sim.images); //Regular mesh for interpolation

    int i; //Loop indexing
    double dx, dy, distance;

    for(i = 0; i < sim.images; i++)
    {
        if(i == 0)
        {
            arc_lengths[i] = 0.0;
        }
        else
        {
            dx = swarm_string[i].x - swarm_string[i-1].x;
            dy = swarm_string[i].y - swarm_string[i-1].y;
            distance = sqrtl(dx*dx+dy*dy); //Euclidean norm of z[i] - z[i-1]
            arc_lengths[i] = arc_lengths[i-1] + distance;
            //cout << arc_lengths[i] << endl;
        }
    }

    for(i = 0; i < sim.images; i++)
    {
        alpha_star[i] = arc_lengths[i] / arc_lengths[sim.images - 1]; //Normalize arc lengths to get alpha star
    }

    tk::spline sx; //Cubic spline in x
    tk::spline sy; //Cubic spline in y

    vector<double> image_x = vector<double>(sim.images); //Can only interpolate in one dimension at a time
    vector<double> image_y = vector<double>(sim.images);

    for(i = 0; i < sim.images; i++)
    {
        image_x[i] = swarm_string[i].x; 
        image_y[i] = swarm_string[i].y;
        //cout << image_x[i] << endl;
    }
    //exit(EXIT_FAILURE); 

    sx.set_points(alpha_star, image_x); //Give s alphas to generate images
    sy.set_points(alpha_star, image_y);

    for(i = 0; i < sim.images; i++)
    {//Distribute points on even mesh using cubic spline interpolation
        double alpha_i = (double)i / (sim.images - 1); //Regular mesh points
        //cout << alpha_i << endl; 
        double x = sx(alpha_i);
        double y = sy(alpha_i);
        swarm_string[i].x = x;
        swarm_string[i].y = y;
        //cout << x << " " << y << endl;
    }
    //exit(EXIT_FAILURE); 

    //Enforce a piecewise linear interpolation such that all points are linearly equidistant (in each collective variable)

    /*int i; //Loop indexing
    int k; //For reparametrization algorithm
    double dx, dy, dr2; //X distance, y distance, and square distance

    //Reparameterize
    vector<double> length = vector<double>(sim.images); //The length function is the length, up to each image, of the piecewise linear curve interpolated through the images
    length[0] = 0; //By definition
    for(i = 1; i < sim.images; i++)
    {
        dx = swarm_string[i].x - swarm_string[i-1].x;
        dy = swarm_string[i].y - swarm_string[i-1].y; 
        dr2 = dx*dx+dy*dy;
        length[i] = length[i-1] + sqrtl(dr2); //Build up the lengthx vector
    }
    vector<double> sm = vector<double>(sim.images);//sm vector defined in the algorithm of eq 49 and 50 of Maragliano et al (2006)
    sm[0] = 0; //Don't use this value
    for(i = 1; i < sim.images; i++)
    {
        sm[i] = (i)*length[sim.images-1] / (sim.images - 1); //Build up the sm vector
    }
    //Reparametrize points, excluding the ends
    for(i = 1; i < sim.images - 1; i++)
    {
        //Select k such that L(k-1) < s(m) <= L(k)
        k = 1; 
        while(!(length[k-1] < sm[i] && sm[i] <= length[k]))
        {
            k++; //Find the correct value of k
        }
        //Reparametrize
        dx = swarm_string[k].x - swarm_string[k-1].x;
        dy = swarm_string[k].y - swarm_string[k-1].y;
        dr2 = dx*dx+dy*dy;
        swarm_string[i].x = swarm_string[k-1].x + (sm[i]-length[k-1])*((swarm_string[k].x - swarm_string[k-1].x)/(sqrtl(dr2))); //Linear reparametrization algorithm
        swarm_string[i].y = swarm_string[k-1].y + (sm[i]-length[k-1])*((swarm_string[k].y - swarm_string[k-1].y)/(sqrtl(dr2))); 
    }

    //Reparametrize y
    vector<double> lengthy = vector<double>(sim.images); //The length function is the length, up to each image, of the piecewise linear curve interpolated through the images
    lengthy[0] = 0; //By definition
    for(i = 1; i < sim.images; i++)
    {
        lengthy[i] = lengthy[i-1] + (swarm_string[i].y - swarm_string[i-1].y); //Build up the lengthx vector
    }
    vector<double> sm_y = vector<double>(sim.images);//sm vector defined in the algorithm of eq 49 and 50 of Maragliano et al (2006)
    sm_y[0] = 0; //Don't use this value
    for(i = 1; i < sim.images; i++)
    {
        sm_y[i] = (i)*lengthy[sim.images-1] / (sim.images - 1); //Build up the sm vector
    }
    //Reparametrize points, excluding the ends
    for(i = 1; i < sim.images - 1; i++)
    {
        //Select k such that L(k-1) < s(m) <= L(k)
        k = 1; 
        while(!(lengthy[k-1] < sm_y[i] && sm_y[i] <= lengthy[k]))
        {
            k++; //Find the correct value of k
        }
        //Reparametrize
        swarm_string[i].y = swarm_string[k-1].y + (sm_y[i]-lengthy[k-1])*((swarm_string[k].y - swarm_string[k-1].y)/(abs(swarm_string[k].y - swarm_string[k-1].y))); //Linear reparametrization algorithm
    }
    */
}

void Swarm_method::write_log()
{//Write a simulation log
    ofstream logout ("simul_log.out"); //Stream for writing
    int days, hours, minutes;
    double seconds, rem;
    rem = ((double)(end-start)) / CLOCKS_PER_SEC; //For tracking time
    days = (int)floor(rem/86400);
    rem = rem - days*86400;
    hours = (int)floor(rem/3600);
    rem = rem - 3600*hours; 
    minutes = (int)floor(rem/60);
    seconds = rem - 60*minutes;

    logout << sim.images << " Images" << endl;
    logout << iter_counter << " Iterations" << endl;
    logout << sim.n_trajectories << " Trajectories per swarm" << endl;
    logout << sim.l_trajectories << " Time steps per trajectory" << endl;

    logout << "Langevin dynamics, with: " << endl;
    logout << sim.friction << " Friction coefficient" << endl;
    logout << sim.time_step << " Time step" << endl;
    logout << sim.mass << " Particle mass" << endl;
    logout << sim.kT << " kT" << endl;

    if(sim.potential_flag == 1)
    {
        logout << "Mueller potential" << endl;
    }

    logout << "Elapsed time was: " << days << " Days " << hours << " Hours " << minutes << " Minutes " << seconds << " Seconds " << endl; 
}
