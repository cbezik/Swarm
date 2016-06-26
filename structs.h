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
};

struct image
{//The images along the string
    double x, y; //For Mueller potential the collective variables are just the x and y coordinates
};


#endif