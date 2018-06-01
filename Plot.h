#include <stdlib.h>
#include <stdio.h>

/* 

TODO:
3) Create function for phase space maybe

DESCRIPTION:
This file contains convenient functions for plotting data using gnuplot.

In order for these  functions to work, you must have install gnuplot on your computer and 
modify the enviroment variables of your computer so that the default search path variable
on your system includes the bin directory of gnuplot.

*/ 

void plotFunction(char * commandsForGnuplot[], int num_commands, double xvals[], double yvals[], int num_points, char * output_file )
{   /*
    
    This function takes a set of commands, data, and a filename and produces a 2D plot using gnuplot
    based on these data

    INPUT: 
    commandsForGnuplot -- array of strings, element of which is a valid gnuplot command
    num_commands -- an integer that gives the number of strings in the commandsForGnuplot array
    xvals -- array of doubles that gives the x coordinate of each data point
    yvals -- array of doubles that gives the y coordinate of each data point
    num_points -- an integer which gives the number of elements in xvals and yvals (should be the same)
    output_file -- filename with extension included that gives the name of the desired file 

    OUTPUT:
    Returns nothing but generates the desired plot in a new window

    EXAMPLE INPUTS:
    char * commandsForGnuplot[] = {"set title \"TITLE\"", "plot [0:6] [0:7] 'data.temp'"};
    double xvals[NUM_POINTS] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double yvals[NUM_POINTS] = {5.0 ,3.0, 1.0, 3.0, 5.0};
    
    */

    FILE * temp = fopen(output_file, "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
    
    int i;
    for (i=0; i < num_points; i++)
    {
    fprintf(temp, "%e %e\n", xvals[i], yvals[i]); //Write the data to a temporary file
    }

    for (i=0; i < num_commands; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }

}

void plotHarmonicEnergy(char * method, double xvals[], double position[], double velocity[], int num_points){
    /*
    Plots the total energy at each point in time for a harmonic oscillator where mass and spring constant are unity

    method -- string indicating the name of the method used so a title can be given to the plot
    xvals --  array of doubles that gives the x coordinate of each point in time
    position -- array the gives the position of the oscillator at each point of discretised time
    velocity -- array that gives the velocity of the oscillator at teach point of discretised time
    num_points -- integer that gives the numeber points in each array (should be the same number of elements in each array)

    */

    char energyCommand[80];
    char title[80];
    strcpy(energyCommand, "plot 'energy.temp' with lines");
    sprintf(title, "set title \"Undamped Harmonic Oscillator (%s) - Energy Loss\"", method);
    char * energyCommands[] = {title, "set xlabel \"Time\"", "show xlabel", "set ylabel \"Energy\"", "show ylabel" , energyCommand};
    int num_commands = 6;

    double evals[num_points];

    for (int i = 0; i < num_points; ++i)
        evals[i] = pow(position[i], 2) + pow(velocity[i], 2); 

    plotFunction(energyCommands, num_commands, xvals, evals, num_points, "energy.temp");

}



// void plotSurface(char * commandsForGnuplot[], int num_commands, double solution[][], int spatial_points, int time_points, char * output_file, double delta_x, double delta_t ){


//     FILE * temp = fopen(output_file, "w");
//     /*Opens an interface that one can use to send commands as if they were typing into the
//      *     gnuplot command line.  "The -persistent" keeps the plot open even after your
//      *     C program terminates.
//      */
//     FILE * gnuplotPipe = _popen("gnuplot -persistent", "w");
    

//     for (int i=0; i <= spatial_points; i++)
//     {
//         for (int j = 0; j < time_points; ++j){
//             fprintf(temp, "%lf %lf %lf \n", delta_x * i, delta_t * j , u_solution[i][j]); //Write the data to a temporary file
//         }
        
//         fprintf(temp, "\n");
//     }

//     fclose(temp);


//     for (int i=0; i < num_commands; i++)
//     {
//         fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
//     }
// }