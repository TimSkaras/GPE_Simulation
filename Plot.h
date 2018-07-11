#include <stdlib.h>
#include <stdio.h>

/* 

DESCRIPTION:
This file contains convenient functions for plotting data using gnuplot.

In order for these  functions to work, you must have install gnuplot on your computer and 
modify the enviroment variables of your computer so that the default search path variable
on your system includes the bin directory of gnuplot.

*/ 

struct plot_information
{
    int num_commands;
    char * output_file;
    double x_length;
    double y_length;
    double length;
    double T;
};

struct plot_information3D
{
    int num_commands;
    char * output_file;
    double x_length;
    double y_length;
    double z_length;
    double T;
};

void plotFunction(char * commandsForGnuplot[], int num_commands, double xvals[], double yvals[], int num_points, char * output_file ){   
    /*
    
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

void plotSurface2DFull(char * commandsForGnuplot[], struct plot_information info, int spx, int spy, int arg_time_points, double solution[spx][spx][arg_time_points], int time_index ){
    /*
    This function gives a 3-dimensional plot of data to 2+1D problem. Plots position versus time. Solution is specified using a 2-dimensional array where each row is a position 
    in space and each column is a point in time. Converts the solution from a matrix specifying the z-value at each point and converts it into a file listing each point by 
    its x,y,z coordinate

    INPUT: 
    commandsForGnuplot -- array of strings, element of which is a valid gnuplot command
    info -- struct containing relevant information to the plotter
    space_points -- integer describing the number of points in space, which should equal the number of rows in solution
    arg_time_points -- integer describing the number of points in time, which should equal the number of columns in solution
    solution -- 2-dimensional array matching the description provided above

    OUTPUT:
    Returns nothing but generates the desired plot in a new window

    EXAMPLE INPUTS:
    char * commandsForGnuplot[] = {"set title \"TITLE\"", "splot 'data.temp' with lines"};
    double solution[][2] = {{1,2,3}, {2,4,6}};
    */

    int index_y = (int) ceil(spy/50.);
    int new_y = (int) floor(spy /  ceil(spy/50.));
    int index_x = (int) ceil(spx/500.);
    int new_x = (int) floor(spx /  ceil(spx/500.));
    double reduced_solution[new_x][new_y];

    double HX = info.x_length/spx;
    double HY = info.y_length/spy;
    

    // If solution has too many points, we need to take some out so that the plot has perceptible contours
    // There should always be more time points than space points

    // printf("index_y: %d \n T_new: %d \n", index_y, new_y);
    if (spx * spy > 1500)
    {   

        for (int i = 0; i < new_x; ++i){

            for (int p = 0; p < new_y; ++p)
                reduced_solution[i][p] = solution[i * index_x][p * index_y][time_index];
        }
    }

    FILE * temp = fopen(info.output_file, "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
    

    if (spx * spy <= 1500)
    {
        for (int i=0; i < spx; i++)
        {
            for (int j = 0; j < spy; ++j)
                    fprintf(temp, "%e %e %e \n", HX * i, HY * j , solution[i][j][time_index]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }
    else{

        for (int i=0; i < new_x; i++){

            for (int j = 0; j < new_y; ++j)
                fprintf(temp, "%e %e %e \n", index_x * HX * i, index_y * HY * j , reduced_solution[i][j]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }

    fclose(temp);

    for (int i=0; i < info.num_commands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
}

void plotSurface2DReduced(char * commandsForGnuplot[], struct plot_information info, int space_points, int arg_time_points, double solution[][arg_time_points] ){
    /*
    This function gives a 3-dimensional plot of data to 2+1D problem. Plots position versus time. Solution is specified using a 2-dimensional array where each row is a position 
    in space and each column is a point in time. Converts the solution from a matrix specifying the z-value at each point and converts it into a file listing each point by 
    its x,y,z coordinate

    INPUT: 
    commandsForGnuplot -- array of strings, element of which is a valid gnuplot command
    info -- struct containing relevant information to the plotter
    space_points -- integer describing the number of points in space, which should equal the number of rows in solution
    arg_time_points -- integer describing the number of points in time, which should equal the number of columns in solution
    solution -- 2-dimensional array matching the description provided above

    OUTPUT:
    Returns nothing but generates the desired plot in a new window

    EXAMPLE INPUTS:
    char * commandsForGnuplot[] = {"set title \"TITLE\"", "splot 'data.temp' with lines"};
    double solution[][2] = {{1,2,3}, {2,4,6}};
    */

    double Dt = info.T / arg_time_points;
    double H = info.length/space_points;

    double time_samples = 100.;
    double space_samples = 200.;

    int index_t = (int) ceil(arg_time_points/time_samples);
    int new_time_points = (int) floor(arg_time_points /  ceil(arg_time_points/time_samples));
    int index_s = (int) ceil(space_points/space_samples);
    int new_spatial_points = (int) floor(space_points /  ceil(space_points/space_samples));
    double reduced_solution[new_spatial_points + 1][new_time_points];
    

    // If solution has too many points, we need to take some out so that the plot has perceptible contours
    // There should always be more time points than space points

    // printf("index_t: %d \n T_new: %d \n", index_t, new_time_points);
    if (space_points * arg_time_points > 1500)
    {   
        for (int i = 0; i < new_spatial_points; ++i){

            for (int p = 0; p < new_time_points; ++p)
                reduced_solution[i][p] = solution[i * index_s][p * index_t];
        }
    }

    // Add last row for cosmetic purposes
    for (int j = 0; j < new_time_points; ++j)
    {
        reduced_solution[new_spatial_points][j] = reduced_solution[0][j];
    }

    FILE * temp = fopen(info.output_file, "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
    

    if (space_points * arg_time_points <= 1500)
    {
        for (int i=0; i < space_points; i++)
        {
            for (int j = 0; j < arg_time_points; ++j)
                    fprintf(temp, "%e %e %.8e \n", H * i, Dt * j , solution[i][j]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }
    else{

        for (int i=0; i < new_spatial_points + 1; i++){

            for (int j = 0; j < new_time_points; ++j)
                fprintf(temp, "%e %e %.8e \n", index_s * H * i, index_t * Dt * j , reduced_solution[i][j]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }

    fclose(temp);

    for (int i=0; i < info.num_commands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
}

void plotSurface3DReduced(char * commandsForGnuplot[], struct plot_information3D info, int space_points, int arg_time_points, double solution[][arg_time_points] ){
    /*
    This function gives a 3-dimensional plot of data to 3+1D problem. Plots position versus time. Solution is specified using a 2-dimensional array where each row is a position 
    in space and each column is a point in time. Converts the solution from a matrix specifying the z-value at each point and converts it into a file listing each point by 
    its x,y,z coordinate

    INPUT: 
    commandsForGnuplot -- array of strings, element of which is a valid gnuplot command
    num_commands -- an integer that gives the number of strings in the commandsForGnuplot array
    solution -- 2-dimensional array matching the description provided above
    space_points -- integer describing the number of points in space, which should equal the number of rows in solution
    arg_time_points -- integer describing the number of points in time, which should equal the number of columns in solution
    T -- double representing length of time simulated    
    output_file -- filename with extension included that gives the name of the desired file 

    OUTPUT:
    Returns nothing but generates the desired plot in a new window

    EXAMPLE INPUTS:
    char * commandsForGnuplot[] = {"set title \"TITLE\"", "splot 'data.temp' with lines"};
    double solution[][2] = {{1,2,3}, {2,4,6}};
    */

    double Dt = info.T / arg_time_points;
    double HZ = info.z_length/space_points;

    double time_samples = 70.;
    double space_samples = 400.;

    int index_t = (int) ceil(arg_time_points/time_samples);
    int new_time_points = (int) floor(arg_time_points /  ceil(arg_time_points/time_samples));
    int index_s = (int) ceil(space_points/space_samples);
    int new_spatial_points = (int) floor(space_points /  ceil(space_points/space_samples));
    double reduced_solution[new_spatial_points + 1][new_time_points];
    

    // If solution has too many points, we need to take some out so that the plot has perceptible contours
    // There should always be more time points than space points

    // printf("index_t: %d \n T_new: %d \n", index_t, new_time_points);
    if (space_points * arg_time_points > 1500)
    {   
        for (int i = 0; i < new_spatial_points; ++i){

            for (int p = 0; p < new_time_points; ++p)
                reduced_solution[i][p] = solution[i * index_s][p * index_t];
        }
    }

    // Add last row for cosmetic purposes
    for (int j = 0; j < new_time_points; ++j)
    {
        reduced_solution[new_spatial_points][j] = reduced_solution[0][j];
    }

    FILE * temp = fopen(info.output_file, "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen("gnuplot -persistent", "w");
    

    if (space_points * arg_time_points <= 1500)
    {
        for (int i=0; i < space_points; i++)
        {
            for (int j = 0; j < arg_time_points; ++j)
                    fprintf(temp, "%e %e %.8e \n", HZ * i, Dt * j , solution[i][j]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }
    else{

        for (int i=0; i < new_spatial_points + 1; i++){

            for (int j = 0; j < new_time_points; ++j)
                fprintf(temp, "%e %e %.8e \n", index_s * HZ * i, index_t * Dt * j , reduced_solution[i][j]); //Write the data to a temporary file

            fprintf(temp, "\n");
        }
    }

    fclose(temp);

    for (int i=0; i < info.num_commands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
}