/*
This program numerically finds the thomas fermi ground state, solves the nonlinear schrodinger equation with harmonic potential in 2D using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant W is the harmonic potential constant omega

NOTES:


IMPORTANT RESULTS:
- When nonlinearity is turned off, the change in groundstate (30000 iterations) under real time propagation is less than 1 in 10^8
- When nonlinearity is small (G~.01), oscillation in g.s. is not bad (1 in e-4) and normalization err is not bad either (1 in e-10) (time step = e-5)
- When nonlinearity is small and (time step = 5e-6), norm error is order of magnitude better but oscillation in g.s. is no different
- Now try case where we run imag time prop algorithm for twice as many iterations --- same result
- Now try larger time step for imag time prop. Increased time step for ITP by factor of 4 -- result is unstable
- Smaller time step for ITP by factor of 2, amplitude of oscillation is about half -- result shows promise


TODO:
1) Change the ITP from RK4 to standard FD -- Done
2) Modify RTP to have 2 RHS functions instead of 4
3) Transfer RTP code into its own function
4) Create stepTwo function for 2D

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 0.00125
#define XLENGTH 6.0 //4.0
#define YLENGTH 0.15 //0.1
#define TIME_POINTS 250 	//number of time points
#define SPACE_POINTS 150	//number of spatial points
#define SPX 800 //600
#define SPY 20 //20

#define Dxx(array, x, y, pee) ( (-1. * array[mod(x + 2, SPX)][y][pee] + 16.* array[mod(x + 1, SPX)][y][pee] - \
		30. * array[mod(x , SPX)][y][pee] + 16. * array[mod(x - 1, SPX)][y][pee] +  -1 * array[mod(x - 2, SPX)][y][pee]) / (12. * pow(HX, 2)) )
;

#define Dyy(array, x, y, pee) ( (-1. * array[x][mod(y + 2, SPY)][pee] + 16.* array[x][mod(y + 1, SPY)][pee] - \
		30. * array[x][mod(y , SPY)][pee] + 16. * array[x][mod(y - 1, SPY)][pee] +  -1 * array[x][mod(y - 2, SPY)][pee]) / (12. * pow(HY, 2)) )
;

const double HX = XLENGTH / (SPX);
const double HY = YLENGTH / (SPY);
const double Dt = TIME / (TIME_POINTS);
const double WX = 7 * 2. * PI;
const double WY = 476 * 2. * PI;
const double G = 1000.;

// How many time steps do we want
const int time_points = TIME_POINTS;

// How many spatial points are there
const int spx = SPX;
const int spy = SPY;

// ------------------ Auxiliary Functions ----------------

int mod(int a, int b){
	// This is a standard modulus function -- always outputs a positive number
	// E.g., mod(-1, 3) = 2

    int r = a % b;
    return r < 0 ? r + b : r;
}

void plotSurface(char * commandsForGnuplot[], int num_commands, double solution[SPX][SPY][TIME_POINTS], int spx, int spy, int time_index, char * output_file ){
	/*
	This function gives a 3-dimensional plot of data to 1+1D problem. Plots position versus time. Solution is specified using a 2-dimensional array where each row is a position 
	in space and each column is a point in time. Converts the solution from a matrix specifying the z-value at each point and converts it into a file listing each point by 
	its x,y,z coordinate

    INPUT: 
    commandsForGnuplot -- array of strings, element of which is a valid gnuplot command
    num_commands -- an integer that gives the number of strings in the commandsForGnuplot array
   	solution -- 2-dimensional array matching the description provided above
	spx -- integer describing the number of points in x dimension
	spy -- integer describing the number of points in y dimension
	time_index -- integer describing the desired point in time to be plotted  
	output_file -- filename with extension included that gives the name of the desired file 

    OUTPUT:
    Returns nothing but generates the desired plot in a new window

    EXAMPLE INPUTS:
    char * commandsForGnuplot[] = {"set title \"TITLE\"", "splot 'data.temp' with lines"};
    double solution[][2] = {{1,2,3}, {2,4,6}};
	*/

    int index_y = (int) ceil(spy/50.);
	int new_y = (int) floor(spy /  ceil(spy/50.));
	int index_x = (int) ceil(spx/90.);
	int new_x = (int) floor(spx /  ceil(spx/90.));
    double reduced_solution[new_x][new_y];
    

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

    FILE * temp = fopen(output_file, "w");
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

    for (int i=0; i < num_commands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
}

void plotNormalization(double full_solution[SPX][SPY][TIME_POINTS]){
	/*
	* This function takes the squared solution matrix, finds the normalization at each point in time and plots the error (normalization should equal unity)
	*/

	double norm[TIME_POINTS];
	double xvals[TIME_POINTS];

	for(int p = 0; p < TIME_POINTS; ++p){
		norm[p] = 0.0;
		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j)
				norm[p] = norm[p] + full_solution[i][j][p];
		norm[p] = norm[p] * HX * HY;
	}

	// the norm should be unity, so we can subtract off 1 to get the error at each point in time
	for(int i = 0; i < TIME_POINTS; ++i)
		norm[i] = norm[i] - 1.0;

	// Now we plot the normalization array and print the norm loss rate
	printf("Normalization Loss Rate: %e\n", norm[TIME_POINTS-1]/TIME);

	// Create array of x values to plot then norm against
	for (int p = 0; p < TIME_POINTS; ++p)
		xvals[p] = p * Dt;
	
	char * commandsForGnuplot[] = {"set title \"Normalization Error vs. Time\"", "set xlabel \"Time\"", "set ylabel \"Normalization Error\"" , "plot 'norm.temp' with lines"};
	int num_commands = 4;
	plotFunction(commandsForGnuplot, num_commands, xvals, norm, time_points, "norm.temp");
}

// ----------------- Real Time Propagation Functions ---------------


double f(double imag_temp[SPX][SPY][4], double real_temp[SPX][SPY][4], int i, int j, int p){
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));

	return -.5 * (Dxx(imag_temp, i, j, p) + Dyy(imag_temp, i, j, p)) + V * imag_part + G * (pow(real_part, 2) + pow(imag_part, 2)) * imag_part;
}

double g(double real_temp[SPX][SPY][4], double imag_temp[SPX][SPY][4], int i, int j, int p){
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));

	return .5 * (Dxx(real_temp, i, j, p) + Dyy(real_temp, i, j, p)) - V * real_part - G * (pow(real_part, 2) + pow(imag_part, 2)) * real_part;
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = .000005; // This is the time step just used by the imaginary time propagation method

double PsiNorm(double Array[SPX][SPY][4]){
	// Integrates down the first x-y matrix assuming spatial square has area HX * HY

	double integral = 0.0;

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			integral = integral + pow(Array[i][j][0], 2);

	integral = HX * HY * integral;

	return integral;
}

double initialGuess(double x, double y){

	return x * (XLENGTH - x) * y * (YLENGTH - y);
}

void normalize(double Array[SPX][SPY][4]){
	// Takes an array with a function in the first column and normalizes that function back to 1

	double norm = sqrt( PsiNorm(Array) );

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			Array[i][j][0] = Array[i][j][0] / norm;
}

double fgs(double Psi_real[SPX][SPY][4], int i, int j, int p){
	// f is the function that gives the derivative for psi_real, which is why it is a function of imag_solution

		double real_part = Psi_real[i][j][p];
	
		return 1./2 * (Dxx(Psi_real, i, j, p) + Dyy(Psi_real, i, j, p))  - 1/2. * (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH/2.0) * WY, 2)) * real_part - G * pow(real_part , 2) * real_part;
}

void findGroundState(double real_solution[SPX][SPY][4], int iterations) {

	// Declare initial variables
	double real_temp[SPX][SPY][3];

	// Assign initial guess
	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			real_solution[i][j][0] = initialGuess(i * HX, j * HY);

	normalize(real_solution);

	for (int p = 1; p < iterations; ++p)
	{

		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j)
				real_solution[i][j][1] = real_solution[i][j][0] + Delta_t * fgs(real_solution, i, j, 0);

		// Move solution back an index so we can repeat the process
		for (int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
			{
				real_solution[i][j][0] = real_solution[i][j][1];
				real_solution[i][j][1] = 0.0;
			}


		// Now normalize the solution back to 1
		normalize(real_solution);
	}
}

int main(){

	// Declare initial variables
	static double real_solution[SPX][SPY][TIME_POINTS]; // This will be solution matrix where each column is a discrete point in time and each row a discrete point in space
	static double imag_solution[SPX][SPY][TIME_POINTS]; // Same as real solution but for the imaginary component of Psi
	double real_temp[SPX][SPY][4];
	double imag_temp[spx][spy][4];
	double real_temp8[SPX][SPY][4];
	double imag_temp8[spx][spy][4];
	double K[spx][spy][3];
	double L[spx][spy][3];
	double K8[spx][spy][3];
	double L8[spx][spy][3];
	double k_38, l_38;
	double k_3, l_3;


	static double full_solution[SPX][SPY][TIME_POINTS]; // This will have the psi square solution so we can plot the results
	

	// Initialize arrays to clear junk out of arrays
	for (int i = 0; i < spx; ++i){
		for (int j = 0; j < spy; ++j){

			for(int p = 0; p < time_points; ++p){
				real_solution[i][j][p] = 0;
				imag_solution[i][j][p] = 0;
				full_solution[i][j][p] = 0.0;
			}

			for (int p = 0; p < 3; ++p){
				real_temp[i][j][p] = 0;
				imag_temp[i][j][p] = 0;
				K[i][j][p] = 0;
				L[i][j][p] = 0;
			}
		}
	}

	// find the ground state
	findGroundState(real_temp, 10000);

	// Save ground state in a separate matrix for comparison and assign the ground state as initial condition to real soltuion
		// If we save the ground state, we can run real time propagation and determine how much the initial condition changes
		// If it changes very little then this would confirm that we have found the ground state
	double ground_state[SPX][SPY];
	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j){
			ground_state[i][j] = real_temp[i][j][0];
			real_solution[i][j][0] = real_temp[i][j][0];
		}

	double third = 1./3.;

	double ratio = 0.17293573625975;

	// run real time propagation algorithm
	for (int p = 1; p < time_points; ++p){

		// Load solution into first column of temp matrices
		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j){
				real_temp8[i][j][0] = real_solution[i][j][p - 1];
				imag_temp8[i][j][0] = imag_solution[i][j][p - 1];
				real_temp[i][j][0] = real_solution[i][j][p - 1];
				imag_temp[i][j][0] = imag_solution[i][j][p - 1];
		}

		// Forward Euler Step
		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K8[i][j][0] = f(imag_temp8, real_temp8, i, j, 0);
				L8[i][j][0] = g(real_temp8, imag_temp8, i, j, 0);
				K[i][j][0] = f(imag_temp, real_temp, i, j, 0);
				L[i][j][0] = g(real_temp, imag_temp, i, j, 0);
				real_temp8[i][j][1] = real_solution[i][j][p - 1] + third * Dt * K8[i][j][0]; 
				imag_temp8[i][j][1] = imag_solution[i][j][p - 1] + third * Dt * L8[i][j][0];
				real_temp[i][j][1] = real_solution[i][j][p - 1] + 0.5 * Dt * K[i][j][0]; 
				imag_temp[i][j][1] = imag_solution[i][j][p - 1] + 0.5 * Dt * L[i][j][0];
			}
		}

		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K8[i][j][1] = f(imag_temp8, real_temp8, i, j, 1);
				L8[i][j][1] = g(real_temp8, imag_temp8, i, j, 1);
				K[i][j][1] = f(imag_temp, real_temp, i, j, 1);
				L[i][j][1] = g(real_temp, imag_temp, i, j, 1);
				real_temp8[i][j][2] = real_solution[i][j][p - 1] + Dt * (-third * K8[i][j][0] + K8[i][j][1]);
				imag_temp8[i][j][2] = imag_solution[i][j][p - 1] + Dt * (-third * L8[i][j][0] + L8[i][j][1]);
				real_temp[i][j][2] = real_solution[i][j][p - 1] + 0.5 * Dt * K[i][j][1];
				imag_temp[i][j][2] = imag_solution[i][j][p - 1] + 0.5 * Dt * L[i][j][1];
			}
		}

		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K8[i][j][2] = f(imag_temp8, real_temp8, i, j, 2);
				L8[i][j][2] = g(real_temp8, imag_temp8, i, j, 2);
				K[i][j][2] = f(imag_temp, real_temp, i, j, 2);
				L[i][j][2] = g(real_temp, imag_temp, i, j, 2);
				real_temp8[i][j][3] = real_solution[i][j][p - 1] + Dt * (K8[i][j][0] - K8[i][j][1] + K8[i][j][2]);
				imag_temp8[i][j][3] = imag_solution[i][j][p - 1] + Dt * (L8[i][j][0] - L8[i][j][1] + L8[i][j][2]);
				real_temp[i][j][3] = real_solution[i][j][p - 1] + Dt * K[i][j][2];
				imag_temp[i][j][3] = imag_solution[i][j][p - 1] + Dt * L[i][j][2];
			}	
		}

		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				k_38 = f(imag_temp8, real_temp8, i, j, 3);
				l_38 = g(real_temp8, imag_temp8, i, j, 3);
				k_3 = f(imag_temp, real_temp, i, j, 3);
				l_3 = g(real_temp, imag_temp, i, j, 3);
				real_solution[i][j][p] = real_solution[i][j][p - 1] + \
				ratio * 1./8 * Dt * (K8[i][j][0] + 3 * K8[i][j][1] + 3 * K8[i][j][2] + k_38) + \
				(1 - ratio) * 1./6 * Dt * (K[i][j][0] + 2 * K[i][j][1] + 2 * K[i][j][2] + k_3);
				imag_solution[i][j][p] = imag_solution[i][j][p - 1] + \
				ratio * 1./8 * Dt * (L8[i][j][0] + 3 * L8[i][j][1] + 3 * L8[i][j][2] + l_38) + \
				(1 - ratio) * 1./6 * Dt * (L[i][j][0] + 2 * L[i][j][1] + 2 * L[i][j][2] + l_3);
			}
		}
	}

	// Find the psi squared solution at the end time so we can plot the results
	for (int i = 0; i < spx; ++i)
		for (int j = 0; j < spy; ++j)
			for(int p = 0; p < time_points; ++p)
				full_solution[i][j][p] = pow(real_solution[i][j][p], 2) + pow(imag_solution[i][j][p], 2);


	// Now lets plot our results

	// Set up commands for gnuplot
	char title[100];
	sprintf(title, "set title \"Schrodinger Equation Harmonic Trap (Imag Time) (Runge Kutta)");
	char * commandsForGnuplot[] = { title, "set xlabel \"X-Axis\"", "show xlabel", "set ylabel \"Y-Axis\"", "show ylabel" , "splot 'sol.temp' with lines"};
	int num_commands = 6;


	// plotSurface(commandsForGnuplot, num_commands, full_solution, spx, spy, TIME_POINTS - 1, "sol.temp");

	printf("Space Step: %e \nTime Step: %e %e \nRatio: %f\n", HX, Dt, 0, HX / Dt);		

	// We also want to compare how much the ground state has changed

	// First form a difference matrix
	double difference_matrix[SPX][SPY][TIME_POINTS];
	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			difference_matrix[i][j][0] = pow(ground_state[i][j], 2) - full_solution[i][j][time_points - 1];

	char * commandsForGnuplot2[] = {"set autoscale", "set title \"Difference in Ground State after Evolution\"", "set xlabel \"X-Axis\"", "show xlabel", "set ylabel \"Y-Axis\"", "show ylabel" , "splot 'difference.temp' with lines"};
	num_commands = 7;

	// plotSurface(commandsForGnuplot2, num_commands, difference_matrix, spx, spy, 0, "difference.temp");

	plotNormalization(full_solution);

	return 0;
}