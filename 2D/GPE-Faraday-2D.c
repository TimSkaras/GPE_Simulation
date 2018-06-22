/*
This program numerically finds the thomas fermi ground state, solves the nonlinear schrodinger equation with harmonic potential in 2D using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant W is the harmonic potential constant omega

TODO:
1) Create stepTwo function for 2D -- postpone
2) Fix memory problems -- postpone
3) Implement addNoise function -- done
4) Modify RTP to support time dependent equations -- done
5) Add integration to plotter so that time dependence can be included


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 10.0
#define XLENGTH 300.0 //4.0
#define YLENGTH 10.0 //0.1
#define TIME_POINTS 1500 	//number of time points
#define SPX 800 //600
#define SPY 20 //20
#define NOISE_VOLUME 0.02

#define Dxx(array, x, y, pee) ( (-1. * array[mod(x + 2, SPX)][y][pee] + 16.* array[mod(x + 1, SPX)][y][pee] - \
		30. * array[mod(x , SPX)][y][pee] + 16. * array[mod(x - 1, SPX)][y][pee] +  -1 * array[mod(x - 2, SPX)][y][pee]) / (12. * pow(HX, 2)) )
;

#define Dyy(array, x, y, pee) ( (-1. * array[x][mod(y + 2, SPY)][pee] + 16.* array[x][mod(y + 1, SPY)][pee] - \
		30. * array[x][mod(y , SPY)][pee] + 16. * array[x][mod(y - 1, SPY)][pee] +  -1 * array[x][mod(y - 2, SPY)][pee]) / (12. * pow(HY, 2)) )
;

const double HX = XLENGTH / (SPX);
const double HY = YLENGTH / (SPY);
const double WX = 7./476;
const double WY = 1.;
const double G = 1000.;
const double OMEGA = 2.;
const double EPS = 0.;
const double WAVENUMBER_INPUT = .11;

// How many time steps do we want
const int time_points = TIME_POINTS;

// How many spatial points are there
const int spx = SPX;
const int spy = SPY;

struct plot_settings{
	/*
	*	This struct is used to pass the realTimeProp function settings for the plotting capabalities
	*
	*	MEMBERS:
	*	plot3D -- integer value either 0 or 1, 1 turns plotting feature on and 0 turns it off
	*	title -- string that tells the plotter what to display for the title of the solution plot
	*	plot_normalization -- tells plotter whether to plot the normalization curve as well
	*
	*/
	
	int plot3D;
	char * title;
	int plot_normalization;
};

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
	int index_x = (int) ceil(spx/500.);
	int new_x = (int) floor(spx /  ceil(spx/500.));
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

void plotNormalization(double full_solution[SPX][SPY][TIME_POINTS], int arg_time_points, double T){
	/*
	* This function takes the squared solution matrix, finds the normalization at each point in time and plots the error (normalization should equal unity)
	*/

	double norm[arg_time_points];
	double xvals[arg_time_points];
	double Dt = T/arg_time_points;

	for(int p = 0; p < arg_time_points; ++p){
		norm[p] = 0.0;
		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j)
				norm[p] = norm[p] + full_solution[i][j][p];
		norm[p] = norm[p] * HX * HY;
	}

	// the norm should be unity, so we can subtract off 1 to get the error at each point in time
	for(int i = 0; i < arg_time_points; ++i)
		norm[i] = norm[i] - 1.0;

	// Now we plot the normalization array

	// Create array of x values to plot then norm against
	for (int p = 0; p < arg_time_points; ++p)
		xvals[p] = p * Dt;
	
	char * commandsForGnuplot[] = {"set title \"Normalization Error vs. Time\"", "set xlabel \"Time\"", "set ylabel \"Normalization Error\"" , "plot 'norm.temp' with lines"};
	int num_commands = 4;
	plotFunction(commandsForGnuplot, num_commands, xvals, norm, time_points, "norm.temp");
}

void loadMatrix(double matrix[SPX][SPY][4], char * filename){
	/*
	Takes matrix data in filename and loads it into matrix
	*/

	FILE * f = fopen(filename, "r");

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			matrix[i][j][0] = 0;

	fclose(f);
}

// ----------------- Real Time Propagation Functions ---------------

double f(double imag_temp[SPX][SPY][4], double real_temp[SPX][SPY][4], int i, int j, int p, double t){
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));
	double new_G = G * (1. +  EPS * sin(OMEGA * t));

	return -.5 * (Dxx(imag_temp, i, j, p) + Dyy(imag_temp, i, j, p)) + V * imag_part + new_G * (pow(real_part, 2) + pow(imag_part, 2)) * imag_part;
}

double g(double real_temp[SPX][SPY][4], double imag_temp[SPX][SPY][4], int i, int j, int p, double t){
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));
	double new_G = G * (1. +  EPS * sin(OMEGA * t));


	return .5 * (Dxx(real_temp, i, j, p) + Dyy(real_temp, i, j, p)) - V * real_part - new_G * (pow(real_part, 2) + pow(imag_part, 2)) * real_part;
}

void addNoise(double initialCondition[SPX][SPY][4], double waveNumber){
	/*
	This function adds noise to the initial condition matrix
	*/

	// First we need to find a wavelength of noise that has an integer number of wavelengths within the space 
	//  	of the trap and also is close to the assigned wavenumber in the function input
	double num = waveNumber * XLENGTH / PI / 2.;
	int n = (int) floor (num);

	if (mod(n, 2) != 0)
		n = (int) ceil(num);

	double kp_eff = n * 2 * PI / XLENGTH;

	// Noise should be small relative to homogeneous solution
	double noise_volume = NOISE_VOLUME;

	// Loop through initial condition and add noise
	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			initialCondition[i][j][0] = initialCondition[i][j][0] * (1 + noise_volume * sin(kp_eff * i * HX));
}

void realTimeProp(double initialCondition[SPX][SPY][4], double T, int arg_time_points, double real_solution[][SPY][arg_time_points], double imag_solution[][SPY][arg_time_points], struct plot_settings plot){
	/*
	*	This function takes uses the array initialCondition as an initial condition and propagates that profile forward 
	*	in real time using RK4. The function outputs the solution matrix, which has dimensions specified by space_points
	*	constant and the arg_time_points argument parameter. This solver assumes periodic boundary conditions
	*
	*	INPUTS:
	*	initialCondition -- array of length SPACE_POINTS
	*	T -- double indicating the total desired length of simulated time
	*	arg_time_points -- integer indicating the desired number of discretized points in time
	*	real_solution -- empty two dimensional array of doubles that will be used to store the real part of the solution
	*	imag_solution -- empty two dimensional array of doubles that will be used to store the imaginary part of the solution
	*	plot -- struct carrying information to tell the plotter what to do. See declaration of struct for more information
	*
	*	OUTPUT:
	*	Plots result if plot_on is turned on and returns solution solution matrix with dimensions SPACE_POINTS x arg_time_points
	*
	*/

	double Dt = T/arg_time_points;

	double real_temp[SPX][SPY][4];
	double imag_temp[SPX][SPY][4];
	double K[SPX][SPY][3];
	double L[SPX][SPY][3];
	double k_3, l_3;

	double full_solution[SPX][SPY][arg_time_points]; // This will have the psi square solution so we can plot the results
	
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

	// Add the noise to the ground state
	addNoise(initialCondition, WAVENUMBER_INPUT);

	// Assign initial conditions
	for (int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			real_solution[i][j][0] = initialCondition[i][j][0];


	double t;
	// Real Time Propagation
	for (int p = 1; p < arg_time_points; ++p)
	{
		// Load solution into first column of temp matrices
		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j){
				real_temp[i][j][0] = real_solution[i][j][p - 1];
				imag_temp[i][j][0] = imag_solution[i][j][p - 1];
		}


		t = Dt * (p - 1); //Forward Euler step
		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K[i][j][0] = f(imag_temp, real_temp, i, j, 0, t);
				L[i][j][0] = g(real_temp, imag_temp, i, j, 0, t);
				real_temp[i][j][1] = real_solution[i][j][p - 1] + .5 * Dt * K[i][j][0]; 
				imag_temp[i][j][1] = imag_solution[i][j][p - 1] + .5 * Dt * L[i][j][0];
			}
		}

		t = t + .5 * Dt; //Add half a time step 
		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K[i][j][1] = f(imag_temp, real_temp, i, j, 1, t);
				L[i][j][1] = g(real_temp, imag_temp, i, j, 1, t);
				real_temp[i][j][2] = real_solution[i][j][p - 1] + .5 * Dt * K[i][j][1];
				imag_temp[i][j][2] = imag_solution[i][j][p - 1] + .5 * Dt * L[i][j][1];
			}
		}

		// t does not change for this step
		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				K[i][j][2] = f(imag_temp, real_temp, i, j, 2, t);
				L[i][j][2] = g(real_temp, imag_temp, i, j, 2, t);
				real_temp[i][j][3] = real_solution[i][j][p - 1] + Dt * K[i][j][2];
				imag_temp[i][j][3] = imag_solution[i][j][p - 1] + Dt * L[i][j][2];
			}	
		}

		t = Dt * p; //Add full step for Backward Euler step
		for (int i = 0; i < spx; ++i){
			for (int j = 0; j < spy; ++j){

				k_3 = f(imag_temp, real_temp, i, j, 3, t);
				l_3 = g(real_temp, imag_temp, i, j, 3, t);
				real_solution[i][j][p] = real_solution[i][j][p - 1] + 1./6 * Dt * (K[i][j][0] + 2 * K[i][j][1] + 2 * K[i][j][2] + k_3);
				imag_solution[i][j][p] = imag_solution[i][j][p - 1] + 1./6 * Dt * (L[i][j][0] + 2 * L[i][j][1] + 2 * L[i][j][2] + l_3);
			}
		}
	}

	if (plot.plot3D == 1){
		// Find the psi squared solution so we can plot the results
		for (int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
				for(int k = 0; k < arg_time_points; ++k)
					full_solution[i][j][k] = pow(real_solution[i][j][k], 2) + pow(imag_solution[i][j][k], 2);

		// Now lets plot our results

		// Set up commands for gnuplot
		char plotCommand[100];
		char title[100];
		strcpy(plotCommand,"splot 'sol.temp' with lines");
		sprintf(title, "set title \"");
		strcat(title, plot.title);
		strcat(title, "\"");
		char * commandsForGnuplot[] = { title, "set xlabel \"X-Axis\"", "show xlabel", "set ylabel \"Y-Axis\"", "show ylabel" , plotCommand};
		int num_commands = 6;


		plotSurface(commandsForGnuplot, num_commands, full_solution, spx, spy, arg_time_points - 1,  "sol.temp");


		// We also want to compare how much the ground state has changed

		// First form a difference matrix
		double difference_matrix[SPX][SPY][TIME_POINTS];
		for (int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
				difference_matrix[i][j][0] = pow(initialCondition[i][j][0], 2) - full_solution[i][j][time_points - 1];

		char * commandsForGnuplot2[] = {"set autoscale", "set title \"Difference in Ground State after Evolution\"", "set xlabel \"X-Axis\"", "show xlabel", "set ylabel \"Y-Axis\"", "show ylabel" , "splot 'difference.temp' with lines"};
		num_commands = 7;

		plotSurface(commandsForGnuplot2, num_commands, difference_matrix, spx, spy, 0, "difference.temp");

	}
	if (plot.plot_normalization == 1){

		plotNormalization(full_solution,arg_time_points, T);
	}
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = .0002; // This is the time step just used by the imaginary time propagation method

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

	// Write the g.s. to a separate file
	FILE * f = fopen("groundState.txt","w");

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j){

			if (j < SPY - 1)
				fprintf(f, "%e ", real_solution[i][j][0]);
			else
				fprintf(f, "%e\n", real_solution[i][j][0]);
		}
	fclose(f);
}

int main(){

	// Declare initial variables
	static double real_solution[SPX][SPY][TIME_POINTS]; // This will be solution matrix where each column is a discrete point in time and each row a discrete point in space
	static double imag_solution[SPX][SPY][TIME_POINTS]; // Same as real solution but for the imaginary component of Psi
	double initialCondition[SPX][SPY][4];

	struct plot_settings plot_solution = {.plot3D = 1, .title = "Real Time Solution.", .plot_normalization = 1};


	// find the ground state
	findGroundState(initialCondition, 80000);

	// Run RTP
	realTimeProp(initialCondition, TIME, TIME_POINTS, real_solution, imag_solution, plot_solution);


	return 0;
}