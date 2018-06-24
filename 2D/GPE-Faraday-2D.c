/*
This program numerically finds the thomas fermi ground state, solves the nonlinear schrodinger equation with harmonic potential in 2D using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant W is the harmonic potential constant omega

TODO:
1) Create stepTwo function for 2D -- postpone
2) Fix memory problems -- postpone


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 15.0
#define XLENGTH 300.0 //4.0
#define YLENGTH 10.0 //0.1
#define TIME_POINTS 4000 //number of time points
#define SPX 1000 //600
#define SPY 30 //20
#define NOISE_VOLUME 0.06

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
const double EPS = 0.1;
const double WAVENUMBER_INPUT = 1.5;

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

void plotTimeDependence(double solution[SPX][SPY][TIME_POINTS], int spx, int spy, struct plot_settings plot){

	double solution1D[SPX][TIME_POINTS];

	double integral_holder = 0.0;
	// Integrate over the y-axis
	for(int p = 0; p < TIME_POINTS; ++p)	
		for(int i = 0; i < SPX; ++i){
			
			integral_holder = 0.0;

			for (int j = 0; j < SPY; ++j)
				integral_holder = integral_holder + solution[i][j][p] * HY;	

			solution1D[i][p] = integral_holder;
		}

	// Set up commands for gnuplot
	char plotCommand[100];
	char title[100];
	// sprintf(plotCommand, "plot [0:%d] [-2:2] 'sol.temp' ", TIME);
	strcpy(plotCommand,"splot '1Dsol.temp' with lines");
	sprintf(title, "set title \"");
	strcat(title, plot.title);
	strcat(title, "\"");
	char * commands[15] = { title, "set xlabel \"Position\"", "show xlabel", "set ylabel \"Time\"", "show ylabel" , plotCommand};
	int num_commands = 6;

	struct plot_information info = { .num_commands = num_commands, \
		.output_file = "1Dsol.temp", .x_length = XLENGTH, .y_length = YLENGTH, .T = TIME};

	plotSurface2DReduced(commands, info , spx, TIME_POINTS, solution1D);
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
	for(int i = 1; i < arg_time_points; ++i)
		norm[i] = norm[i] - norm[0];

	norm[0] = norm[0] - norm[0];
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

		plotTimeDependence(full_solution, SPX, SPY, plot);

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

	printf("test\n");
	// Declare initial variables
	static double real_solution[SPX][SPY][TIME_POINTS]; // This will be solution matrix where each column is a discrete point in time and each row a discrete point in space
	static double imag_solution[SPX][SPY][TIME_POINTS]; // Same as real solution but for the imaginary component of Psi
	double initialCondition[SPX][SPY][4];

	struct plot_settings plot_solution = {.plot3D = 1, .title = "Real Time Solution", .plot_normalization = 1};


	// find the ground state
	findGroundState(initialCondition, 800);

	// Run RTP
	realTimeProp(initialCondition, TIME, TIME_POINTS, real_solution, imag_solution, plot_solution);


	return 0;
}