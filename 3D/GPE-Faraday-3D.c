/*
This program numerically finds the thomas fermi ground state, solves the nonlinear schrodinger equation with harmonic potential in 2D using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant W is the harmonic potential constant omega

RESULTS:
1) 80 length units along x and y direction --> g.s. has enough room along these directions

TODO:
1) Create stepTwo function for 2D -- postpone
2) why are there jumps on the g.s.?

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 0.00000250
#define XLENGTH 90.0 
#define YLENGTH 90.0 
#define ZLENGTH 105.0
#define TIME_POINTS 4 //number of time points
#define SPX 50 //600
#define SPY 50 //20
#define SPZ 300
#define NOISE_VOLUME 0.0

#define Dxx(array, x, y, z, pee) ( (-1. * array[mod(x + 2, SPX)][y][z][pee] + 16.* array[mod(x + 1, SPX)][y][z][pee] - \
		30. * array[mod(x , SPX)][y][z][pee] + 16. * array[mod(x - 1, SPX)][y][z][pee] +  -1 * array[mod(x - 2, SPX)][y][z][pee]) / (12. * pow(HX, 2)) )
;

#define Dyy(array, x, y, z, pee) ( (-1. * array[x][mod(y + 2, SPY)][z][pee] + 16.* array[x][mod(y + 1, SPY)][z][pee] - \
		30. * array[x][mod(y , SPY)][z][pee] + 16. * array[x][mod(y - 1, SPY)][z][pee] +  -1 * array[x][mod(y - 2, SPY)][z][pee]) / (12. * pow(HY, 2)) )
;

#define Dzz(array, x, y ,z , pee) ( (-1. * array[x][y][mod(z + 2, SPZ)][pee] + 16.* array[x][y][mod(z + 1, SPZ)][pee] - \
		30. * array[x][y][mod(z , SPZ)][pee] + 16. * array[x][y][mod(z - 1, SPZ)][pee] +  -1 * array[x][y][mod(z - 2, SPZ)][pee]) / (12. * pow(HZ, 2)) )
;

const double HX = XLENGTH / (SPX);
const double HY = YLENGTH / (SPY);
const double HZ = ZLENGTH / (SPZ);
const double WX = 1.;
const double WY = 1.;
const double WZ = 7./476;
const double G = 1000.;
const double OMEGA = 2.;
const double EPS = 0.0;
const double WAVENUMBER_INPUT = 1.1;

// How many time steps do we want
const int time_points = TIME_POINTS;

// How many spatial points are there
const int spx = SPX;
const int spy = SPY;
const int spz = SPZ;

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

void plotTimeDependence(double solution1D[SPZ][TIME_POINTS], int spz, struct plot_settings plot){

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

	struct plot_information3D info = { .num_commands = num_commands, \
		.output_file = "1Dsol.temp", .x_length = XLENGTH, .y_length = YLENGTH, .z_length = ZLENGTH, .T = TIME};

	plotSurface3DReduced(commands, info , spz, TIME_POINTS, solution1D);
}

void plotNormalization(int arg_time_points, double solution1D[SPZ][arg_time_points], double T){
	/*
	* This function takes the squared solution matrix, finds the normalization at each point in time and plots the error (normalization should equal unity)
	*/

	double norm[arg_time_points];
	double xvals[arg_time_points];
	double Dt = T/arg_time_points;

	for(int p = 0; p < arg_time_points; ++p){
		norm[p] = 0.0;
		for(int i = 0; i < SPZ; ++i)
			norm[p] = norm[p] + solution1D[i][p];
		norm[p] = norm[p] * HZ;
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

double contraction(double solution[SPX][SPY][SPZ][4], int z_index ){
	/*
	This function will take a 3D matrix and contract across the y dimension at a given z- index (integer z_index) -- assumed to be on index 0 for 4th-axis
	*/

	double integral_holder = 0.0;
	// Integrate over the x and y-axis

	for(int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			integral_holder = integral_holder + solution[i][j][z_index][0];	

	return integral_holder * HY * HX;
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

double f(double imag_temp[SPX][SPY][SPZ][4], double real_temp[SPX][SPY][SPZ][4], int i, int j, int k, int p, double t){
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][k][p];
	double real_part = real_temp[i][j][k][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2) + pow((k * HZ - ZLENGTH / 2.0) * WZ, 2));
	double new_G = G * (1. +  EPS * sin(OMEGA * t));

	return -.5 * (Dxx(imag_temp, i, j, k, p) + Dyy(imag_temp, i, j, k, p) + Dzz(imag_temp, i ,j ,k ,p)) + V * imag_part + new_G * (pow(real_part, 2) + pow(imag_part, 2)) * imag_part;
}

double g(double real_temp[SPX][SPY][SPZ][4], double imag_temp[SPX][SPY][SPZ][4], int i, int j, int k, int p, double t){
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][k][p];
	double real_part = real_temp[i][j][k][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2) + pow((k * HZ - ZLENGTH / 2.0) * WZ, 2));
	double new_G = G * (1. +  EPS * sin(OMEGA * t));


	return .5 * (Dxx(real_temp, i, j, k, p) + Dyy(real_temp, i, j, k, p) + Dzz(imag_temp, i ,j ,k ,p)) - V * real_part - new_G * (pow(real_part, 2) + pow(imag_part, 2)) * real_part;
}

void addNoise(double initialCondition[SPX][SPY][SPZ][4], double waveNumber){
	/*
	This function adds noise to the initial condition matrix
	*/

	// First we need to find a wavelength of noise that has an integer number of wavelengths within the space 
	//  	of the trap and also is close to the assigned wavenumber in the function input
	double num = waveNumber * ZLENGTH / PI / 2.;
	int n = (int) floor (num);

	if (mod(n, 2) != 0)
		n = (int) ceil(num);

	double kp_eff = n * 2 * PI / ZLENGTH;

	// Noise should be small relative to homogeneous solution
	double noise_volume = NOISE_VOLUME;

	// Loop through initial condition and add noise
	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				initialCondition[i][j][k][0] = initialCondition[i][j][k][0] * (1 + noise_volume * sin(kp_eff * k * HZ));
}

void realTimeProp(double initialCondition[SPX][SPY][SPZ][4], double T, int arg_time_points, double real_solution[][SPY][SPZ][4], double imag_solution[][SPY][SPZ][4], struct plot_settings plot){
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

	double real_temp[SPX][SPY][SPZ][4]; // This matrix and the one below will be used to store the state of the system at intermediary time steps -- time steps 
	double imag_temp[SPX][SPY][SPZ][4]; // that are need to do RK4 but are not essential otherwise and can therefore be discarded after each iteration
	double K[SPX][SPY][SPZ][3];
	double L[SPX][SPY][SPZ][3];
	double k_3, l_3;

	double full_solution[SPX][SPY][SPZ][4]; // This matrix stores the psi squared solution at each time slice to make contraction more convenient
	double solution1D[SPZ][arg_time_points]; //This matrix contracts the full solution along the radial 'skinny' dimension because it is not essential to observe dynamics
	
	// Initialize arrays to clear junk out of arrays
	for (int i = 0; i < spx; ++i){
		for (int j = 0; j < spy; ++j){
			for (int k = 0; k < spz; ++k){
				for (int p = 0; p < 3; ++p){
					real_solution[i][j][k][p] = 0;
					imag_solution[i][j][k][p] = 0;
					full_solution[i][j][k][p] = 0.0;

					real_temp[i][j][k][p] = 0;
					imag_temp[i][j][k][p] = 0;
					K[i][j][k][p] = 0;
					L[i][j][k][p] = 0;
				}
			}
		}
	}

	// Add the noise to the ground state
	addNoise(initialCondition, WAVENUMBER_INPUT);

	// Assign initial conditions
	for (int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k){
				real_solution[i][j][k][0] = initialCondition[i][j][k][0];
				full_solution[i][j][k][0] = pow(initialCondition[i][j][k][0], 2);
			}
	
	for(int i = 0; i < SPZ; ++i)
		solution1D[i][0] = contraction(full_solution, i);

	double t;
	// Real Time Propagation
	for (int p = 1; p < arg_time_points; ++p)
	{
		// Load solution into first column of temp matrices
		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j)
				for(int k = 0; k < spz; ++k){
					real_temp[i][j][k][0] = real_solution[i][j][k][0];
					imag_temp[i][j][k][0] = imag_solution[i][j][k][0];
				}


		t = Dt * (p - 1); //Forward Euler step
		for (int i = 0; i < spx; ++i)
			for (int j = 0; j < spy; ++j)
				for (int k = 0; k < spz; ++k){
					K[i][j][k][0] = f(imag_temp, real_temp, i, j, k, 0, t);
					L[i][j][k][0] = g(real_temp, imag_temp, i, j, k, 0, t);
					real_temp[i][j][k][1] = real_solution[i][j][k][0] + .5 * Dt * K[i][j][k][0]; 
					imag_temp[i][j][k][1] = imag_solution[i][j][k][0] + .5 * Dt * L[i][j][k][0];
				}

		t = t + .5 * Dt; //Add half a time step 
		for (int i = 0; i < spx; ++i)
			for (int j = 0; j < spy; ++j)
				for(int k = 0; k < spz; ++k){
					K[i][j][k][1] = f(imag_temp, real_temp, i, j, k, 1, t);
					L[i][j][k][1] = g(real_temp, imag_temp, i, j, k, 1, t);
					real_temp[i][j][k][2] = real_solution[i][j][k][0] + .5 * Dt * K[i][j][k][1];
					imag_temp[i][j][k][2] = imag_solution[i][j][k][0] + .5 * Dt * L[i][j][k][1];
				}


		// t does not change for this step
		for (int i = 0; i < spx; ++i)
			for (int j = 0; j < spy; ++j)
				for(int k = 0; k < spz; ++k){
					K[i][j][k][2] = f(imag_temp, real_temp, i, j, k, 2, t);
					L[i][j][k][2] = g(real_temp, imag_temp, i, j, k, 2, t);
					real_temp[i][j][k][3] = real_solution[i][j][k][0] + Dt * K[i][j][k][2];
					imag_temp[i][j][k][3] = imag_solution[i][j][k][0] + Dt * L[i][j][k][2];
				}	

		t = Dt * p; //Add full step for Backward Euler step
		for (int i = 0; i < spx; ++i)
			for (int j = 0; j < spy; ++j)
				for(int k = 0; k < spz; ++k){
					k_3 = f(imag_temp, real_temp, i, j, k, 3, t);
					l_3 = g(real_temp, imag_temp, i, j, k, 3, t);
					real_solution[i][j][k][1] = real_solution[i][j][k][0] + 1./6 * Dt * (K[i][j][k][0] + 2 * K[i][j][k][1] + 2 * K[i][j][k][2] + k_3);
					imag_solution[i][j][k][1] = imag_solution[i][j][k][0] + 1./6 * Dt * (L[i][j][k][0] + 2 * L[i][j][k][1] + 2 * L[i][j][k][2] + l_3);
				}

		// Add new iteration to full solution
		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j)
				for(int k = 0; k < SPZ; ++k)
					full_solution[i][j][k][0] = pow(real_solution[i][j][k][1], 2) + pow(imag_solution[i][j][k][1], 2);

		// Contract full solution and save it to solution1D
		for(int i = 0; i < SPZ; ++i)
			solution1D[i][p] = contraction(full_solution, i);

		// Move new iteration in real/imag solution back an index to restart process
		for(int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
				for(int k = 0; k < SPZ; ++k)
				{
					real_solution[i][j][k][0] = real_solution[i][j][k][1];
					imag_solution[i][j][k][0] = imag_solution[i][j][k][1];
				}


	}

	if (plot.plot3D == 1){

		// Now lets plot our results
		plotTimeDependence(solution1D, SPZ, plot);
	}
	if (plot.plot_normalization == 1){

		plotNormalization(arg_time_points, solution1D, T);
	}
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = .000025; // This is the time step just used by the imaginary time propagation method

double PsiNorm(double Array[SPX][SPY][SPZ][4]){
	// Integrates down the first x-y matrix assuming spatial square has area HX * HY

	double integral = 0.0;

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				integral = integral + pow(Array[i][j][k][0], 2);

	integral = HX * HY * HZ * integral;

	return integral;
}

double initialGuess(double x, double y, double z){

	return x * (XLENGTH - x) * y * (YLENGTH - y) * z * (ZLENGTH - z);
}

void normalize(double Array[SPX][SPY][SPZ][4]){
	// Takes an array with a function in the first column and normalizes that function back to 1

	double norm = sqrt( PsiNorm(Array) );

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				Array[i][j][k][0] = Array[i][j][k][0] / norm;
}

double fgs(double Psi_real[SPX][SPY][SPZ][4], int i, int j, int k, int p){
	// f is the function that gives the derivative for psi_real, which is why it is a function of imag_solution

		double real_part = Psi_real[i][j][k][p];
		double V = pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH/2.0) * WY, 2) + pow((k * HZ - ZLENGTH/2.0) * WZ, 2);
	
		return 1./2 * (Dxx(Psi_real, i, j, k, p) + Dyy(Psi_real, i, j, k, p) + Dzz(Psi_real, i, j, k, p))  - 1/2. * V * real_part - G * pow(real_part , 2) * real_part;
}

void findGroundState(double real_solution[SPX][SPY][SPZ][4], int iterations) {

	// Assign initial guess
	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				real_solution[i][j][k][0] = initialGuess(i * HX, j * HY, k * HZ);

	normalize(real_solution);

	for (int p = 1; p < iterations; ++p){

		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j)
				for(int k = 0; k < SPZ; ++k)
					real_solution[i][j][k][1] = real_solution[i][j][k][0] + Delta_t * fgs(real_solution, i, j, k, 0);

		for(int i = 0; i < spx; ++i)
			for(int j = 0; j < spy; ++j)
				for(int k = 0; k < SPZ; ++k)
					real_solution[i][j][k][2] = real_solution[i][j][k][1] + Delta_t * fgs(real_solution, i, j, k, 1);

		// Move solution back an index so we can repeat the process
		for (int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
				for(int k = 0; k < SPZ; ++k)
				{
					real_solution[i][j][k][0] = real_solution[i][j][k][2];
					real_solution[i][j][k][1] = 0.0;
					real_solution[i][j][k][2] = 0.0;
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
	static double real_solution[SPX][SPY][SPZ][4]; // This will be the solution matrix -- it will only hold the state of the system at one time slice
	static double imag_solution[SPX][SPY][SPZ][4]; // Same as real solution but for the imaginary component of Psi
	double initialCondition[SPX][SPY][SPZ][4];

	struct plot_settings plot_solution = {.plot3D = 1, .title = "Real Time Solution", .plot_normalization = 1};


	// find the ground state
	findGroundState(initialCondition, 400);

	for(int i = 0; i < SPX; ++i)
		printf("%e ", initialCondition[i][25][150][0]);

	// Run RTP
	realTimeProp(initialCondition, TIME, TIME_POINTS, real_solution, imag_solution, plot_solution);


	return 0;
}