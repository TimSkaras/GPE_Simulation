/*
This file has the same functionality as GPE-Faraday-3D.c but I have changed several operations, functions, or structural elements to try and make it faster.
Code that has been optimized is often less maintainable or less readable, hence the intention of keeping the original code and the optimized code separate.

This program numerically finds the thomas fermi ground state, solves the nonlinear schrodinger equation with harmonic potential in 2D using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant W is the harmonic potential constant omega

I have implemented parallelization capabilities for RTP. To make use of it, you must use the following line to compile and run the code

gcc -Wl,--stack,1073741824 GPE-Faraday-3D.c -fopenmp -o GPE-Faraday-3D.exe %% GPE-Faraday-3D.exe

the "Wl" flag increases the stack size and the -fopenmp flag turns on parallel funcitonality

You can set the number of cores to use in the preprocessor constants 

RESULTS:
1) 80 length units along x and y direction --> g.s. has enough room along these directions
TODO:
1) Create stepTwo function for 2D -- postpone

160 -- 106 seconds


*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#include <sys/time.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 0.00250
#define XLENGTH 12.0 
#define YLENGTH 12.0 
#define ZLENGTH 210.0
#define TIME_POINTS 15 //number of time points
#define SPX 32 //600
#define SPY 32 //20
#define SPZ 512
#define NOISE_VOLUME 0.06
#define NUM_THREADS 2

const double xdivisor = 1./ (12. * pow(XLENGTH / (SPX), 2));
const double ydivisor = 1./ (12. * pow(YLENGTH / (SPY), 2));
const double zdivisor = 1./ (12. * pow(ZLENGTH / (SPZ), 2));
// Full Fourth Order Spatial Derivative
#define Dxx(array, x, y, z, pee) ( (-1. * array[mod(x + 2, SPX)][y][z][pee] + 16.* array[mod(x + 1, SPX)][y][z][pee] - \
		30. * array[mod(x , SPX)][y][z][pee] + 16. * array[mod(x - 1, SPX)][y][z][pee] +  -1 * array[mod(x - 2, SPX)][y][z][pee]) * xdivisor)
;

#define Dyy(array, x, y, z, pee) ( (-1. * array[x][mod(y + 2, SPY)][z][pee] + 16.* array[x][mod(y + 1, SPY)][z][pee] - \
		30. * array[x][mod(y , SPY)][z][pee] + 16. * array[x][mod(y - 1, SPY)][z][pee] +  -1 * array[x][mod(y - 2, SPY)][z][pee]) * ydivisor )
;

#define Dzz(array, x, y ,z , pee) ( (-1. * array[x][y][mod(z + 2, SPZ)][pee] + 16.* array[x][y][mod(z + 1, SPZ)][pee] - \
		30. * array[x][y][mod(z , SPZ)][pee] + 16. * array[x][y][mod(z - 1, SPZ)][pee] +  -1 * array[x][y][mod(z - 2, SPZ)][pee]) * zdivisor )
;


/*const double xdivisor = 1./ (2. * pow(XLENGTH / (SPX), 2));
const double ydivisor = 1./ (2. * pow(YLENGTH / (SPY), 2));
const double zdivisor = 1./ (2. * pow(ZLENGTH / (SPZ), 2));
// Full Second Order Spatial Derivative 
#define Dxx(array, x, y, z, pee) ( (array[mod(x + 1, SPX)][y][z][pee] - 2. * array[mod(x , SPX)][y][z][pee] + array[mod(x - 1, SPX)][y][z][pee]) * xdivisor)
;
#define Dyy(array, x, y, z, pee) ( (array[x][mod(y + 1, SPY)][z][pee] - 2. * array[x][mod(y , SPY)][z][pee] +  array[x][mod(y - 1, SPY)][z][pee]) * ydivisor)
;
#define Dzz(array, x, y ,z , pee) ( ( array[x][y][mod(z + 1, SPZ)][pee] - 2.* array[x][y][mod(z , SPZ)][pee] + array[x][y][mod(z - 1, SPZ)][pee] ) * zdivisor)
;*/

const double HX = XLENGTH / (SPX);
const double HY = YLENGTH / (SPY);
const double HZ = ZLENGTH / (SPZ);
const double WX = 1.; // 1
const double WY = 1.; // 1
const double WZ = 7./476; // 7./476
const double G = 1000.;
const double OMEGA = 2.;
const double EPS = 0.15;
const double WAVENUMBER_INPUT = 1.1;
const double T_MOD = 5 * PI; // Amount of time the scattering length is modulated
const int RED_COEFF = 1;

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
	char * output_file;
	int plot_normalization;
};

// ------------------ Auxiliary Functions ----------------

int mod(int a, int b){
	// This is a standard modulus function -- always outputs a positive number
	// E.g., mod(-1, 3) = 2

    int r = a % b;
    return r < 0 ? r + b : r;
}

void plotTimeDependence(int sp, int arg_time_points, double solution1D[sp][arg_time_points], struct plot_settings plot){

	// Set up commands for gnuplot
	char plotCommand[100];
	char title[100];
	// sprintf(plotCommand, "plot [0:%d] [-2:2] 'sol.temp' ", TIME);
	sprintf(plotCommand,"splot '%s' with lines", plot.output_file);
	sprintf(title, "set title \"");
	strcat(title, plot.title);
	strcat(title, "\"");
	char * commands[15] = { title, "set xlabel \"Position\"", "show xlabel", "set ylabel \"Time\"", "show ylabel" , plotCommand};
	int num_commands = 6;
	
	struct plot_information3D info = { .num_commands = num_commands, \
		.output_file = plot.output_file, .length = (sp == SPX) ? XLENGTH:ZLENGTH , .T = TIME, .x_length = XLENGTH, .y_length = YLENGTH, .z_length = ZLENGTH};

	plotSurface3DReduced(commands, info , sp, arg_time_points, solution1D);
}

double contraction(double solution[SPX][SPY][SPZ][4], int contraction_axis, int contraction_index ){
	/*
	This function will take a 3D matrix and contract dimensions perpendicular to the contraction_axis
	INPUTS:
	solution -- 4D matrix holding the 3D solution at a single slice in time
	contraction_axis -- integer indicating which axis should stay (e.g., if this var == 3, then the x and y dimensions will be integrated out)
	*/

	double integral_holder = 0.0;
	// Integrate over the x and y-axis

	if (contraction_axis == 3){

		for(int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
				integral_holder += solution[i][j][contraction_index][0];

		integral_holder *= HX * HY;	

	}else if (contraction_axis == 2){

		for(int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPZ; ++j)
				integral_holder += solution[i][contraction_index][j][0];

		integral_holder *= HX * HY;	

	}else if (contraction_axis == 1){
		
		for(int i = 0; i < SPY; ++i)
			for (int j = 0; j < SPZ; ++j)
				integral_holder += solution[contraction_index][i][j][0];

		integral_holder *= HX * HY;
	}

	return integral_holder;
}

void plotZDensity3D(double initialCondition[SPX][SPY][SPZ][4]){
    // Plots the number density for the initial condition/g.s. along the z axis

    double zvals[SPZ];
    double zdensity[SPZ];
    double full_solution[SPX][SPY][SPZ][4];

    for(int i = 0; i < SPX; ++i)
    	for(int j = 0; j < SPY; ++j)
    		for(int k = 0; k < SPZ; ++k)
    			full_solution[i][j][k][0] = pow(initialCondition[i][j][k][0], 2);

    for (int i = 0; i < SPZ; ++i)
    {
        zvals[i] = HZ * i;
        zdensity[i] = contraction(full_solution, 3, i);
    }

    char * zcommandsForGnuplot[] =  {"set title \"Z-density\"", "set xlabel \"z-axis\"", "plot 'zdensity.temp' with lines"};
    int num_commands = 3;
    plotFunction(zcommandsForGnuplot, num_commands, zvals, zdensity, SPZ, "zdensity.temp");
}

void plotXDensity3D(double initialCondition[SPX][SPY][SPZ][4]){
    // Plots the number density for the initial condition/g.s. along the z axis

    double xvals[SPX];
    double xdensity[SPX];
    double full_solution[SPX][SPY][SPZ][4];

    for(int i = 0; i < SPX; ++i)
    	for(int j = 0; j < SPY; ++j)
    		for(int k = 0; k < SPZ; ++k)
    			full_solution[i][j][k][0] = pow(initialCondition[i][j][k][0], 2);

    for (int i = 0; i < SPX; ++i)
    {
        xvals[i] = HX * i;
        xdensity[i] = contraction(full_solution, 1, i);
    }

    char * xcommandsForGnuplot[] =  {"set title \"X-density\"", "set xlabel \"x-axis\"", "plot 'xdensity.temp' with lines"};
    int num_commands = 3;
    plotFunction(xcommandsForGnuplot, num_commands, xvals, xdensity, SPX, "xdensity.temp");
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
	plotFunction(commandsForGnuplot, num_commands, xvals, norm, TIME_POINTS, "norm.temp");
}

void saveSolution(int sp, int arg_time_points, double matrix[sp][arg_time_points], char * filename){

	FILE * f = fopen(filename, "w");
	fprintf(f, "%d %d %d\n\n", SPX, SPY, SPZ);
	fprintf(f, "%f %f %f\n", XLENGTH, YLENGTH, ZLENGTH);

	for(int j = 0; j < arg_time_points; ++j)
		for(int i = 0; i < sp; ++i)
			fprintf(f, "%e\n", matrix[i][j]);

	fclose(f);
}

void loadSolution(int sp, int arg_time_points, double matrix[sp][arg_time_points], char * filename){

	FILE * f = fopen(filename, "r");
	int rows, columns, lines;
	double xlen, ylen, zlen, buff; 

	fscanf(f, "%d %d %d", &rows, &columns, &lines);

	if (rows != SPX || columns != SPY || lines != SPZ){			
		printf("File matrix does not match argument matrix in dimension\n");
		fclose(f);
		return;
	}

	fscanf(f, "%lf %lf %lf", &xlen, &ylen, &zlen);

	for(int j = 0; j < arg_time_points; ++j)
		for(int i = 0; i < sp; ++i){
			fscanf(f, "%lf", &buff);
			matrix[i][j] = buff;
		}

	fclose(f);
}

void saveMatrix(double matrix[SPX][SPY][SPZ][4], char * filename){
	/*
	This function saves the numbers in matrix to single column in txt file "filename" -- the first line
	of the output file gives the dimensions of the matrix "row, column"
	*/

	FILE * f = fopen(filename, "w");
	fprintf(f, "%d %d %d\n\n", SPX, SPY, SPZ);
	fprintf(f, "%f %f %f\n", XLENGTH, YLENGTH, ZLENGTH);

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				fprintf(f, "%e\n", matrix[i][j][k][0]);

	fclose(f);
}

void loadMatrix(double matrix[SPX][SPY][SPZ][4], char * filename){
	/*
	Takes matrix data in filename and loads it into matrix
	*/

	FILE * f = fopen(filename, "r");
	int rows, columns, lines;
	double xlen, ylen, zlen, buff; 

	printf("test\n");

	fscanf(f, "%d %d %d", &rows, &columns, &lines);

	if (rows != SPX || columns != SPY || lines != SPZ){			
		printf("File matrix does not match argument matrix in dimension\n");
		fclose(f);
		return;
	}

	fscanf(f, "%lf %lf %lf", &xlen, &ylen, &zlen);

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k){
				fscanf(f, "%lf", &buff);
				matrix[i][j][k][0] = buff;
			}

	fclose(f);
}

// ----------------- Real Time Propagation Functions ---------------

/*double f(double imag_temp[SPX][SPY][SPZ][4], double real_temp[SPX][SPY][SPZ][4], int i, int j, int k, int p, double t){
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][k][p];
	double real_part = real_temp[i][j][k][p];
	double V = .5 *  (pow((i * HX - XLENGTH * 0.5) * WX, 2) + pow((j * HY - YLENGTH * 0.5) * WY, 2) + pow((k * HZ - ZLENGTH * 0.5) * WZ, 2));
	double new_G = (t < T_MOD) ? G * (1. +  EPS * sin(OMEGA * t)) : G;

	return -.5 * (Dxx(imag_temp, i, j, k, p) + Dyy(imag_temp, i, j, k, p) + Dzz(imag_temp, i ,j ,k ,p)) + V * imag_part + new_G * (pow(real_part, 2) + pow(imag_part, 2)) * imag_part;
}*/
/*
double g(double real_temp[SPX][SPY][SPZ][4], double imag_temp[SPX][SPY][SPZ][4], int i, int j, int k, int p, double t){
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][k][p];
	double real_part = real_temp[i][j][k][p];
	double V = .5 *  (pow((i * HX - XLENGTH * 0.5) * WX, 2) + pow((j * HY - YLENGTH * 0.5) * WY, 2) + pow((k * HZ - ZLENGTH * 0.5) * WZ, 2));
	double new_G = (t < T_MOD) ? G * (1. +  EPS * sin(OMEGA * t)) : G;


	return .5 * (Dxx(real_temp, i, j, k, p) + Dyy(real_temp, i, j, k, p) + Dzz(real_temp, i ,j ,k ,p)) - V * real_part - new_G * (pow(real_part, 2) + pow(imag_part, 2)) * real_part;
}*/

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
	int reduction_coeff = RED_COEFF;

	double full_solution[SPX][SPY][SPZ][4]; // This matrix stores the psi squared solution at each time slice to make contraction more convenient
	double solutionXD[SPX][arg_time_points/reduction_coeff];
	double solutionZD[SPZ][arg_time_points/reduction_coeff]; //This matrix contracts the full solution along the radial 'skinny' dimension because it is not essential to observe dynamics

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
		solutionZD[i][0] = contraction(full_solution, 3, i);

	for(int i = 0; i < SPX; ++i)
		solutionXD[i][0] = contraction(full_solution, 1, i);

	double t;

	#define V(i, j, k) (.5*( pow((i * HX - XLENGTH * 0.5) * WX, 2) + pow((j * HY - YLENGTH * 0.5) * WY, 2) + pow((k * HZ - ZLENGTH * 0.5) * WZ, 2)))
	;
	#define NEW_G ( (t < T_MOD) ? G * (1. +  EPS * sin(OMEGA * t)) : G )
	;
	#define f(imag_temp, real_temp, i, j, k, p, t) (-.5*(Dxx(imag_temp, i, j, k, p) + Dyy(imag_temp, i, j, k, p) + Dzz(imag_temp, i ,j ,k ,p)) + \
	V(i, j, k) * imag_temp[i][j][k][p] + NEW_G * (pow(real_temp[i][j][k][p], 2) + pow(imag_temp[i][j][k][p], 2)) * imag_temp[i][j][k][p])
	;
	#define g(real_temp, imag_temp, i, j, k, p, t) (.5 * (Dxx(real_temp, i, j, k, p) + Dyy(real_temp, i, j, k, p) + Dzz(real_temp, i ,j ,k ,p)) - \
	V(i, j, k) * real_temp[i][j][k][p] - NEW_G * (pow(real_temp[i][j][k][p], 2) + pow(imag_temp[i][j][k][p], 2)) * real_temp[i][j][k][p])
	;

	int i, j, k;

	long start, end;
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	// Real Time Propagation
	for (int p = 1; p < arg_time_points; ++p)
	{

		t = Dt * (p - 1); //Forward Euler step
		#pragma omp parallel for collapse(3) num_threads(NUM_THREADS)
		for (i = 0; i < SPX; ++i)
			for (j = 0; j < SPY; ++j)
				for (k = 0; k < SPZ; ++k){
					K[i][j][k][0] = f(imag_solution, real_solution, i, j, k, (p - 1) % 4, t);
					L[i][j][k][0] = g(real_solution, imag_solution, i, j, k, (p - 1) % 4, t);
					real_temp[i][j][k][1] = real_solution[i][j][k][(p - 1) % 4] + .5 * Dt * K[i][j][k][0]; 
					imag_temp[i][j][k][1] = imag_solution[i][j][k][(p - 1) % 4] + .5 * Dt * L[i][j][k][0];
				}

		t = t + .5 * Dt; //Add half a time step 
		#pragma omp parallel for collapse(3) num_threads(NUM_THREADS)
		for (i = 0; i < SPX; ++i)
			for (j = 0; j < SPY; ++j)
				for(k = 0; k < SPZ; ++k){
					K[i][j][k][1] = f(imag_temp, real_temp, i, j, k, 1, t);
					L[i][j][k][1] = g(real_temp, imag_temp, i, j, k, 1, t);
					real_temp[i][j][k][2] = real_solution[i][j][k][(p - 1) % 4] + .5 * Dt * K[i][j][k][1];
					imag_temp[i][j][k][2] = imag_solution[i][j][k][(p - 1) % 4] + .5 * Dt * L[i][j][k][1];
				}


		// t does not change for this step
		#pragma omp parallel for collapse(3) num_threads(NUM_THREADS)
		for (i = 0; i < SPX; ++i)
			for (j = 0; j < SPY; ++j)
				for(k = 0; k < SPZ; ++k){
					K[i][j][k][2] = f(imag_temp, real_temp, i, j, k, 2, t);
					L[i][j][k][2] = g(real_temp, imag_temp, i, j, k, 2, t);
					real_temp[i][j][k][3] = real_solution[i][j][k][(p - 1) % 4] + Dt * K[i][j][k][2];
					imag_temp[i][j][k][3] = imag_solution[i][j][k][(p - 1) % 4] + Dt * L[i][j][k][2];
				}	

		t = Dt * p; //Add full step for Backward Euler step
		#pragma omp parallel for collapse(3) num_threads(NUM_THREADS) private(k_3, l_3)
		for (i = 0; i < SPX; ++i){
			for (j = 0; j < SPY; ++j){
				for(k = 0; k < SPZ; ++k){
					k_3 = f(imag_temp, real_temp, i, j, k, 3, t);
					l_3 = g(real_temp, imag_temp, i, j, k, 3, t);
					real_solution[i][j][k][p % 4] = real_solution[i][j][k][(p - 1) % 4] + 1./6 * Dt * (K[i][j][k][0] + 2 * K[i][j][k][1] + 2 * K[i][j][k][2] + k_3);
					imag_solution[i][j][k][p % 4] = imag_solution[i][j][k][(p - 1) % 4] + 1./6 * Dt * (L[i][j][k][0] + 2 * L[i][j][k][1] + 2 * L[i][j][k][2] + l_3);
				}
			}
		}

		if((p % reduction_coeff) == 0){
			// Add new iteration to full solution
			#pragma omp parallel for collapse(3) num_threads(NUM_THREADS)
			for(int i = 0; i < SPX; ++i)
				for(int j = 0; j < SPY; ++j)
					for(int k = 0; k < SPZ; ++k)
						full_solution[i][j][k][0] = pow(real_solution[i][j][k][p % 4], 2) + pow(imag_solution[i][j][k][p % 4], 2);

			// Contract full solution and save it to solution1D
			for(int i = 0; i < SPZ; ++i)
				solutionZD[i][p/reduction_coeff] = contraction(full_solution, 3, i);

			for(int i = 0; i < SPX; ++i)
				solutionXD[i][p/reduction_coeff] = contraction(full_solution, 1, i);
		}
	}

	gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%ld milliseconds elapsed\n", (end - start));

	saveSolution(SPX, arg_time_points/reduction_coeff, solutionXD, "xsolution.txt");
	saveSolution(SPZ, arg_time_points/reduction_coeff, solutionZD, "zsolution.txt");

	if (plot.plot3D == 1){

		struct plot_settings plotX = {.plot3D = 1, .title = "Real Time Solution Linear X Density", .plot_normalization = 1, .output_file = "solXD.temp"};
		struct plot_settings plotZ = {.plot3D = 1, .title = "Real Time Solution Linear Z Density", .plot_normalization = 1, .output_file = "solZD.temp"};

		// Now lets plot our results
		plotTimeDependence(SPX, arg_time_points/reduction_coeff, solutionXD, plotX);
		plotTimeDependence(SPZ, arg_time_points/reduction_coeff, solutionZD, plotZ);
	}
	if (plot.plot_normalization == 1){

		plotNormalization(arg_time_points/reduction_coeff, solutionZD, T);
	}
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = .01; // This is the time step just used by the imaginary time propagation method

double PsiNorm(double Array[SPX][SPY][SPZ][4]){
	// Integrates down the first x-y matrix assuming spatial square has area HX * HY

	double integral = 0.0;

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				integral += pow(Array[i][j][k][0], 2);

	integral *= HX * HY * HZ;

	return integral;
}

double initialGuess(double x, double y, double z){

	return x * (XLENGTH - x) * y * (YLENGTH - y) * z * (ZLENGTH - z);
}

void normalize(double Array[SPX][SPY][SPZ][4]){
	// Takes an array with a function in the first column and normalizes that function back to 1

	double norm = 1 / sqrt( PsiNorm(Array) );

	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				Array[i][j][k][0] = Array[i][j][k][0] * norm;
}

void advance(double real_solution[SPX][SPY][SPZ][4], int p ){
	// Advance the solution matrix one index forward in time where p is the solution index of the currently solved time

	int next = (p + 1) % 4;
	double real_part;
	double V;
	double lhs;

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k){
				real_part = real_solution[i][j][k][p];
				V = pow((i * HX - XLENGTH / 2) * WX, 2) + pow((j * HY - YLENGTH / 2) * WY, 2) + pow((k * HZ - ZLENGTH / 2) * WZ, 2);
				lhs = (Dxx(real_solution, i, j, k, p) + Dyy(real_solution, i, j, k, p) + Dzz(real_solution, i, j, k, p)) / 2  - V * real_part / 2 - G * pow(real_part , 3);
				real_solution[i][j][k][next] = real_solution[i][j][k][p] + Delta_t * lhs;
			}
}

void findGroundState(double real_solution[SPX][SPY][SPZ][4], int iterations) {

	// Assign initial guess
	for (int i = 0; i < SPX; ++i)
		for (int j = 0; j < SPY; ++j)
			for(int k = 0; k < SPZ; ++k)
				real_solution[i][j][k][0] = initialGuess(i * HX, j * HY, k * HZ);

	normalize(real_solution);

	loadMatrix(real_solution, "groundState3D.txt");

	#define fgs(real_temp, i, j, k, p) ((Dxx(real_temp, i, j, k, p) + Dyy(real_temp, i, j, k, p) + Dzz(real_temp, i, j, k, p)) / 2  - V(i, j, k) * real_temp[i][j][k][p] - G * pow(real_temp[i][j][k][p] , 3))
	;
	
	long start, end;
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

	for (int p = 1; p < iterations; ++p){

		#pragma omp parallel for collapse(3) num_threads(NUM_THREADS)
		for(int i = 0; i < SPX; ++i)	
			for(int j = 0; j < SPY; ++j)
				for(int k = 0; k < SPZ; ++k){
					real_solution[i][j][k][p % 4] = real_solution[i][j][k][(p - 1) % 4] + Delta_t * fgs(real_solution, i, j, k, (p - 1) % 4);
				}

		// Now normalize the solution back to 1
		if(p % 20 == 0)	
			normalize(real_solution);
	}

	gettimeofday(&timecheck, NULL);
    end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

    printf("%ld milliseconds elapsed\n", (end - start));

	normalize(real_solution);

	saveMatrix(real_solution, "groundState3D.txt");
}

int main(){

	// Declare initial variables
	double real_solution[SPX][SPY][SPZ][4]; // This will be the solution matrix -- it will only hold the state of the system at one time slice
	double imag_solution[SPX][SPY][SPZ][4]; // Same as real solution but for the imaginary component of Psi
	double initialCondition[SPX][SPY][SPZ][4];

	struct plot_settings plot_solution = {.plot3D = 0, .title = "Real Time Solution", .plot_normalization = 1};

	// Load the groundstate
	loadMatrix(initialCondition, "groundState3D.txt");

	// find the ground state
	// findGroundState(initialCondition, 1000);

	// plotZDensity3D(initialCondition);
	// plotXDensity3D(initialCondition);
	// Run RTP
	realTimeProp(initialCondition, TIME, TIME_POINTS, real_solution, imag_solution, plot_solution);


	return 0;
}