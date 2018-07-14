/*

This program has the capability to find the ground state for arbitrary Hamiltonian using the Imaginary 
Time Propagation Method (ITP), advance that ground state in real time using RK4, and plot the resulting time evolution.
The wave number of the perturbation is selected by the parameter K_P. The program will add noise 
with a wave number having an integer number of wave lengths in the well.

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant G is a parameter of the system proportional to the scattering length

Important Results:
1) spx = 600, SPY = 40, EPS = 0.2, noise = on -- the wave function becomes very noisy after 100 s or so
2) spx = 600, SPY = 40, EPS = 0, noise = on -- the wave function stays coherent but noise waves travel to outside over time
3) spx = 900, SPY = 45 (Ylength = 15), EPS = .01, noise on -- small oscillations visible in y density, noise dissipates, no instability observed
4) spx = 900, SPY = 45 (Ylength = 15), EPS = .08, noise on (changed K: 0.7 -> 0.5) -- noise dissipates, small oscillation in y-density
5) spx = 900, SPY = 45 (Ylength = 15), EPS = .08, noise on (k = 0.5) -- noise dissipates, medium oscillation in y-denisty, new wavelength starts to emerge at end of simulation with k = 0.6444
6) spx = 900, SPY = 45 (Ylength = 15), EPS = .08, noise on (changed K: 0.5 -> 0.65) -- wave noise grows a good deal, nothing of note in y-density, not too unstable
7) spx = 1050 (300 --> 350), SPY = 45 (Ylength = 15), EPS = .08, noise on (changed K: 0.5 -> 0.65) -- 

TODO:
1) Create stepTwo function for 2D -- postpone
2) Find way to save ground state -- done
3) Include fft and associated analysis

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 0.00000200000
#define XLENGTH 350.0 //4.0
#define YLENGTH 15.0 //0.1
#define TIME_POINTS 10 //number of time points
#define SPX 1050 //600
#define SPY 45 //20
#define NOISE_VOLUME 0.06 // 0.06

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
const double EPS = 0.12;
const double WAVENUMBER_INPUT = 0.65;
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

void plotTimeDependence( int sp, int arg_time_points, double solution1D[sp][arg_time_points], struct plot_settings plot){

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
	
	struct plot_information info = { .num_commands = num_commands, \
		.output_file = plot.output_file, .length = (sp == SPX) ? XLENGTH:YLENGTH , .T = TIME, .x_length = XLENGTH, .y_length = YLENGTH};

	plotSurface2DReduced(commands, info , sp, arg_time_points, solution1D);
}

void plotNormalization(int arg_time_points, double solution1D[SPX][arg_time_points], double T){
	/*
	* This function takes the squared solution matrix, finds the normalization at each point in time and plots the error (normalization should equal unity)
	*/

	double norm[arg_time_points];
	double xvals[arg_time_points];
	double Dt = T/arg_time_points;

	for(int p = 0; p < arg_time_points; ++p){
		norm[p] = 0.0;
		for(int i = 0; i < SPX; ++i)
			norm[p] = norm[p] + solution1D[i][p];
		norm[p] = norm[p] * HX;
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
	plotFunction(commandsForGnuplot, num_commands, xvals, norm, arg_time_points, "norm.temp");
}

double contraction(double solution[SPX][SPY][4], int contraction_axis, int i ){
	/*
	This function will take a 3D matrix and contract across the axis perpendicular to contraction_axis
	*/

	double integral = 0.0;
	// Integrate over the y-axis

	if( contraction_axis == 1){
		for (int j = 0; j < SPY; ++j)
			integral = integral + solution[i][j][0];

		integral = integral * HY;
	}else if(contraction_axis == 2){

		for (int j = 0; j < SPX; ++j)
			integral = integral + solution[j][i][0];

		integral = integral * HX;
	}	

	return integral;
}

void saveMatrix(double matrix[SPX][SPY][4], char * filename){
	/*
	This function saves the numbers in matrix to single column in txt file "filename" -- the first line
	of the output file gives the dimensions of the matrix "row, column"
	*/

	FILE * f = fopen(filename, "w");
	fprintf(f, "%d %d\n\n", SPX, SPY);
	fprintf(f, "%f %f\n", XLENGTH, YLENGTH);

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j)
			fprintf(f, "%e\n", matrix[i][j][0]);

	fclose(f);
}

void loadMatrix(double matrix[SPX][SPY][4], char * filename){
	/*
	Takes matrix data in filename and loads it into matrix
	*/

	FILE * f = fopen(filename, "r");
	int rows, columns;
	double xlen, ylen, buff; 

	printf("test\n");

	fscanf(f, "%d %d", &rows, &columns);

	if (rows != SPX || columns != SPY){			
		printf("File matrix does not match argument matrix in dimension\n");
		fclose(f);
		return;
	}

	fscanf(f, "%lf %lf", &xlen, &ylen);

	for(int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j){
			fscanf(f, "%lf", &buff);
			matrix[i][j][0] = buff;
		}

	fclose(f);
}

// ----------------- Real Time Propagation Functions ---------------

double f(double imag_temp[SPX][SPY][4], double real_temp[SPX][SPY][4], int i, int j, int p, double t){
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));
	double new_G = (t < T_MOD) ? G * (1. +  EPS * sin(OMEGA * t)) : G;

	return -.5 * (Dxx(imag_temp, i, j, p) + Dyy(imag_temp, i, j, p)) + V * imag_part + new_G * (pow(real_part, 2) + pow(imag_part, 2)) * imag_part;
}

double g(double real_temp[SPX][SPY][4], double imag_temp[SPX][SPY][4], int i, int j, int p, double t){
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double imag_part = imag_temp[i][j][p];
	double real_part = real_temp[i][j][p];
	double V = .5 *  (pow((i * HX - XLENGTH/2.0) * WX, 2) + pow((j * HY - YLENGTH / 2.0) * WY, 2));
	double new_G = (t < T_MOD) ? G * (1. +  EPS * sin(OMEGA * t)) : G;

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
			initialCondition[i][j][0] = initialCondition[i][j][0] * (1 + noise_volume * sin(kp_eff * i * HX) * exp(- .001 * pow(i * HX - XLENGTH/2., 2) ));
}

void realTimeProp(double initialCondition[SPX][SPY][4], double T, int arg_time_points, double real_solution[][SPY][4], double imag_solution[][SPY][4], struct plot_settings plot){
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

	double real_temp[SPX][SPY][4]; // This matrix and the one below will be used to store the state of the system at intermediary time steps -- time steps 
	double imag_temp[SPX][SPY][4]; // that are need to do RK4 but are not essential otherwise and can therefore be discarded after each iteration
	double K[SPX][SPY][3];
	double L[SPX][SPY][3];
	double k_3, l_3;
	int reduction_coeff = RED_COEFF;

	// printf("%d\n", arg_time_points/reduction_coeff);
	
	double full_solution[SPX][SPY][4]; // This matrix stores the psi squared solution at each time slice to make contraction more convenient
	double solutionXD[SPX][arg_time_points/reduction_coeff]; //This matrix contracts the full solution along the radial 'skinny' dimension
	double solutionYD[SPY][arg_time_points/reduction_coeff]; //This matrix contracts the full solution along the radial 'axial' dimension

	// Add the noise to the ground state
	addNoise(initialCondition, WAVENUMBER_INPUT);

	// Assign initial conditions
	for (int i = 0; i < SPX; ++i)
		for(int j = 0; j < SPY; ++j){
			real_solution[i][j][0] = initialCondition[i][j][0];
			full_solution[i][j][0] = pow(initialCondition[i][j][0], 2);
		}
	
	for(int i = 0; i < SPX; ++i)
		solutionXD[i][0] = contraction(full_solution, 1, i);
	
	for(int i = 0; i < SPY; ++i)
		solutionYD[i][0] = contraction(full_solution, 2, i);

	
	double t;
	// Real Time Propagation
	for (int p = 1; p < arg_time_points; ++p)
	{
		
		// Load solution into first column of temp matrices
		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j){
				real_temp[i][j][0] = real_solution[i][j][0];
				imag_temp[i][j][0] = imag_solution[i][j][0];
		}

		
		t = Dt * (p - 1); //Forward Euler step
		for (int i = 0; i < SPX; ++i){
			for (int j = 0; j < SPY; ++j){

				K[i][j][0] = f(imag_temp, real_temp, i, j, 0, t);
				L[i][j][0] = g(real_temp, imag_temp, i, j, 0, t);
				real_temp[i][j][1] = real_solution[i][j][0] + .5 * Dt * K[i][j][0]; 
				imag_temp[i][j][1] = imag_solution[i][j][0] + .5 * Dt * L[i][j][0];
			}
		}

		t = t + .5 * Dt; //Add half a time step 
		for (int i = 0; i < SPX; ++i){
			for (int j = 0; j < SPY; ++j){

				K[i][j][1] = f(imag_temp, real_temp, i, j, 1, t);
				L[i][j][1] = g(real_temp, imag_temp, i, j, 1, t);
				real_temp[i][j][2] = real_solution[i][j][0] + .5 * Dt * K[i][j][1];
				imag_temp[i][j][2] = imag_solution[i][j][0] + .5 * Dt * L[i][j][1];
			}
		}

		// t does not change for this step
		for (int i = 0; i < SPX; ++i){
			for (int j = 0; j < SPY; ++j){

				K[i][j][2] = f(imag_temp, real_temp, i, j, 2, t);
				L[i][j][2] = g(real_temp, imag_temp, i, j, 2, t);
				real_temp[i][j][3] = real_solution[i][j][0] + Dt * K[i][j][2];
				imag_temp[i][j][3] = imag_solution[i][j][0] + Dt * L[i][j][2];
			}	
		}

		t = Dt * p; //Add full step for Backward Euler step
		for (int i = 0; i < SPX; ++i){
			for (int j = 0; j < SPY; ++j){

				k_3 = f(imag_temp, real_temp, i, j, 3, t);
				l_3 = g(real_temp, imag_temp, i, j, 3, t);
				real_solution[i][j][1] = real_solution[i][j][0] + 1./6 * Dt * (K[i][j][0] + 2 * K[i][j][1] + 2 * K[i][j][2] + k_3);
				imag_solution[i][j][1] = imag_solution[i][j][0] + 1./6 * Dt * (L[i][j][0] + 2 * L[i][j][1] + 2 * L[i][j][2] + l_3);
			}
		}


		// Contract full solution and save it to solution1D
		if (mod(p, reduction_coeff) == 0){

			// Add new iteration to full solution
			for(int i = 0; i < SPX; ++i)
				for(int j = 0; j < SPY; ++j)
					full_solution[i][j][0] = pow(real_solution[i][j][1], 2) + pow(imag_solution[i][j][1], 2);
			for(int i = 0; i < SPX; ++i)
				solutionXD[i][p/reduction_coeff] = contraction(full_solution, 1, i);

			for(int i = 0; i < SPY; ++i)
				solutionYD[i][p/reduction_coeff] = contraction(full_solution, 2, i);
		}
	
		// Move new iteration in real/imag solution back an index to restart process
		for(int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
			{
				real_solution[i][j][0] = real_solution[i][j][1];
				imag_solution[i][j][0] = imag_solution[i][j][1];
			}
		
	}

	if (plot.plot3D == 1){

		struct plot_settings plotX = {.plot3D = 1, .title = "Real Time Solution Linear X Density", .plot_normalization = 1, .output_file = "solXD.temp"};
		struct plot_settings plotY = {.plot3D = 1, .title = "Real Time Solution Linear Y Density", .plot_normalization = 1, .output_file = "solYD.temp"};

		// Now lets plot our results
		plotTimeDependence( SPX, arg_time_points/reduction_coeff, solutionXD, plotX);
		plotTimeDependence(SPY, arg_time_points/reduction_coeff, solutionYD, plotY);

	}
	if (plot.plot_normalization == 1){

		plotNormalization(arg_time_points/reduction_coeff, solutionXD, T);
	}
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = .0001; // This is the time step just used by the imaginary time propagation method

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

	for (int p = 1; p < iterations; ++p){

		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j)
				real_solution[i][j][1] = real_solution[i][j][0] + Delta_t * fgs(real_solution, i, j, 0);

		for(int i = 0; i < SPX; ++i)
			for(int j = 0; j < SPY; ++j)
				real_solution[i][j][2] = real_solution[i][j][1] + Delta_t * fgs(real_solution, i, j, 1);

		// Move solution back an index so we can repeat the process
		for (int i = 0; i < SPX; ++i)
			for (int j = 0; j < SPY; ++j)
			{
				real_solution[i][j][0] = real_solution[i][j][2];
				real_solution[i][j][1] = 0.0;
				real_solution[i][j][2] = 0.0;
			}


		// Now normalize the solution back to 1
		normalize(real_solution);
	}

	// Write the g.s. to a separate file
	saveMatrix(real_solution, "groundState2D.txt");
}

int main(){

	// Declare initial variables
	static double real_solution[SPX][SPY][4]; // This will be the solution matrix -- it will only hold the state of the system at one time slice
	static double imag_solution[SPX][SPY][4]; // Same as real solution but for the imaginary component of Psi
	double initialCondition[SPX][SPY][4];

	struct plot_settings plot_solution = {.plot3D = 1, .title = "Real Time Solution", .plot_normalization = 1, .output_file = "sol.temp"};

	loadMatrix(initialCondition, "groundState2D.txt");

	// find the ground state
	// findGroundState(initialCondition, 30000);

	// Run RTP
	realTimeProp(initialCondition, TIME, TIME_POINTS, real_solution, imag_solution, plot_solution);

	return 0;
}