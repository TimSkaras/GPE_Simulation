/*
This program numerically solves the non-linear schrodinger equation using a 4th order Runge Kutta Scheme
and then plots the solution using a surface plot command on gnuplot

Important constants are defined as preprocessor constants (like the desired length of time for the solution)
The preprocessor constant G is a parameter of the system proportional to the scattering length

This program numerically solves homogeneous initial conditions with perturbations. The wave number of the 
perturbation is selected by the parameter K_P. The program will add noise with a wave number having
an integer number of wave lengths in the well.

Note: We expect exponential growth whenever 0 < K_P < sqrt(4 * |G|) * A (where G < 0) and stability otherwise.
I.e., the perturbations will always be stable if G is positive

TO DO:
1) Measure energy to make sure energy is in fact decreasing for each iteration of findgroundstate

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "..\Plot.h"

#define PI M_PI
#define TIME 3.0
#define LENGTH 15.0
#define TIME_POINTS 12000 	//Time step
#define SPACE_POINTS 450 	//spatial step
#define K_P 0.25 // Wave number of the perturbation


// This is just a macro for finding the fourth ourder 
#define Dxx(array, index, pee) ((-1. * array[mod((index) + 2, space_points)][(pee)] + 16.* array[mod((index) + 1, space_points)][(pee)] - \
		30. * array[mod((index) , space_points)][(pee)] + 16. * array[mod((index) - 1, space_points)][(pee)] +  -1 * array[mod((index) - 2, space_points)][(pee)])/(12. * pow(H, 2)))
;
const double H = LENGTH / (SPACE_POINTS);
const double G = 250.; // Scattering Length parameter
const double W = 1.5; //Trapping frequency
const double OMEGA = 1.; //Perturbation frequency
const double EPS = 0.0; //Perturbation amplitude

// How many time steps do we want
const int time_points = TIME_POINTS;

// How many spatial points are there
const int space_points = SPACE_POINTS;

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

// ---------------Auxiliary Functions-----------------

int mod(int a, int b){
	// This is a standard modulus function -- always outputs a positive number
	// E.g., mod(-1, 3) = 2

    int r = a % b;
    return r < 0 ? r + b : r;
}

void plotSurface(char * commandsForGnuplot[], int num_commands, int space_points, int arg_time_points, double solution[][arg_time_points], double T, char * output_file ){
	/*
	This function gives a 3-dimensional plot of data to 1+1D problem. Plots position versus time. Solution is specified using a 2-dimensional array where each row is a position 
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

    double Dt = T / arg_time_points;

	double time_samples = 70.;
	double space_samples = 70.;

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

    FILE * temp = fopen(output_file, "w");
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

    for (int i=0; i < num_commands; i++)
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
}

void plotNormalization(int space_points, int arg_time_points, double Array[SPACE_POINTS + 1][arg_time_points], double T){
	/*
	
	This function plots the normalization of the numerical solution at each point in time using a midpoint integration algorithm

	The last row in Array is just a copy of the first row to make the plot look correct (there is a discrepancy due to the periodic boundary conditions)
	As a result, the last row should not be included when calculating the norm because it has no physical significance.

	*/
	double Dt = T / arg_time_points;

	double xvals[arg_time_points];
	double norm[arg_time_points];
	double integral, baseline;

	baseline = 0.;

	for (int i = 0; i < space_points; ++i)
		baseline = baseline +  H * (Array[i][0]);

	for (int p = 0; p < arg_time_points; ++p)
	{
		integral = 0.;

		for (int i = 0; i < space_points; ++i)
			integral = integral + (Array[i][p]);
		
		integral = H * integral;

		norm[p] = integral - baseline;
	}

	// Create array of x values to plot then norm against
	for (int p = 0; p < arg_time_points; ++p)
	{
		xvals[p] = p * Dt;
	}

	char * commandsForGnuplot[] = {"set title \"Normalization Error G=100-GS vs. Time\"", "set xlabel \"Time\"", "set ylabel \"Normalization Error\"" , "plot 'norm.temp' with lines"};
	int num_commands = 4;
	plotFunction(commandsForGnuplot, num_commands, xvals, norm, arg_time_points, "norm.temp");
}

// -----------------Real Time Propagation Functions-----------------

double f(double imag_temp[SPACE_POINTS][4], double real_temp[SPACE_POINTS][4], int i, int p, double t){
	
	// f gives the derivative for psi_real but for the K matrices because they need a function with a different argument
	
	double psi_i = imag_temp[i][p];
	double psi_r = real_temp[i][p];
	double new_G = G * (1. + EPS * sin(OMEGA * t));
	return -1 * Dxx(imag_temp, i, p) / 2. + 1/2. * pow((i * H - LENGTH / 2.0) * W, 2) * psi_i + new_G * (pow(psi_r, 2) + pow(psi_i, 2)) * psi_i;
}

double g(double real_temp[SPACE_POINTS][4], double imag_temp[SPACE_POINTS][4], int i, int p, double t){
	
	// g gives the derivative for imag_temp but for the L matrices because they need a function with a different argument
	
	double psi_i = imag_temp[i][p];
	double psi_r = real_temp[i][p];
	double new_G = G * (1. + EPS * sin(OMEGA * t));
	return Dxx(real_temp, i, p ) / 2.  - 1/2. * pow((i * H - LENGTH / 2.0) * W, 2) * psi_r - new_G * (pow(psi_r, 2) + pow(psi_i, 2)) * psi_r;
}

void addNoise(double Array[][4], double waveNumber){
	// Adds noise to the first column of Array with the desired wave number

	// First we need to find noise with integer number of wavelengths
	double num = waveNumber * LENGTH / PI;
	int n = (int) floor(num);

	if (mod(n, 2) != 0)
		n = (int) ceil(num);

	double kp_eff = n * PI / LENGTH;

	// Noise should be small relative to homogeneous solution
	double noise_volume = .00;

	// Loop through array and add noise
	for (int i = 0; i < space_points; ++i)
		Array[i][0] = Array[i][0] + noise_volume * sin(kp_eff * i * H);
}

void realTimeProp(double initialCondition[SPACE_POINTS], double T, int arg_time_points, double real_solution[][arg_time_points], double imag_solution[][arg_time_points], struct plot_settings plot){
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

	double real_temp[space_points][4];
	double imag_temp[space_points][4];
	double K[space_points][3];
	double L[space_points][3];
	double k_3, l_3;

	double full_solution[space_points + 1][arg_time_points]; // This will have the psi square solution so we can plot the results
	
	// Initialize arrays because apparently C needs you to do this
	for (int i = 0; i < space_points; ++i){

		for (int j = 0; j < arg_time_points; ++j)
		{
			real_solution[i][j] = 0;
			imag_solution[i][j] = 0;
		}

		for (int j = 0; j < 3; ++j)
		{
			real_temp[i][j] = 0;
			imag_temp[i][j] = 0;
			K[i][j] = 0;
			L[i][j] = 0;
		}
	}

	// Assign initial conditions
	for (int i = 0; i < space_points; ++i)
		real_solution[i][0] = initialCondition[i];

	double t;
	// Real time evolution
	for (int p = 1; p < arg_time_points; ++p){

		// Load solution into first column of temp matrices
		for(int i = 0; i < space_points; ++i){
			real_temp[i][0] = real_solution[i][p - 1];
			imag_temp[i][0] = imag_solution[i][p - 1];
		}

		t = Dt * (p - 1); //Forward Euler step
		for (int i = 0; i < space_points; ++i)
		{
			K[i][0] = f(imag_temp, real_temp, i, 0, t);
			L[i][0] = g(real_temp, imag_temp, i, 0, t);
			real_temp[i][1] = real_solution[i][p - 1] + .5 * Dt * K[i][0];
			imag_temp[i][1] = imag_solution[i][p - 1] + .5 * Dt * L[i][0];
		}

		t = t + .5 * Dt; //Add half a time step 
		for (int i = 0; i < space_points; ++i)
		{
			K[i][1] = f(imag_temp, real_temp, i, 1, t);
			L[i][1] = g(real_temp, imag_temp, i, 1, t);
			real_temp[i][2] = real_solution[i][p - 1] + .5 * Dt * K[i][1];
			imag_temp[i][2] = imag_solution[i][p - 1] + .5 * Dt * L[i][1];
		}

		// t does not change for this step
		for (int i = 0; i < space_points; ++i)
		{
			K[i][2] = f(imag_temp, real_temp, i, 2, t);
			L[i][2] = g(real_temp, imag_temp, i, 2, t);
			real_temp[i][3] = real_solution[i][p - 1] + Dt * K[i][2];
			imag_temp[i][3] = imag_solution[i][p - 1] + Dt * L[i][2];
		}

		t = Dt * p; //Add full step for Backward Euler step
		for (int i = 0; i < space_points; ++i)
		{
			k_3 = f(imag_temp, real_temp, i, 3, t);
			l_3 = g(real_temp, imag_temp, i, 3, t);
			real_solution[i][p] = real_solution[i][p - 1] + 1./6 * Dt * (K[i][0] + 2 * K[i][1] + 2 * K[i][2] + k_3);
			imag_solution[i][p] = imag_solution[i][p - 1] + 1./6 * Dt * (L[i][0] + 2 * L[i][1] + 2 * L[i][2] + l_3);
		}
		
	}

	if (plot.plot3D == 1){
		// Find the psi squared solution so we can plot the results
		for (int i = 0; i < space_points; ++i)
			for (int j = 0; j < arg_time_points; ++j)
				full_solution[i][j] = pow(real_solution[i][j], 2) + pow(imag_solution[i][j], 2);


		// Find the difference from g.s. for diagnostic purposes
		for(int i = 0; i < space_points; ++i)
			for(int j = 0; j < arg_time_points; ++j)
				full_solution[i][j] = full_solution[i][j] - initialCondition[i] * initialCondition[i];

		
		printf("%e\n", full_solution[20][10]);
		// Add the last row for cosmetic purposes
		for (int i = 0; i < arg_time_points; ++i)
			full_solution[SPACE_POINTS][i] = full_solution[0][i];

		// Now lets plot our results

		// Set up commands for gnuplot
		char plotCommand[100];
		char title[100];
		// sprintf(plotCommand, "plot [0:%d] [-2:2] 'sol.temp' ", TIME);
		strcpy(plotCommand,"splot 'sol.temp' with lines");
		sprintf(title, "set title \"");
		strcat(title, plot.title);
		strcat(title, "\"");
		char * commandsForGnuplot[] = { title, "set xlabel \"Position\"", "show xlabel", "set ylabel \"Time\"", "show ylabel" , plotCommand};
		int num_commands = 6;


		plotSurface(commandsForGnuplot, num_commands, space_points + 1, arg_time_points, full_solution, T,  "sol.temp");
	}
	if (plot.plot_normalization == 1){

		plotNormalization(space_points, arg_time_points, full_solution, T);
	}
}

// -------------Imaginary Time Propagation Functions------------

const double Delta_t = (LENGTH / SPACE_POINTS)/600.;

double PsiNorm(double Array[][4]){
	// Integrates down the first column assuming spatial step is H

	double integral = 0.0;

	for (int i = 1; i < SPACE_POINTS; ++i)
		integral = integral + .5 * (pow(Array[i - 1][0], 2) + pow(Array[i][0], 2)) * H; 

	return integral;
}

double initialGuess(double x){

	return x * (LENGTH - x);
}

void normalize(double Array[][4]){
	// Takes an array with a function in the first column and normalizes that function back to 1

	double norm = PsiNorm(Array);
	norm = sqrt(norm);

	for (int i = 0; i < SPACE_POINTS; ++i)
		Array[i][0] = Array[i][0] / norm;
}

double fgs(double Psi_real[][4], int i, int p){
	// f is the function that gives the derivative for psi_real, which is why it is a function of imag_solution
	return 0.5 * Dxx(Psi_real, i, p) - 0.5 * pow(W * (i * H - LENGTH/2.), 2) * Psi_real[i][p] - G * pow(Psi_real[i][p], 2) * Psi_real[i][p];
}

void findGroundState(double real_solution[][4], int space_points, int iterations) {

	// Declare initial variables
	double real_temp[space_points][4];
	double ground_state[space_points];

	// Assign initial gaussian profile
	for (int i = 0; i < space_points; ++i)
		real_solution[i][0] = initialGuess(i * H);

	normalize(real_solution);

	for (int p = 1; p < iterations; ++p)
	{

		for(int i = 0; i < space_points; ++i)
			real_solution[i][1] = real_solution[i][0] + Delta_t * fgs(real_solution, i, 0);

		// Move solution back an index so we can repeat the process
		for (int i = 0; i < space_points; ++i)
		{
			real_solution[i][0] = real_solution[i][1];
			real_solution[i][1] = 0.0;
		}

		// Now normalize the solution back to 1
		if (mod(p, 5) == 0)
			normalize(real_solution);
	}
	normalize(real_solution);
}

void stepTwo(double real_solution[][4], int space_points){
	/*
	*	This function takes an output from findGroundstate that is sufficiently close to the actual ground state and runs
	*	this output through real time propagation and averages the solution over real time to remove any further error from g.s.
	*
	*	The method works by propagating the g.s. through a fixed period of time, averaging over each value the solution obtains,
	*	and then repeating this process by iterting that g.s. through the same fixed period of time
	*
	*/

	int iterative_period = 6000;
	int trials = 30;
	double T = 6.0;

	double initialCondition[space_points];

	double real_sample[space_points][iterative_period];
	double imag_sample[space_points][iterative_period];

	for(int i = 0; i < space_points; ++i)
		initialCondition[i] = real_solution[i][0];

	struct plot_settings trialPlot = {.plot3D = 0, .title = "Step Two Iterative Period", .plot_normalization = 0};
	struct plot_settings stepTwoPlot = {.plot3D = 1, .title = "Twenty Averaging Trials", .plot_normalization = 0};

	for (int n = 0; n < trials; ++n)
	{
		realTimeProp(initialCondition, T, iterative_period, real_sample, imag_sample, trialPlot);

		for (int i = 0; i < space_points; ++i)	
			for (int j = 1; j < iterative_period; ++j)
				initialCondition[i] = initialCondition[i] + sqrt(pow(real_sample[i][j], 2) + pow(imag_sample[i][j], 2));

		for(int i = 0; i < space_points; ++i)
			initialCondition[i] = initialCondition[i]/iterative_period;
			
	}


	for(int i = 0; i < space_points; ++i)
		real_solution[i][0] = initialCondition[i];

	normalize(real_solution);

	realTimeProp(initialCondition, T, iterative_period, real_sample, imag_sample, stepTwoPlot );
}

int main() {

	double Dt = TIME / TIME_POINTS;

	// Declare initial variables
	double real_solution[SPACE_POINTS][TIME_POINTS]; // This will be solution matrix where each column is a discrete point in time and each row a discrete point in space
	double imag_solution[SPACE_POINTS][TIME_POINTS]; // Same as real solution but for the imaginary component of Psi
	double ground_state[SPACE_POINTS][4];
	double initial_condition[SPACE_POINTS];

	findGroundState(ground_state, space_points, 40000);

	stepTwo(ground_state, space_points);

	for(int i = 0; i < space_points; ++i)
		initial_condition[i] = ground_state[i][0];

	struct plot_settings plot_solution = {.plot3D = 1, .title = "Difference Between Real Time Solution and ITP G.S.", .plot_normalization = 1};

	// realTimeProp(initial_condition, TIME, time_points, real_solution, imag_solution, plot_solution);

	double psi_r, psi_i;
	for(int i = 0; i < space_points; ++i){
		psi_r = real_solution[i][time_points - 1];
		psi_i = imag_solution[i][time_points - 1];
		ground_state[i][0] = sqrt(psi_r * psi_r + psi_i * psi_i);
	}
		
	// stepTwo(ground_state, space_points);

	printf("%d %d \nSpace Step: %.5e \nTime Step: %.5e \nRatio: %f", space_points, time_points, H, Dt,  H * H / Dt);
	
	return 0;
}

