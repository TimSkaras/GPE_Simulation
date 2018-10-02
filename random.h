#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

void seed(){

	// Choose pseudo-random seed
    int time_number;
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    time_number = (int)timecheck.tv_sec * 1000 + (int)timecheck.tv_usec / 1000;

    int seed = mod(time_number, 10000);
    srand(seed);
}

double generateGaussianNoise(double mu, double sigma){
	/*
	*	This function samples a random variable with a gaussian distribution
	*
	*	INPUTS:
	*	mu -- double describing average of desired distribution
	*	sigma -- double describing standard deviation of desired distribution
	*
	*
	*	NOTE: a seed must be picked prior to using this function, otherwise srand(1) is assumed, i.e., seed = 1
	*/


    double two_pi = 2.0 * M_PI;
    double u1, u2;

    for (int i = 0; i < 5; ++i){
       u1 = rand() * (1.0 / RAND_MAX);
       u2 = rand() * (1.0 / RAND_MAX);
    }

    double z0;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    return z0 * sigma + mu;
}
