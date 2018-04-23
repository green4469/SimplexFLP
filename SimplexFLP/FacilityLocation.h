#ifndef __FL__
#define __FL__
#include <iostream>
#include <cstdlib>
#include <ilcplex/ilocplex.h>
#include <random>
#include <string>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <numeric>

using namespace std;

class FacilityLocation {
public:
	/* The number of facilities and clients */
	int n_facilities = 0;
	int n_clients = 0;

	/* Rounded solution's objective cost */
	double rounded_cost = 0.0;

	/* Original problem's optimal objective cost */
	double optimal_cost = 0.0;

	/* Input of LP-solver */
	double *opening_cost = NULL;			// n_facilities
	double **connection_cost = NULL;		// n_facilities x n_clients

	/* Output of LP-solver */
	double **clients_coordinate = NULL;		// n_clients x n_facilities

	/* exponential clocks of facilities */
	double *exponential_clocks = NULL;		// n_facilities

	/* Output of rounding algorithm */
	bool *opening_table = NULL;				// n_facilities
	bool **connection_table = NULL;			// n_facilities x n_clients

	/* Constructor */
	FacilityLocation(int argc, char *argv[], int dim=0);

	/* Deconstructor */
	~FacilityLocation();

	/* LP-solver */
	double LP_solve();
	double get_optimal();

	/* Rounding algorithm */
	void rounding_alg();

	/* Use competing exponential clocks of facilities */
	double rounding_alg_exp();

	/* Use distortion of the simplex */
	double rounding_alg_dist();

	/* Triangular inequality */
	void triangular_inequality();
};
#endif
