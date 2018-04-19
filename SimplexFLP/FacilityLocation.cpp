#include "FacilityLocation.h"

/* Considering floating point error, */
int CompareDoubleUlps(double x, double y, int ulpsTolerance = 4)
{
	double diff = x - y;

	__int64 nx = *((__int64*)&x);
	__int64 ny = *((__int64*)&y);

	if ((nx & 0x8000000000000000) != (ny & 0x8000000000000000))
	{
		if (x == y)
			return 0;

		return (diff > 0) ? 1 : -1;
	}

	__int64 ulpsDiff = nx - ny;
	if ((ulpsDiff >= 0 ? ulpsDiff : -ulpsDiff) <= ulpsTolerance)
		return 0;

	return (diff > 0) ? 1 : -1;
}

/* Constructor definition */
FacilityLocation::FacilityLocation(int argc, char *argv[])
{
	ifstream in(argv[1]);
	in >> n_facilities;
	in >> n_clients;

	/* Memory allocation of opening_cost(f), connection_cost(d) */
	opening_cost = new double[n_facilities];
	for (int i = 0; i < n_facilities; i++) {
		in >> opening_cost[i];
	}

	connection_cost = new double*[n_facilities];
	for (int i = 0; i < n_facilities; i++) {
		connection_cost[i] = new double[n_clients];
		for (int j = 0; j < n_clients; j++) {
			connection_cost[i][j] = DBL_MAX;
		}
	}

	int f, c;
	while (in) {
		in >> f >> c;
		in >> connection_cost[f][c];
	}

	in.close();

	/* Make the input graph satisfy triangular inequality constraint */
	triangular_inequality();

	/* Memory allocation of clients_coordinate(j) */
	clients_coordinate = new double*[n_clients];
	for (int j = 0; j < n_clients; j++) {
		clients_coordinate[j] = new double[n_facilities];
		for (int i = 0; i < n_facilities; i++) {
			clients_coordinate[j][i] = 0.0;
		}
	}

	/* Memory allocation of the exponential clocks of the facilities(Z) */
	exponential_clocks = new double[n_facilities];

	/* Memory allocation of opening_table, connection_table (Output of rounding alg) */
	opening_table = new bool[n_facilities];
	connection_table = new bool*[n_facilities];
	for (int i = 0; i < n_facilities; i++) {
		connection_table[i] = new bool[n_clients];
		for (int j = 0; j < n_clients; j++) {
			connection_table[i][j] = false;
		}
	}

}

/* Deconstructor definition */
FacilityLocation::~FacilityLocation()
{
	/* Memory deallocation of opening_cost(f), connection_cost(d) */
	delete opening_cost;
	for (int i = 0; i < n_facilities; i++) {
		delete[] connection_cost[i];
	}
	delete[] connection_cost;

	/* Memory deallocation of clients_coordinate(j) */
	for (int i = 0; i < n_clients; i++) {
		delete[] clients_coordinate[i];
	}
	delete[] clients_coordinate;

	/* Memory deallocation of the expoential clocks of the facilities */
	delete exponential_clocks;

	/* Memory deallocation of opening_table, connection_table */
	delete opening_table;
	for (int i = 0; i < n_facilities; i++) {
		delete[] connection_table[i];
	}
	delete[] connection_table;
}

/* LP-solver */
double FacilityLocation::LP_solve()
{
	IloEnv env;

	/* set and initialize connection variables */
	/*
	IloNumVar * x = new IloNumVar[n_facilities * n_clients];
	//IloNumVar x[n_clients * n_facilities];
	for (int i = 0; i < n_facilities; ++i) {
	for (int j = 0; j < n_clients; ++j) {
	x[i*n_clients + j] = IloNumVar(env, 0, IloInfinity);
	}
	}
	*/
	IloNumVarArray x(env, (IloInt)(n_facilities * n_clients), (IloNum)0, IloInfinity);


	/* set and initialize opening variables */
	/*
	IloNumVar *y = new IloNumVar[n_facilities];
	//IloNumVar y[n_facilities];
	for (int i = 0; i < n_facilities; ++i)
	y[i] = IloNumVar(env, 0, IloInfinity);
	*/
	IloNumVarArray y(env, (IloInt)n_facilities, (IloNum)0, IloInfinity);


	/* set ranges of sums of connection variables (1 <= sum_i(x_ij) <= 1 for all j) */
	/*
	IloExpr * sum_expr = new IloExpr[n_clients];
	IloRange * sum_condition = new IloRange[n_clients];
	//IloExpr sum_expr[n_clients];
	//IloRange sum_condition[n_clients];
	for (int j = 0; j < n_clients; ++j) {
	sum_expr[j] = IloExpr(env);
	for (int i = 0; i < n_facilities; ++i) {
	sum_expr[j] += x[i*n_clients + j];
	}
	sum_condition[j] = (sum_expr[j] == 1.0);
	}
	*/
	IloExprArray sum_expr(env);
	IloRangeArray sum_condition(env);
	for (int j = 0; j < n_clients; ++j) {
		sum_expr.add(IloExpr(env));
		for (int i = 0; i < n_facilities; ++i) {
			sum_expr[j] += x[i*n_clients + j];
		}
		sum_condition.add(sum_expr[j] == 1.0);
	}



	/* set ranges of connection variables and opening variables (-inf <= x_ij - y_i <= 0 for all i, j) */
	/*
	IloRange * x_range = new IloRange[n_clients * n_facilities];
	//IloRange x_range[n_clients * n_facilities];
	for (int i = 0; i < n_facilities; ++i) {
	for (int j = 0; j < n_clients; ++j) {
	x_range[i * n_clients + j] = IloRange(env, -IloInfinity, 0);
	x_range[i * n_clients + j].setLinearCoef(x[i*n_clients + j], 1);
	x_range[i * n_clients + j].setLinearCoef(y[i], -1);
	}
	}
	*/
	IloRangeArray x_range(env);
	for (int i = 0; i < n_facilities; ++i) {
		for (int j = 0; j < n_clients; ++j) {
			x_range.add(IloRange(env, -IloInfinity, 0));
			x_range[i * n_clients + j].setLinearCoef(x[i*n_clients + j], 1);
			x_range[i * n_clients + j].setLinearCoef(y[i], -1);
		}
	}

	/* set the obj fct (minimize sum_i(y_i*f_i) + sum_{i,j}(x_ij*d(i,j)) )*/
	IloObjective obj = IloMinimize(env, 0);
	for (int i = 0; i < n_facilities; ++i) {
		obj.setLinearCoef(y[i], this->opening_cost[i]);
		for (int j = 0; j < n_clients; ++j) {
			obj.setLinearCoef(x[i*n_clients + j], this->connection_cost[i][j]);
		}
	}

	/* compile the model */
	/*
	IloModel model(env);
	for (int j = 0; j < n_clients; ++j) {
	//model.add(sum_range[j]);
	model.add(sum_condition[j]);
	for (int i = 0; i < n_facilities; ++i) {
	model.add(x_range[i*n_clients + j]);
	}
	}
	*/
	IloModel model(env);
	model.add(sum_condition);
	model.add(x_range);
	model.add(obj);

	/* solve the model */
	IloCplex solver(model);
	try {
		solver.solve();
	}
	catch (IloException &ex) {
		cerr << ex << endl;
	}
	/* save results*/
	for (int i = 0; i < n_facilities; ++i) {
		for (int j = 0; j < n_clients; ++j) {
			this->clients_coordinate[j][i] = solver.getValue(x[i*n_clients + j]);
		}
	}
	return solver.getObjValue();

}

/* Rounding algorithm */
void FacilityLocation::rounding_alg()
{
	double r = (double)rand() / RAND_MAX;

	if (r <= 2.0 / 3.0) {
		rounded_cost = rounding_alg_exp();
	}
	else {
		rounded_cost = rounding_alg_dist();
	}
}

double FacilityLocation::rounding_alg_exp()
{
	/* Initialize exponential clocks */
	std::default_random_engine generator;  // memory leakage
	for (int i = 0; i < n_facilities; i++) {
		std::exponential_distribution<double> distribution(1.0);
		exponential_clocks[i] = distribution(generator);
	}

	/* Traversing all clients on simplex, find minimum Zi/ui. Then assign the client to ith facility */
	for (int j = 0; j < n_clients; j++) {
		double min = DBL_MAX;
		int selected_facility;
		for (int i = 0; i < n_facilities; i++) {
			if (CompareDoubleUlps(exponential_clocks[i] / clients_coordinate[j][i], min) == -1) {
				min = exponential_clocks[i] / clients_coordinate[j][i];
				selected_facility = i;
			}
		}
		connection_table[selected_facility][j] = true;
	}

	/* If at least one of all clients are connected ith facility, open that facility */
	for (int i = 0; i < n_facilities; i++) {
		bool ifOpen = false;
		for (int j = 0; j < n_clients; j++) {
			if (connection_table[i][j] == true) {
				ifOpen = true;
			}
		}
		opening_table[i] = ifOpen;
	}

	/* Calculate cost */

	double result = 0.0;

	for (int i = 0; i < n_facilities; i++) {
		if (opening_table[i] == true)
			result += opening_cost[i];
	}

	for (int i = 0; i < n_facilities; i++) {
		for (int j = 0; j < n_clients; j++) {
			if (connection_table[i][j] == true)
				result += connection_cost[i][j];
		}
	}

	return result;
}

double FacilityLocation::rounding_alg_dist()
{
	double r = (double)rand() / RAND_MAX;
	int *assigned_facility = new int[n_clients];
	
	/* Initial assigned facility is -1 (None) */
	for (int j = 0; j < n_clients; j++) {
		assigned_facility[j] = -1;
	}

	for (int i = 0; i < n_facilities - 1; i++) {
		for (int j = 0; j < n_clients; j++) {
			if (assigned_facility[j] == -1 && r < clients_coordinate[j][i] * clients_coordinate[j][i]) {
				assigned_facility[j] = i;
			}
		}
	}

	for (int j = 0; j < n_clients; j++) {
		if (assigned_facility[j] == -1) {
			assigned_facility[j] = n_facilities - 1;
		}
	}

	for (int i = 0; i < n_facilities; i++) {
		opening_table[i] = false;
		for (int j = 0; j < n_clients; j++) {
			if (assigned_facility[j] == i) {
				connection_table[i][j] = true;
				opening_table[i] = true;
			}
		}
	}

	/* Calculate cost */
	double result = 0.0;

	for (int i = 0; i < n_facilities; i++) {
		if (opening_table[i] == true)
			result += opening_cost[i];
	}

	for (int i = 0; i < n_facilities; i++) {
		for (int j = 0; j < n_clients; j++) {
			if (connection_table[i][j] == true)
				result += connection_cost[i][j];
		}
	}

	return result;
}

int minDistance(double dist[], bool sptSet[], unsigned int V)
{
	// Initialize min value
	double min = DBL_MAX;
	int min_index;

	for (int v = 0; v < V; v++)
		if (sptSet[v] == false && dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

void FacilityLocation::triangular_inequality(void) {
	/* Redefine the bipartite graph */
	double **graph = new double*[n_clients + n_facilities];
	for (int i = 0; i < n_clients + n_facilities; ++i) {
		graph[i] = new double[n_clients + n_facilities];
		for (int j = 0; j < n_clients + n_facilities; ++j) {
			graph[i][j] = 1e+200;
		}
	}
	/* draw the graph */
	for (int i = 0; i < n_facilities; ++i) {
		for (int j = 0; j < n_clients; ++j) {
			graph[i + n_clients][j] = connection_cost[i][j];
			graph[j][i + n_clients] = connection_cost[i][j];
		}
	}
	/* for each client */
	for (int j = 0; j < n_clients; ++j) {
		/* get shortest paths for each facility */
		int src = j;
		double *dist = new double[n_clients + n_facilities];     // The output array.  dist[i] will hold the shortest
																 // distance from src to i

		bool *sptSet = new bool[n_clients + n_facilities]; // sptSet[i] will true if vertex i is included in shortest
														   // path tree or shortest distance from src to i is finalized

														   // Initialize all distances as INFINITE and stpSet[] as false
		for (int i = 0; i < n_clients + n_facilities; i++)
			dist[i] = DBL_MAX, sptSet[i] = false;

		// Distance of source vertex from itself is always 0
		dist[src] = 0;

		// Find shortest path for all vertices
		for (int count = 0; count < n_clients + n_facilities - 1; count++)
		{
			// Pick the minimum distance vertex from the set of vertices not
			// yet processed. u is always equal to src in first iteration.
			int u = minDistance(dist, sptSet, n_clients + n_facilities);

			// Mark the picked vertex as processed
			sptSet[u] = true;

			// Update dist value of the adjacent vertices of the picked vertex.
			for (int v = 0; v < n_clients + n_facilities; v++)

				// Update dist[v] only if is not in sptSet, there is an edge from 
				// u to v, and total weight of path from src to  v through u is 
				// smaller than current value of dist[v]

				if (!sptSet[v] && graph[u][v] && CompareDoubleUlps(dist[u], DBL_MAX) != 0
					&& CompareDoubleUlps(dist[u] + graph[u][v], dist[v]) < 0)
					dist[v] = dist[u] + graph[u][v];
		}

		/* assign min(dist[u], e(i, src)) to connection_cost[i][src] for each facility i */
		for (int i = n_clients; i < this->n_facilities + this->n_clients; ++i) {
			if (dist[i] < this->connection_cost[i - n_clients][src])
				this->connection_cost[i - n_clients][src] = dist[i];
		}
		delete dist;
		delete sptSet;

	}
	for (int i = 0; i < n_clients + n_facilities; ++i) {
		delete graph[i];
	}
	delete graph;
}