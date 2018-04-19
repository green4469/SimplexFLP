#include "FacilityLocation.h"

int main(int argc, char *argv[])
{
	FacilityLocation *FL = new FacilityLocation(argc, argv);

	cout << "LP-solution: " << FL->LP_solve() << endl;
	FL->rounding_alg();
	cout << "Rounded solution: " << FL->rounded_cost << endl;

	delete FL;
	return 0;
}