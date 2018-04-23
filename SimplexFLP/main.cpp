#include "FacilityLocation.h"

string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
);

string replace_all(
	__in const std::string &message,
	__in const std::string &pattern,
	__in const std::string &replace
) {

	std::string result = message;
	std::string::size_type pos = 0;
	std::string::size_type offset = 0;

	while ((pos = result.find(pattern, offset)) != std::string::npos)
	{
		result.replace(result.begin() + pos, result.begin() + pos + pattern.size(), replace);
		offset = pos + replace.size();
	}

	return result;
}


int main(int argc, char *argv[])
{
	if (argc == 1) {
		int n_iters = 0;
		int dim = 0;
		int j = 0;

		cout << "몇 번 반복?: ";
		cin >> n_iters;

		cout << "n_facilities = n_clients = ?: ";
		cin >> dim;

		for (int i = 0; i < n_iters; i++) {
			std::srand(unsigned(std::time(NULL)) + ++j * 10);
			FacilityLocation *FL = new FacilityLocation(argc, argv, dim);
			cout << "LP-solution: " << FL->LP_solve() << endl;
			FL->rounding_alg();
			cout << "Rounded solution: " << FL->rounded_cost << endl;

			delete FL;
		}
	}

	if (argc == 2) {
		FacilityLocation *FL = new FacilityLocation(argc, argv);
		double sol = FL->LP_solve();
		cout << "LP-solution: " << sol << endl;
		FL->rounding_alg();
		cout << "Rounded solution: " << FL->rounded_cost << endl;

		string out_file;
		out_file = argv[1];
		out_file = replace_all(out_file, "IN", "OUT");
		ofstream out(out_file);
		out << sol << endl;
		out << FL->rounded_cost << endl;
		bool* opening_table = FL->opening_table;
		for (int i = 0; i < FL->n_facilities; i++) {
			if (opening_table[i] == true) {
				out << i << " ";
			}
		}
		out << endl;
		bool ** connection_table = FL->connection_table;
		for (int i = 0; i < FL->n_facilities; i++) {
			for (int j = 0; j < FL->n_clients; j++) {
				if (connection_table[i][j] == true) {
					out << i << "-" << j << endl;
				}
			}
		}



		delete FL;
	}
	return 0;
}