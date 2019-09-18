#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "json.hxx"
#include <stdexcept> //For error-throwing

#include "RateCoefficient.hxx"
#include "sharedFunctions.hxx"

#include <algorithm> //for upper/lower_bound
#include <typeinfo>

using namespace std;

// for convenience
using json = nlohmann::json;

RateCoefficient::RateCoefficient(const string& filename){
	// # Create an instance of RateCoefficient by reading an OpenADAS JSON file
	
	json data_dict = retrieveFromJSON(filename);

	atomic_number   = data_dict["charge"];
	element         = data_dict["element"];
	adf11_file      = filename;

	vector<vector< vector<double> > > extract_log_coeff = data_dict["log_coeff"];
	vector<double> extract_log_temperature = data_dict["log_temperature"];
	vector<double> extract_log_density = data_dict["log_density"];
	// Doing this as a two-step process - since the first is casting JSON data into the stated type.
	// The second copies the value to the corresponding RateCoefficient attribute
	log_coeff = extract_log_coeff;
	log_temperature = extract_log_temperature;
	log_density = extract_log_density;

};
ostream& operator<<(ostream& os, const RateCoefficient& RC){  
    os << "RateCoefficient object from " << RC.adf11_file << endl;
    return os;  
}
double RateCoefficient::call0D(const int k, const double eval_Te, const double eval_Ne){

	// """Evaluate the ionisation/recombination coefficients of
	// 	k'th atomic state at a given temperature and density.
	// 	Args:
	// 		k  (int): Ionising or recombined ion stage,
	// 			between 0 and k=Z-1, where Z is atomic number.
	// 		Te (double): Temperature in [eV].
	// 		ne (double): Density in [m-3].
	// 	Returns:
	// 		c (double): Rate coefficent in [m3/s].
	
	// Perform a basic interpolation based on linear distance
	// values to search for
  
  double eval_log10_Te = log10(eval_Te < 1e-5 ? 1e-5 : eval_Te);
  double eval_log10_Ne = log10(eval_Ne < 1e-5 ? 1e-5 : eval_Ne);
	// Look through the log_temperature and log_density attributes of RateCoefficient to find nearest (strictly lower)
	// Subtract 1 from answer to account for indexing from 0
	int low_Te = lower_bound(log_temperature.begin(), log_temperature.end(), eval_log10_Te) - log_temperature.begin() - 1;
	int low_Ne = lower_bound(log_density.begin(), log_density.end(), eval_log10_Ne) - log_density.begin() - 1;

	// Bounds checking -- make sure you haven't dropped off the end of the array
        // An easy error to make is supplying the function arguments already having taken the log10
        if (low_Te >= static_cast<int>(log_temperature.size())-1) {
          // Te out of bounds on high side

          if (!warned_te_range) {
            // Print warning the first time this occurs
            std::cerr << "WARNING (Atomicpp::RateCoefficient): log Te too high (" <<  eval_log10_Te << " > " << *log_temperature.rbegin() << ")\n";
            std::cerr << "Te, Ne: " << eval_Te << ", " << eval_Ne << endl;
            warned_te_range = true;
          }
          eval_log10_Te = *log_temperature.rbegin(); // Last element
          low_Te = log_temperature.size()-2;
          
        } else if (low_Te <= -1) {
          // Te out of bounds on low side

          if (!warned_te_range) {
            std::cerr << "WARNING (Atomicpp::RateCoefficient): log Te too low (" <<  eval_log10_Te << " < " << *log_temperature.begin() << ")\n";
            std::cerr << "Te, Ne: " << eval_Te << ", " << eval_Ne << endl;
            warned_te_range = true;
          }
          eval_log10_Te = *log_temperature.begin();
          low_Te = 0;
        }
        
	if (low_Ne >= static_cast<int>(log_density.size())-1) {
          // Ne out of bounds on high side
          if (!warned_ne_range) {
            std::cerr << "WARNING (Atomicpp::RateCoefficient): log Ne too high (" <<  eval_log10_Ne << " > " << *log_density.rbegin() << ")\n";
            std::cerr << "Te, Ne: " << eval_Te << ", " << eval_Ne << endl;
            warned_ne_range = true;
          }
          eval_log10_Ne = *log_density.rbegin(); // Last element
          low_Ne = log_density.size()-2;
          
        } else if (low_Ne <= -1) {
          // Ne out of bounds on low side
          if (!warned_ne_range) {
            std::cerr << "WARNING (Atomicpp::RateCoefficient): log Ne too low (" <<  eval_log10_Ne << " < " << *log_density.begin() << ")\n";
            std::cerr << "Te, Ne: " << eval_Te << ", " << eval_Ne << endl;
            warned_ne_range = true;
          }
          eval_log10_Ne = *log_density.begin();
          low_Ne = 0;
	}

	int high_Te = low_Te + 1;
	int high_ne = low_Ne + 1;

	double Te_norm = 1/(log_temperature[high_Te] - log_temperature[low_Te]); //Spacing between grid points
	double ne_norm = 1/(log_density[high_ne] - log_density[low_Ne]); //Spacing between grid points

	double x = (eval_log10_Te - log_temperature[low_Te])*Te_norm;
	double y = (eval_log10_Ne - log_density[low_Ne])*ne_norm;
	
	// // Construct the simple interpolation grid
	// // Find weightings based on linear distance
	// // w01 ------ w11    ne -> y
	// //  | \     / |      |
	// //  |  w(x,y) |    --/--Te -> x
	// //  | /     \ |      |
	// // w00 ------ w10  

	double eval_log10_coeff =
	(log_coeff[k][low_Te][low_Ne]*(1-y) + log_coeff[k][low_Te][high_ne]*y)*(1-x)
	+(log_coeff[k][high_Te][low_Ne]*(1-y) + log_coeff[k][high_Te][high_ne]*y)*x;

	double eval_coeff = pow(10,eval_log10_coeff);
	return eval_coeff;
};
int RateCoefficient::get_atomic_number(){
	return atomic_number;
};
string RateCoefficient::get_element(){
	return element;
};
string RateCoefficient::get_adf11_file(){
	return adf11_file;
};
vector<vector< vector<double> > > RateCoefficient::get_log_coeff(){
	return log_coeff;
};
vector<double> RateCoefficient::get_log_temperature(){
	return log_temperature;
};
vector<double> RateCoefficient::get_log_density(){
	return log_density;
};
