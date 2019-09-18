#include <map>
#include <string>
#include <iostream>
#include <stdexcept> //For error-throwing
#include <memory> //For smart pointers
#include <set>

#include "json.hxx"
using json = nlohmann::json;
#include "ImpuritySpecies.hxx"
#include "sharedFunctions.hxx"
#include "RateCoefficient.hxx"

using namespace std;

ImpuritySpecies::ImpuritySpecies(string& impurity_symbol_supplied){
	// Constructor for the ImpuritySpecies class
	// Input is a string, which typically should be of length 1 (i.e. "c" for carbon, "n" for nitrogen)
	// Will return an ImpuritySpecies object
	
	// Convert impurity_symbol_supplied to lowercase (this is what impurity_user_input uses as a key)
	impurity_symbol_supplied[0] = tolower(impurity_symbol_supplied[0]);

	// Declare the path to the JSON file which contains hard-coded parameters for the ADAS data corresponding to the impurities
	// string user_file = get_impurity_user_input();
	string user_file = "impurity_user_input.json";

	// Declare the path the directory which contains the JSON files created from OpenADAS .dat files
	// To convert .dat to .json, see https://github.com/TBody/OpenADAS_to_JSON (should be included as a subdirectory of SD1D)
	// string json_database_path = get_json_database_path();
	string json_database_path = "json_database/";

	// Define a map between a recognisable key for various physics process, and the 3-letter code used by OpenADAS to represent this data
	map<string,string>datatype_abbrevs={
		{"ionisation",           "scd"},
		{"recombination",        "acd"},
		{"cx_recc",              "ccd"},
		{"continuum_power",      "prb"},
		{"line_power",           "plt"},
		{"cx_power",             "prc"},
		{"ionisation_potential", "ecd"} //N.b. ionisation_potential is not a rate-coefficient, but most of the methods are transferable
	};

	// Save from impurity_symbol_supplied and the hard-coded data within user_file
	symbol              = impurity_symbol_supplied;
	// retrieve from JSON and set attributes
	json j_object = retrieveFromJSON(user_file);
	auto check_symbol_in_file = j_object.find(impurity_symbol_supplied);
	if ((check_symbol_in_file != j_object.end())){
		name                = j_object[impurity_symbol_supplied]["name"];
		year                = j_object[impurity_symbol_supplied]["year"];
		has_charge_exchange = j_object[impurity_symbol_supplied]["has_charge_exchange"];
		atomic_number       = j_object[impurity_symbol_supplied]["atomic_number"];
	} else {
		cout << "Error: "<< impurity_symbol_supplied << " not found in " << user_file << "\n";
		throw invalid_argument( "Key not found in user input file" );
	};

	// Add the JSON files associated with this impurity to its .adas_files_dict attribute where the key is the (extended) process name, which maps to a filename (string)
	// Mark the physics processes which rely on charge exchange
	set<string> cx_processes{datatype_abbrevs["cx_recc"],datatype_abbrevs["cx_power"]};
	// Iterate over the key -> value pairs in datatype_abbrevs
	for (auto& kv : datatype_abbrevs) {
		string physics_process = kv.first;
		string filetype_code = kv.second;
		// Check to see if physics_process is a charge exchange process
		bool code_is_cx = find(begin(cx_processes), end(cx_processes), filetype_code) != end(cx_processes);
		// If it is, and the element doesn't have charge exchange data, then skip this process
		if (get_has_charge_exchange() or not(code_is_cx)){
			// Otherwise, add the filenames of the json files corresponding to the rate-coefficients for the physics process to .adas_files_dict
		    addJSONFiles(physics_process, filetype_code, json_database_path);
		};
	}

	// # Use the .adas_file_dict files to generate RateCoefficient objects for each process
	// # Uses the same keys as .adas_file_dict
	makeRateCoefficients();

};
void ImpuritySpecies::addJSONFiles(const string& physics_process, const string& filetype_code, const string& json_database_path){
	// # 1. Make the filename string expected for the json adas file
	// # 2. Check that this file exists in the JSON_database_path/json_data directory
	// # 3. Add this file to the atomic data .adas_files_dict attribute
	string filename;
	string year_to_string = to_string(year);

	filename = json_database_path + "/" + filetype_code + year_to_string.substr(2,4) + "_" + symbol + ".json";

	if (not(test_file_exists(filename))) {
		cout << boolalpha;
		cout << "\nFile found: " << test_file_exists(filename) << "\n";
		throw runtime_error("File not found error");
	}

	adas_files_dict[physics_process] = filename;
};
void ImpuritySpecies::makeRateCoefficients(){
	// # Calls the RateCoefficient constructor method for each entry in the .adas_files_dict
	// # Generates a dictionary smart pointer to RateCoefficient objects as .rate_coefficients
	// See http://umich.edu/~eecs381/handouts/C++11_smart_ptrs.pdf for information on smart pointers (memory-managed)

	for (auto& kv : get_adas_files_dict()) {
		string physics_process = kv.first;
		string filename = kv.second;
		// Make a new RateCoefficient object by calling the RateCoefficient constructor on 'filename'
		// Create a smart pointer 'RC' that points to this object
		shared_ptr<RateCoefficient> RC(new RateCoefficient(filename));
		// Add 'RC' to the rate_coefficients attribute of ImpuritySpecies
		// (n.b. this is a map from a string 'physics_process' to a smart pointer which points to a RateCoefficient object)
		rate_coefficients[physics_process] = RC;
	}
};
// Accessor functions
	string ImpuritySpecies::get_symbol(){
		return symbol;
	};
	string ImpuritySpecies::get_name(){
		return name;
	};
	int ImpuritySpecies::get_year(){
		return year;
	};
	bool ImpuritySpecies::get_has_charge_exchange(){
		return has_charge_exchange;
	};
	int ImpuritySpecies::get_atomic_number(){
		return atomic_number;
	};
	map<string,string> ImpuritySpecies::get_adas_files_dict(){
		return adas_files_dict;
	};
	map<string,shared_ptr<RateCoefficient> > ImpuritySpecies::get_rate_coefficients(){
		return rate_coefficients;
	};
	shared_ptr<RateCoefficient> ImpuritySpecies::get_rate_coefficient(const string& key){
		return rate_coefficients[key];
	};
// Accessing environment variables (shared by any function which calls the ImpuritySpecies.hpp header) -- shared functions
	string get_json_database_path() {
		string json_database_env = "ADAS_JSON_PATH";
		char * environment_variable;
		environment_variable = getenv( json_database_env.c_str() );

		if (environment_variable != NULL) {
		    return environment_variable;
		} else {
			cout << "ADAS_JSON_PATH not found. Setting to default (= json_database/json_data)" << endl;
			return "json_database/json_data";
		}
	}
	string get_impurity_user_input() {
		string json_database_env = "ADAS_JSON_IMPURITY";
		char * environment_variable;
		environment_variable = getenv( json_database_env.c_str() );

		if (environment_variable != NULL) {
		    return environment_variable;
		} else {
			cout << "ADAS_JSON_IMPURITY not found. Setting to default (= c)" << endl;
			return "json_database/json_data";
	}
}
