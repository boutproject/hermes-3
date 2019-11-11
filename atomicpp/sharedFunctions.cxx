#include <string>
#include <fstream>

#include "json.hxx"
using json = nlohmann::json;

using namespace std; //saves having to prepend std:: onto common functions
// for convenience

json retrieveFromJSON(string path_to_file){
	// Do not pass path_to_file by reference - results in error!
	// Reads a .json file given at path_to_file
	// Uses the json module at https://github.com/nlohmann/json/
	// This relies upon the "json.hpp" header which must be included in the same folder as the source
	
	// Open a file-stream at path_to_file
	ifstream json_file(path_to_file);
	// Initialise a json file object at j_object
	json j_object;
	json_file >> j_object;
	return j_object;
};
bool test_file_exists (const string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }   
};
