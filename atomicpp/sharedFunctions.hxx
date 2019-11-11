#ifndef SHAREDFUNCTIONS_H //Preprocessor directives to prevent multiple definitions
#define SHAREDFUNCTIONS_H

	#include <string>
	#include <fstream>

	#include "json.hxx"
	using json = nlohmann::json;

	using namespace std; //saves having to prepend std:: onto common functions
	// for convenience

	/**
	 * @brief Returns a JSON object from a specified file location
	 * @details Reads a .json file given at path_to_file
	 * Uses the json module at https://github.com/nlohmann/json/
	 * This is supplied as a header file "json.hpp" which must be included in the same folder as the source
	 * N.b. do NOT pass path_to_file by reference - results in error (Undefined symbols)!
	 * @param path_to_file string, giving relative or absolute path to file
	 * @return JSON object - basically a map/dictionary which returns the stored value associated with a given (string) key
	 * The JSON object tries to convert data from a custom class to whatever you ask it to assign to - most of the time
	 * it works as intuition would suggest
	 */
	json retrieveFromJSON(string path_to_file);
	/**
	 * @brief A simple program to check if a file exists
	 * 
	 * @param name string, giving relative or absolute path to file
	 * @return boolean - true if file found, false if not
	 */
	bool test_file_exists (const string& name);
#endif
