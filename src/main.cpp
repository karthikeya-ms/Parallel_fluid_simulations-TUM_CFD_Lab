#include <iostream>
#include <string>
#include <algorithm>

#include "Case.hpp"

int main(int argn, char **args) {
    
    std::vector<std::string> cases = {"ChannelWithBFS", "ChannelWithObstacle", "FluidTrap", "LidDrivenCavity", "NaturalConvection", "RayleighBenard"};
    std::string file_name;
    if (argn > 1) {
    	for (int i = 1; i < argn; ++i){
		if (std::find(cases.begin(), cases.end(), args[i]) != cases.end()) {
			std::cout << "Initializing simulation for case: " << args[i] << '\n';
			file_name = "../example_cases/";
			file_name = file_name + args[i] + "/" + args[i] + ".dat";
			Case problem(file_name, argn, args);
			std::cout << "Created new case. Starting Simulation:\n";
			problem.simulate();
			std::cout << "Successfully finished simulation for case: " << args[i] << '\n';
		}
		else {
			std::cout << "Error: Invalid input file name. " << "'" << args[i] << "'" << " is not a valid simulation case." << '\n';
			std::cout << "Please provide one or more of the following option(s) as parameter(s) to simulate a given case and compile again accordingly:" << '\n';
			for (int i = 0; i < cases.size() - 1; ++i){
				std::cout << "- " << cases[i] << "." << '\n';
			}
			std::cout << "- " << cases[cases.size() - 1] << "." << '\n';
			std::cout << "For more information on how to run the code check READM.md file." << std::endl;  
		}
    	}
    }
}
