#include <iostream>
#include <string>
#include <algorithm>

#include "Case.hpp"

int main(int argn, char **args) {

    /*
    std::string option;
    std::vector<std::string> cases = {"ChannelWithBFS", "ChannelWithObstacle", "FluidTrap", "LidDrivenCavity", "NaturalConvection", "RayleighBenard"};

    std::cout << "Please choose one of the following options to simulate a given problem:" << '\n';
    std::cout << "1) " << cases[1] << "." << '\n';
    std::cout << "2) " << cases[2] << "." << '\n';
    std::cout << "3) " << cases[3] << "." << '\n';
    std::cout << "4) " << cases[4] << "." << '\n';
    std::cout << "5) " << cases[5] << "." << '\n';
    std::cout << "6) " << cases[6] << "." << '\n';
    
    std::cin >> option;
    while ((option != "1") and (option != "2") and (option != "3") and (option != "4") and (option != "5") and (option != "6")){
    
    	std::cout << "Please choose a valid option among the available ones:" << '\n';
    	std::cin >> option;
    }
    
    int idx_option = std::stoi(option);
    
    std::string file_name{cases[idx_option]};
    */
    
    std::vector<std::string> cases = {"ChannelWithBFS", "ChannelWithObstacle", "FluidTrap", "LidDrivenCavity", "NaturalConvection", "RayleighBenard"};
    std::string file_name;
    if (argn > 1) {
    	for (int i = 1; i < argn; ++i){
		if (std::find(cases.begin(), cases.end(), args[i]) != cases.end()) {
			std::cout << "Initializing simulation for case: " << args[i] << '\n';
			file_name = "../example_cases/";
			file_name = file_name + args[i] + "/" + args[i] + ".dat";
			Case problem(file_name);
			std::cout << "Created new case. Starting Simulation:\n";
			problem.simulate();
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

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen case_name" << std::endl;
    }
    
}
