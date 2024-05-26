#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "utils.hpp"
#include "consts.hpp"

void write_struc(Simulation sim,std::string file_name,std::string comment){
    double dist_factor = m_to_A;
    std::ofstream outFile(file_name,std::ios::app);
    outFile << sim.mol.size() << std::endl;
    outFile << comment << std::endl;
    for (auto iter = sim.mol.begin(); iter < sim.mol.end() ; iter++){
	    outFile << "Ar" << " " << iter->coord.x*dist_factor << " " << iter->coord.y*dist_factor << " " << iter->coord.z*dist_factor << " "    << std::endl;	
    }
    outFile.close();
}

std::vector<Atom> read_struc(){
    double dist_factor = A_to_m; 
    std::vector<Atom> mol;
    std::ifstream inFile("small.xyz");
    std::string line;
    getline(inFile,line);
    int nAts = std::stoi(line);
    std::cout << "Loading a structure with "<< nAts << " atoms pog"  <<std::endl;
    getline(inFile,line); // now just for reading the line and doing nothing, this should be change for extended xyz
    std::cout << "I am doing nothing with this info " << line << std::endl;
    while (getline(inFile,line)){
	    std::istringstream iss(line);
	    std::string temp;
	    double tx,ty,tz;
	    iss >> temp; // skiping first column, this has to be change when more than arg
	    iss >> tx >> ty >> tz;
        tx = dist_factor*tx;
        ty = dist_factor*ty;
        tz = dist_factor*tz;
	    Atom tempAt = Atom(tx,ty,tz);
	    mol.push_back(tempAt);
    }
    return mol;
    
}


void read_input(Simulation *sim,std::string file_name){
    std::ifstream inFile(file_name);
    std::string line;
    std::cout << "Loading input file\n";
    while (getline(inFile,line)){
        std::istringstream iss(line);
        std::string keyword;
        iss >> keyword;
        if (keyword.compare("pbc") == 0){
            iss >> keyword;
            if (keyword.compare("true") == 0){
                sim->PBC=true;
            }
        else{
            sim->PBC=false;
        }
        }
    }
}