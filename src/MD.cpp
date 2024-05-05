#include <iostream>
#include <sstream>
#include <fstream>
#include "MD.hpp"

void take_step(std::vector<Atom> *mol,double dt){
    for (auto atom = mol->begin(); atom < mol->end() ; atom++){
	atom->coord.x += atom->velocity.x * dt;
	atom->coord.y += atom->velocity.y * dt;
	atom->coord.z += atom->velocity.z * dt;
    }
}

void update_velocity(std::vector<Atom> *mol,double dt){
    for (auto atom = mol->begin(); atom < mol->end() ; atom++){
	atom->velocity.x += atom->force.x * dt;
	atom->velocity.y += atom->force.y * dt;
	atom->velocity.z += atom->force.z * dt;
    }
}    

void write_struc(std::vector<Atom> mol,std::string comment){
    std::ofstream outFile("struc.xyz",std::ios::app);
    outFile << mol.size() << std::endl;
    outFile << comment << std::endl;
    for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	outFile << "Ar" << " " << iter->coord.x << " " << iter->coord.y << " " << iter->coord.z << " "    << std::endl;	
    }
    outFile.close();
}

std::vector<Atom> read_struc(){
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
	float tx,ty,tz;
	iss >> temp; // skiping first column, this has to be change when more than arg
	iss >> tx >> ty >> tz;
	Atom tempAt = Atom(tx,ty,tz);
	mol.push_back(tempAt);
    }
    return mol;
    
}