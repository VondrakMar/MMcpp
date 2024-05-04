#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include "struc.hpp"


float random_number(float r_down,float r_up){
    float r = -r_down + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(r_up-(-r_down))));
    return r;
}

Atom::Atom(double X,double Y, double Z){
	coord.x = X;
	coord.y = Y;
	coord.z = Z;
	force.x = random_number(-0.1,0.1);
	force.y = random_number(-0.1,0.1);
	force.z = random_number(-0.1,0.1);
	element = 18;
}

vec3D random_forces(){
    vec3D rf;
    float r_down = -0.2;
    float r_up = -0.2;
    rf.x = random_number(r_down,r_up);
    rf.y = random_number(r_down,r_up);
    rf.z = random_number(r_down,r_up);
    return rf;
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

double pair_distance(vec3D a1,vec3D a2){
    double r2 = std::pow(a1.x-a2.x,2)+std::pow(a1.y-a2.y,2)+std::pow(a1.z-a2.z,2);
    double r = std::sqrt(r2);
    return r;
}

double LJ_energy_pair(double d,double eps,double sigma){
    double V_ij = 4 * eps * (std::pow((sigma / d),12) - std::pow((sigma / d),6));
    return V_ij;
}
double LJ_forces(double d,double eps,double sigma){
    double V_ij = 4 * eps * (std::pow((sigma / d),12) - std::pow((sigma / d),6));
    return V_ij;
}

double LJ_energy_mol(std::vector<Atom> mol){
    double energy = 0;
    for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	auto current_el = iter->element;
	auto eps1 = element_epsilons.at(current_el);
	auto sigma1 = element_sigma.at(current_el);
	for (auto iter2 = iter+1; iter2 < mol.end() ; iter2++){
	    auto current_el = iter2->element;
	    auto eps2 = element_epsilons.at(current_el);
	    auto sigma2 = element_sigma.at(current_el);
	    auto eps_comb = std::sqrt(eps1*eps2);
	    auto sigma_comb = std::sqrt(sigma1*sigma2);
	    auto d = pair_distance(iter->coord,iter2->coord);
	    energy += LJ_energy_pair(d,eps_comb,sigma_comb);
	}
    }

    return energy;
}
