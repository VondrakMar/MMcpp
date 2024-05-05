#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include "struc.hpp"
#include <random>


double random_number(double r_down,double r_up){
    //double r = -r_down + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(r_up-(-r_down))));
    double r = (double)(rand())/RAND_MAX;
    r = r_down + r*(r_up-r_down);
    return r;
}

Atom::Atom(double X,double Y, double Z){
	coord.x = X;
	coord.y = Y;
	coord.z = Z;
	force.x = random_number(-0.1,0.1);
	force.y = random_number(-0.1,0.1);
	force.z = random_number(-0.1,0.1);
	velocity.x = random_number(-1,1);
	velocity.y = random_number(-1,1);
	velocity.z = random_number(-1,1);
	element = 18;
}

void Atom::print_coord(){
    std::cout << this->coord.x << " " << this->coord.y << " " << this->coord.z << "\n"; 
}
void Atom::print_forces(){
    std::cout << this->force.x << " " << this->force.y << " " << this->force.z << "\n"; 
}
void Atom::print_velocity(){
    std::cout << this->velocity.x << " " << this->velocity.y << " " << this->velocity.z << "\n"; 
}


Simulation::Simulation(std::vector<Atom> Mol){
    mol = Mol;
    PBC = false;
}
Simulation::Simulation(std::vector<Atom> Mol,double x_box,double y_box, double z_box){
    mol = Mol;
    PBC = true;
    box.x = x_box;
    box.y = y_box;
    box.z = z_box;
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




double pair_distance(vec3D a1,vec3D a2){
    double r2 = std::pow(a1.x-a2.x,2)+std::pow(a1.y-a2.y,2)+std::pow(a1.z-a2.z,2);
    double r = std::sqrt(r2);
    return r;
}

double LJ_energy_pair(double d,double eps,double sigma){
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

void LJ_forces(Simulation *sim){
    // this has to be changed omg
    for (auto iter = sim->mol.begin(); iter < sim->mol.end() ; iter++){
	iter->force.x = 0.0;
	iter->force.y = 0.0;
	iter->force.z = 0.0;
    }
    for (auto iter = sim->mol.begin(); iter < sim->mol.end() ; iter++){
	auto current_el = iter->element;
	auto epsilon1 = element_epsilons.at(current_el);
	auto sigma1 = element_sigma.at(current_el);
	
	for (auto iter2 = iter+1; iter2 < sim->mol.end() ; iter2++){
	    auto current_el = iter2->element;
	    auto epsilon2 = element_epsilons.at(current_el);
	    auto sigma2 = element_sigma.at(current_el);
	    auto epsilon_comb = std::sqrt(epsilon1*epsilon2);
	    auto sigma_comb = std::sqrt(sigma1*sigma2);
	    double dx = iter2->coord.x - iter->coord.x;
            double dy = iter2->coord.y - iter->coord.y;
            double dz = iter2->coord.z - iter->coord.z;
	        double r = sqrt(dx*dx + dy*dy + dz*dz);
            double r_inv = sigma_comb / r;
            double r_inv6 = pow(r_inv, 6);
            double r_inv12 = r_inv6 * r_inv6;
            double force_mag = 24 * epsilon_comb * (2 * r_inv12 - r_inv6) * r_inv;
            iter->force.x -= force_mag * dx * r_inv;
            iter->force.y -= force_mag * dy * r_inv;
            iter->force.z -= force_mag * dz * r_inv;
            iter2->force.x += force_mag * dx * r_inv;
            iter2->force.y += force_mag * dy * r_inv;
            iter2->force.z += force_mag * dz * r_inv;
	}
    }
}

