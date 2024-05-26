#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <random>
#include "struc.hpp"
#include "consts.hpp"

double random_number(double r_down,double r_up){
    double r = (double)(rand())/RAND_MAX;
    r = r_down + r*(r_up-r_down);
    return r;
}

Atom::Atom(double X,double Y, double Z){
	coord.x = X;
	coord.y = Y;
	coord.z = Z;
	force.x = 0.0;//random_number(-0.1,0.1);
	force.y = 0.0;//random_number(-0.1,0.1);
	force.z = 0.0;//random_number(-0.1,0.1);
	velocity.x = random_number(-1,1)*430; // 430 for m/s
	velocity.y = random_number(-1,1)*430; // 430 for m/s
	velocity.z = random_number(-1,1)*430; // 430 for m/s
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


Simulation::Simulation(int nAt){
    // This is horrible...
    int count = 0;
    int n_per_dim = 3;
    int x_move = 0;
    int y_move = 0;
    int z_move = 0;
    bool place = true;
    std::vector<Atom> temp_mol;
    while (place){
        Atom temp = Atom(0.1+4*x_move,0.1+4*y_move,0.1+4*z_move);
	temp_mol.push_back(temp);
	x_move++;
	if (x_move == n_per_dim){
	    x_move = 0;
	    y_move++;
	    if (y_move == n_per_dim){
		y_move = 0;
		z_move++;
	    }
	}	
	count++;
	if (count == nAt){
	    place = false;
	}
    }
    mol = temp_mol;
}

Simulation::Simulation(int nAt,double x_box,double y_box, double z_box){
    // This is horrible...
    int count = 0;
    int n_per_dim = 3;
    int x_move = 0;
    int y_move = 0;
    int z_move = 0;
    bool place = true;
    std::vector<Atom> temp_mol;
    while (place){
        Atom temp = Atom(0.1+4*x_move,0.1+4*y_move,0.1+4*z_move);
	temp_mol.push_back(temp);
	x_move++;
	if (x_move == n_per_dim){
	    x_move = 0;
	    y_move++;
	    if (y_move == n_per_dim){
		y_move = 0;
		z_move++;
	    }
	}	
	count++;
	if (count == nAt){
	    place = false;
	}
    }
    mol = temp_mol;
    PBC = true;
    box.x = x_box;
    box.y = y_box;
    box.z = z_box;
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


void Simulation::rescale_v(double T){
    double moment = 0;
    for (auto atom = mol.begin();atom < mol.end() ;atom++){
	double vel = std::sqrt(std::pow(atom->velocity.x,2)+std::pow(atom->velocity.y,2)+std::pow(atom->velocity.z,2));
	moment += vel; 
    }
    std::cout << "moment: " << moment << std::endl;
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
	auto eps1 = element_epsilons.at(current_el)*eV;
	auto sigma1 = element_sigma.at(current_el)*A_to_m;
	for (auto iter2 = iter+1; iter2 < mol.end() ; iter2++){
	    auto current_el = iter2->element;
	    auto eps2 = element_epsilons.at(current_el)*eV;
	    auto sigma2 = element_sigma.at(current_el)*A_to_m;
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
	auto epsilon1 = element_epsilons.at(current_el)*eV;
	auto sigma1 = element_sigma.at(current_el)*A_to_m;
	
	for (auto iter2 = iter+1; iter2 < sim->mol.end() ; iter2++){
	    auto current_el = iter2->element;
	    auto epsilon2 = element_epsilons.at(current_el)*eV;
	    auto sigma2 = element_sigma.at(current_el)*A_to_m;
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

