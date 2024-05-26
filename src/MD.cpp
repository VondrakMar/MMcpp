#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "MD.hpp"

// #include 



double PBC_handlere(double pos, double box) {
    double new_coor = std::fmod(pos, box);
    
    if (new_coor < 0) {
        new_coor += std::abs(box);
    }
    return new_coor;
}


void take_step(Simulation *sim,double dt){
    if (sim->PBC) {
        for (auto atom = sim->mol.begin(); atom < sim->mol.end() ; atom++){
	        atom->coord.x += atom->velocity.x * dt;
	        atom->coord.y += atom->velocity.y * dt;
	        atom->coord.z += atom->velocity.z * dt;
		// atom->coord.x = std::fmod(atom->coord.x,sim->box.x);
	        // atom->coord.y = std::fmod(atom->coord.y,sim->box.y);
	        // atom->coord.z = std::fmod(atom->coord.z,sim->box.z);
		atom->coord.x = PBC_handlere(atom->coord.x,sim->box.x);
	        atom->coord.y = PBC_handlere(atom->coord.y,sim->box.y);
	        atom->coord.z = PBC_handlere(atom->coord.z,sim->box.z);


	}
    }
    else {
        for (auto atom = sim->mol.begin(); atom < sim->mol.end() ; atom++){
	        atom->coord.x += atom->velocity.x * dt;
	        atom->coord.y += atom->velocity.y * dt;
	        atom->coord.z += atom->velocity.z * dt;
        }
    }
}

void update_velocity(Simulation *sim,double dt){
    for (auto atom = sim->mol.begin(); atom < sim->mol.end() ; atom++){
	atom->velocity.x += atom->force.x * dt;
	atom->velocity.y += atom->force.y * dt;
	atom->velocity.z += atom->force.z * dt;
    }
}    


