#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
// #include "struc.hpp"
#include "MD.hpp"



int main(){
    auto mol = read_struc();
    Simulation sim_run = Simulation(mol); 
	for (int i = 0; i <=300; i++){
	LJ_forces(&mol);
	update_velocity(&mol,0.01);
	take_step(&mol,0.01);

	write_struc(mol,"hello seamna");
    }
    return 0;
}
	// auto energy = LJ_energy(mol);

	// std::cout << energy << "hello world\n";
	    // auto rf = random_forces();
	    // iter->coord.x += rf.x;
	    // iter->coord.y += rf.y;
	    // iter->coord.z += rf.z;	    
	// }
	// auto energy = LJ_energy_mol(mol);
	// write_struc(mol,std::to_string(energy));
    
