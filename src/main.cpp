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
    Simulation sim_run = Simulation(mol,3,3,3); 
	for (int i = 0; i <=500; i++){
	LJ_forces(&sim_run);
	update_velocity(&sim_run,0.01);
	take_step(&sim_run,0.01);

	write_struc(sim_run,"hello seamna");
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
    
