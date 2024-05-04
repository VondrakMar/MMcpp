#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include "struc.hpp"




int main(){
    auto mol = read_struc();
    for (int i = 0; i <=100; i++){
	std::cout << "#################################\n";
	LJ_forces(&mol);
	for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	std::cout << iter->force.x << " " << iter->force.y << " " << iter->force.z << std::endl;
	}
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
    
