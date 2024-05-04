#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include "struc.hpp"

// Units are in eV, and A
// parameteres taken from https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405


int main(){
    auto a1 = Atom(1.0,1.0,1.0);
    auto a2 = Atom(1.5,1.3,1.2);
    auto a3 = Atom(1.5,0.8,0.6);
    auto a4 = Atom(1.2,0.3,2.3);
    std::vector<Atom> mol;
    mol.push_back(a1);
    mol.push_back(a2);
    mol.push_back(a3);
    mol.push_back(a4);
    // std::cout << "final_energy " << energy << std::endl;
    for (int i = 0; i <=10000; i++){
	// auto energy = LJ_energy(mol);

	// std::cout << energy << "hello world\n";
	for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	    auto rf = random_forces();
	    iter->coord.x += rf.x;
	    iter->coord.y += rf.y;
	    iter->coord.z += rf.z;	    
	}
	auto energy = LJ_energy_mol(mol);
	write_struc(mol,std::to_string(energy));
    }
    
    return 0;
}
