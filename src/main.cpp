#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include "utils.hpp"
#include "MD.hpp"


int main(){
    // auto mol = read_struc();
    // auto sim_run = Simulation(27);
// auto sim_run = Simulation(7,20,20,20);
    // sim_run.mol[0].coord.y = 10;
    // sim_run.rescale_v(300);
    auto mol = read_struc();
    write_struc(mol,"seamna.xyz","hello seaman");
    // Simulation sim_run = Simulation(mol,3.5,3.5,3.5);
    Simulation sim_run = Simulation(mol);
    read_input(&sim_run,"input.md");
    // for (int i = 0; i <=100; i++){
    // 	LJ_forces(&sim_run);
    // 	update_velocity(&sim_run,1e-15);
    // 	take_step(&sim_run,1e-15);
    // 	// write_struc(sim_run,"Lattice=\"20.0 0.0 0.0 0.0 20.0 0.0 0.0 0.0 20.0\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"");
    // 	write_struc(sim_run,"struc.xyz","why is this happening");
    // }
    return 0;
}
