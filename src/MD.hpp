#include <vector>
#include "struc.hpp"


std::vector<Atom> read_struc();
void write_struc(Simulation sim,std::string comment);
void take_step(Simulation *sim,double dt);
void update_velocity(Simulation *sim,double dt);