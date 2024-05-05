#include <vector>
#include "struc.hpp"


std::vector<Atom> read_struc();
void write_struc(std::vector<Atom> mol,std::string comment);
void take_step(std::vector<Atom> *mol,double dt);
void update_velocity(std::vector<Atom> *mol,double dt);