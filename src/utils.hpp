#include <vector>
#include "struc.hpp"

std::vector<Atom> read_struc();
void read_input(Simulation *sim,std::string file_name);
void write_struc(Simulation sim,std::string file_name,std::string comment);