#include <vector>
#include <string>
#include <map>

// Units are in eV, and A
// parameteres taken from https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405
const std::map<int,float> element_epsilons = {
    {18,0.0103}
};
const std::map<int,float> element_sigma = {
    {18,3.4}
};

double random_number(double r_down,double r_up);

struct vec3D{
    double x,y,z;
};

class Atom{
public:	
    vec3D coord;
    vec3D force;
    vec3D velocity;
    int element;
    Atom(double X,double Y, double Z);
    void print_coord();
    void print_forces();
    void print_velocity();
};

class Simulation{
public:
    std::vector<Atom> mol;
    bool PBC;
    vec3D box;
    Simulation(std::vector<Atom> Mol);
    Simulation(std::vector<Atom> Mol,double x_box,double y_box, double z_box);
};


vec3D random_forces();
double pair_distance(vec3D a1,vec3D a2);
double LJ_energy_pair(double d,double eps,double sigma);
double LJ_energy_mol(std::vector<Atom> mol);
void LJ_forces(Simulation *sim);
