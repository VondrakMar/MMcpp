#include <vector>
#include <string>
#include <map>

const std::map<int,float> element_epsilons = {
    {18,0.0103}
};
const std::map<int,float> element_sigma = {
    {18,3.4}
};

float random_number(float r_down,float r_up);

struct vec3D{
    double x;
    double y;
    double z;
};

class Atom{
public:	
    vec3D coord;
    vec3D force;
    int element;
    Atom(double X,double Y, double Z);
};

vec3D random_forces();

void write_struc(std::vector<Atom> mol,std::string comment);
double pair_distance(vec3D a1,vec3D a2);
double LJ_energy_pair(double d,double eps,double sigma);
double LJ_forces(double d,double eps,double sigma);
double LJ_energy_mol(std::vector<Atom> mol);

