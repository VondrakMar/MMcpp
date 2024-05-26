#include <vector>
#include <string>
#include <cmath>

#ifndef STRUC_H
#define STRUC_H

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
    double E=0;
    Simulation(int nAt);
    Simulation(std::vector<Atom> Mol);
    Simulation(int nAt,double x_box,double y_box, double z_box);
    Simulation(std::vector<Atom> Mol,double x_box,double y_box, double z_box);
    void rescale_v(double T);
};


vec3D random_forces();
double pair_distance(vec3D a1,vec3D a2);
double LJ_energy_pair(double d,double eps,double sigma);
double LJ_energy_mol(std::vector<Atom> mol);
void LJ_forces(Simulation *sim);

#endif