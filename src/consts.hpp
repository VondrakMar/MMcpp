#include <map>
#include <cmath>

// Units are in eV, A, Daltons and ps
//Boltzman constant
const double kb = 1.380649*std::pow(10,-23); // J/K
const double Avogadro = 6.0221409*std::pow(10,26);  
const double eV = 1.60218*std::pow(10,-19); // 1 eV = 1.6*10**-19 J
const double Dalton = 1.66053906660*std::pow(10,-27);
// parameteres taken from https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405



const std::map<int,float> element_epsilons = {
    // eV {18,0.0103}
    {18,1.03}
};
const std::map<int,float> element_sigma = {
    // in Angstroms
    {18,3.4}
};

const std::map<int,float> element_mass = { 
    // in Daltons
    {18,39.95}
};


const double A_to_m = 1e-10;
const double m_to_A = 1e10;