#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <string>
#include <cmath>

// Units are in eV, and A
// parameteres taken from https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.A405
std::map<int,float> element_epsilons = {
    {18,0.0103}
};
std::map<int,float> element_sigma = {
    {18,3.4}
};



float random_number(float r_down,float r_up){
    float r = -r_down + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(r_up-(-r_down))));
    return r;
}


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
    Atom(double X,double Y, double Z){
	coord.x = X;
	coord.y = Y;
	coord.z = Z;
	force.x = random_number(-0.1,0.1);
	force.y = random_number(-0.1,0.1);
	force.z = random_number(-0.1,0.1);
	element = 18;
    }
   
};

vec3D random_forces(){
    vec3D rf;
    float r_down = -0.2;
    float r_up = -0.2;
    rf.x = random_number(r_down,r_up);
    rf.y = random_number(r_down,r_up);
    rf.z = random_number(r_down,r_up);
    return rf;
}

void write_struc(std::vector<Atom> mol,std::string comment){
    std::ofstream outFile("struc.xyz",std::ios::app);
    outFile << mol.size() << std::endl;
    outFile << comment << std::endl;
    for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	outFile << "Ar" << " " << iter->coord.x << " " << iter->coord.y << " " << iter->coord.z << " "    << std::endl;	
    }
    outFile.close();
}

double pair_distance(vec3D a1,vec3D a2){
    double r2 = std::pow(a1.x-a2.x,2)+std::pow(a1.y-a2.y,2)+std::pow(a1.z-a2.z,2);
    double r = std::sqrt(r2);
    return r;
}

double LJ_energy_pair(double d,double eps,double sigma){
    double V_ij = 4 * eps * (std::pow((sigma / d),12) - std::pow((sigma / d),6));
    return V_ij;
}
double LJ_forces(double d,double eps,double sigma){
    double V_ij = 4 * eps * (std::pow((sigma / d),12) - std::pow((sigma / d),6));
    return V_ij;
}


double LJ_energy_mol(std::vector<Atom> mol){
    double energy = 0;
    for (auto iter = mol.begin(); iter < mol.end() ; iter++){
	auto current_el = iter->element;
	auto eps1 = element_epsilons[current_el];
	auto sigma1 = element_sigma[current_el];
	for (auto iter2 = iter+1; iter2 < mol.end() ; iter2++){
	    auto current_el = iter2->element;
	    auto eps2 = element_epsilons[current_el];
	    auto sigma2 = element_sigma[current_el];
	    auto eps_comb = std::sqrt(eps1*eps2);
	    auto sigma_comb = std::sqrt(sigma1*sigma2);
	    auto d = pair_distance(iter->coord,iter2->coord);
	    energy += LJ_energy_pair(d,eps_comb,sigma_comb);
	}
    }

    return energy;
}

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
