#ifndef INITIALIZER_HPP
#define INITIALIZER_HPP

#include "grid.hpp"
#include <string>

struct Simconfig {
    
    // Config Parameters
    
    // Simulation 
    int nx{0};
    int ny{0};
    int nz{0};
    int grid_size{0}; // Computed as nx * ny * nz
    int grid_type{0}; // 0: 2D, 1: 3D
    int max_steps{0};
    int export_freq{};
    std::string output_prefix{""};
    
    // Geometry
    std::string geometry_type{""}; // circle, rectangle, tube, file
    std::string filename{""};

    // Collision Model
    int type_collision{0}; // 0: BGK, 1: TRT, 2: MRT
    double tau_BGK{-1}; // Relaxation time BGK
    double tau_TRT{-1}; // Relaxation time TRT
    std::string tau_MRT_filename{""}; // Filename for MRT relaxation times

    // Forcing
    // Tenho de estudar mais

    // Boundary Conditions
    //planes x0,x1,y0,y1,z0,z1
    int boundary_conditions[6]; //  0: periodic, 1: Dirichlet, 2: Neumann 
    int type_boundary_solver[6]; // 0: bounce-back, 1: Zou-he

    


};

class Initializer{

    public:
        static void initialize(Grid& grid, const Simconfig& config);
};

#endif // INITIALIZER_HPP