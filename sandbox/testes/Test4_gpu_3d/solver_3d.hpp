// solver_3d.hpp 
#ifndef SOLVER_3D_HPP
#define SOLVER_3D_HPP

#include "grid.hpp"
#include "config_init.hpp"

#define Simconfig int // PLACEHOLDER
#define SimParams int // PLACEHOLDER

namespace Solver3D {
    template <typename Lattice>
    void run_simulation(Grid3D<Lattice>& grid, const Simconfig& cfg, const SimParams& params) noexcept;
}

#endif