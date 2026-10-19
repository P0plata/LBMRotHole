#ifndef SOLVER_D2Q9_HPP
#define SOLVER_D2Q9_HPP

#include "grid.hpp"
#include "config_init.hpp"
#include "D_Q_.hpp"
#include "vtk_exporter.hpp"
#include <omp.h>

namespace SolverD2Q9 {
    template <CollisionType ColType>
    void run_simulation(Grid& grid, const Simconfig& cfg);
}

#endif // SOLVER_D2Q9_HPP