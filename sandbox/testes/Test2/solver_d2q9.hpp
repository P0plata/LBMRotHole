#ifndef SOLVER_D2Q9_HPP
#define SOLVER_D2Q9_HPP

#include "grid.hpp"
#include "custom_functions.hpp"
#include "config_init.hpp"
#include "D_Q_.hpp"
#include "vtk_exporter.hpp"
#include <omp.h>

namespace SolverD2Q9 {
    template <CollisionType ColType, bool UseForces>
    void run_simulation(Grid& grid, const Simconfig& cfg, Custom::CustomForces& forces_struct) noexcept;
}

#endif // SOLVER_D2Q9_HPP