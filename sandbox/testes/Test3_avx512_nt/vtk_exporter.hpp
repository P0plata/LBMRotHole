#ifndef VTK_EXPORTER_HPP
#define VTK_EXPORTER_HPP

#include "grid.hpp"
#include "D_Q_.hpp"
#include "custom_functions.hpp"
#include <string>
#include <filesystem>
#include <cstdio>


namespace VTKExporter {
    template <bool UseForces>
    void export_to_vti(const Grid& grid, const Custom::CustomForces& forces_struct, int step, const std::string& prefix = "data/sim");
};

#endif