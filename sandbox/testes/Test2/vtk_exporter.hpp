#ifndef VTK_EXPORTER_HPP
#define VTK_EXPORTER_HPP

#include "grid.hpp"
#include "D_Q_.hpp"
#include <string>
#include <filesystem>
#include <cstdio>


namespace VTKExporter {
    void export_to_vti(const Grid& grid, int step, const std::string& prefix = "data/sim");
};

#endif