#pragma once
#include "lbm.hpp"
#include <string>
#include <filesystem>
#include <cstdio>

class VTKExporter {
public:
    static void export_to_vti(const Grid& grid, int step, const std::string& prefix = "data/sim");
};