
#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>    // <-- Necessário para std::setprecision e std::scientific
#include <cmath>      // <-- Necessário para std::isfinite
#include <set>
#include <filesystem>
#include "grid.hpp"



namespace utils {

    size_t get_system_available_ram();
    void check_memory_availability(size_t max_allowed_bytes);
    size_t get_system_l3_cache_bytes();
    void audit_grid_integrity(const Grid& grid);

}

#endif