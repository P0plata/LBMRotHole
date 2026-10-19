
#ifndef UTILS_HPP
#define UTILS_HPP

#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>


namespace utils {

    size_t get_system_available_ram();
    void check_memory_availability(size_t max_allowed_bytes);



}

#endif