
#ifndef CUSTOM_FUNCTIONS_HPP
#define CUSTOM_FUNCTIONS_HPP

struct Simconfig;

#include "D_Q_.hpp"
#include "grid.hpp"
#include <random>
#include <cmath>
#include <omp.h>
#include <string>


namespace Custom {

    struct CustomForces {
        double* fx_ptr{nullptr}; // Forças externas X (opcional)
        double* fy_ptr{nullptr}; // Forças externas Y (opcional)
        int nx{0};
        int ny{0};
    };

    void initialize_custom_forces(Grid& grid, const Simconfig& cfg, CustomForces& forces_struct);

    void apply_custom_initialization(Grid& grid, const Simconfig& cfg);

    double inline compute_custom_fx(int x, int y, int step, CustomForces& forces_struct) {
    // Exemplo: Força oscilante no eixo X
    // Ajusta a amplitude e frequência conforme necessário
    return forces_struct.fx_ptr[y * forces_struct.nx + x]; // Retorna a força pré-calculada para a célula (x, y)
    
    }

    double inline compute_custom_fy(int x, int y, int step, CustomForces& forces_struct) {
        // Exemplo: Força oscilante no eixo Y
        return forces_struct.fy_ptr[y * forces_struct.nx + x]; // Retorna a força pré-calculada para a célula (x, y)
    }

}

#endif // CUSTOM_FUNCTIONS_HPP