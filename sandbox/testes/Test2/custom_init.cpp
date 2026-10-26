#include "config_init.hpp"
#include "D_Q_.hpp"
#include <random>
#include <cmath>
#include <omp.h>

void Initializer::apply_custom_initialization(Grid& grid, const Simconfig& cfg) {
    const double rho0  = cfg.rho0;
    const double u_base = (std::abs(cfg.ux0) > 1e-9) ? cfg.ux0 : 0.12; // Velocidade base (podes ajustar via ux0 no JSON)
    const int nx       = grid.nx;
    const int ny       = grid.ny;
    const size_t pitch = grid.plane_pitch;
    constexpr double PI = 3.14159265358979323846;

    #pragma omp parallel
    {
        // Gerador gaussiano N(0, 1) por thread
        std::mt19937 gen(1337 + omp_get_thread_num());
        std::normal_distribution<double> randn(0.0, 1.0);

        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const size_t i = static_cast<size_t>(y) * nx + x;

                // Nós sólidos começam em repouso (u = 0)
                if (grid.flag_field_ptr[i] == CellType::SOLID) {
                    grid.f0_ptr[i] = D2Q9::compute_feq0(rho0, 0.0);
                    for (int d = 0; d < 8; ++d) {
                        const double val = D2Q9::compute_feq(d, rho0, 0.0, 0.0, 0.0);
                        grid.f1_ptr[d * pitch + i] = val;
                        grid.f2_ptr[d * pitch + i] = val;
                    }
                    continue;
                }

                // 1. Velocidade macroscópica com a onda cosseno + pequeno ruído gaussiano (1%)
                const double ux = u_base * (1.0 + 0.2 * std::cos(2.0 * PI * static_cast<double>(x) / nx * 4.0))
                                + (0.01 * u_base) * randn(gen);
                const double uy = (0.01 * u_base) * randn(gen);
                const double u_sq = ux * ux + uy * uy;

                // 2. Inicialização exata em equilíbrio com o teu compute_feq0 e compute_feq (Zero choque!)
                grid.f0_ptr[i] = D2Q9::compute_feq0(rho0, u_sq);
                for (int d = 0; d < 8; ++d) {
                    const double feq = D2Q9::compute_feq(d, rho0, ux, uy, u_sq);
                    grid.f1_ptr[d * pitch + i] = feq;
                    grid.f2_ptr[d * pitch + i] = feq;
                }
            }
        }
    }
}