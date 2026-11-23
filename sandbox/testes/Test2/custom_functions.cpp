#include "custom_functions.hpp"
#include "config_init.hpp"

// void Custom::apply_custom_initialization(Grid& grid, const Simconfig& cfg) {
//     const double rho0  = cfg.rho0;
//     const double u_base = (std::abs(cfg.ux0) > 1e-9) ? cfg.ux0 : 0.12; // Velocidade base (podes ajustar via ux0 no JSON)
//     const int nx       = grid.nx;
//     const int ny       = grid.ny;
//     const size_t pitch = grid.plane_pitch;
//     constexpr double PI = 3.14159265358979323846;

//     #pragma omp parallel
//     {
//         // Gerador gaussiano N(0, 1) por thread
//         std::mt19937 gen(1337 + omp_get_thread_num());
//         std::normal_distribution<double> randn(0.0, 1.0);

//         #pragma omp for schedule(static)
//         for (int y = 0; y < ny; ++y) {
//             for (int x = 0; x < nx; ++x) {
//                 const size_t i = static_cast<size_t>(y) * nx + x;

//                 // Nós sólidos começam em repouso (u = 0)
//                 if (grid.flag_field_ptr[i] == CellType::SOLID) {
//                     grid.f0_ptr[i] = D2Q9::compute_feq0(rho0, 0.0);
//                     for (int d = 0; d < 8; ++d) {
//                         const double val = D2Q9::compute_feq(d, rho0, 0.0, 0.0, 0.0);
//                         grid.f1_ptr[d * pitch + i] = val;
//                         grid.f2_ptr[d * pitch + i] = val;
//                     }
//                     continue;
//                 }

//                 // 1. Velocidade macroscópica com a onda cosseno + pequeno ruído gaussiano (1%)
//                 const double ux = u_base * (1.0 + 0.2 * std::cos(2.0 * PI * static_cast<double>(x) / nx * 4.0))
//                                 + (0.01 * u_base) * randn(gen);
//                 const double uy = (0.01 * u_base) * randn(gen);
//                 const double u_sq = ux * ux + uy * uy;

//                 // 2. Inicialização exata em equilíbrio com o teu compute_feq0 e compute_feq (Zero choque!)
//                 grid.f0_ptr[i] = D2Q9::compute_feq0(rho0, u_sq);
//                 for (int d = 0; d < 8; ++d) {
//                     const double feq = D2Q9::compute_feq(d, rho0, ux, uy, u_sq);
//                     grid.f1_ptr[d * pitch + i] = feq;
//                     grid.f2_ptr[d * pitch + i] = feq;
//                 }
//             }
//         }
//     }
// }

// void Custom::initialize_custom_forces(Grid& grid, const Simconfig& cfg, Custom::CustomForces& forces_struct) {
//     // Inicialização de quaisquer parâmetros globais necessários para as forças personalizadas
    
//     if (!cfg.use_forces) {
//         forces_struct.fx_ptr = nullptr;
//         forces_struct.fy_ptr = nullptr;
//         forces_struct.nx = grid.nx;
//         forces_struct.ny = grid.ny;
//         return;
//     }
    
//     const int center_x = cfg.center_x;
//     const int center_y = cfg.center_y;
//     const double g_mass = cfg.g_mass; // Massa do corpo central (ou outro parâmetro relevante)


//     if (!grid.fx_ptr || !grid.fy_ptr) {
//         throw std::runtime_error("Forças externas não alocadas na grelha. Certifique-se de que 'use_forces' está definido como true.");
//     }

//     forces_struct.fx_ptr = grid.fx_ptr;
//     forces_struct.fy_ptr = grid.fy_ptr;
//     forces_struct.nx = grid.nx;
//     forces_struct.ny = grid.ny;

//     #pragma omp parallel for schedule(static)
//     for (int y = 0; y < forces_struct.ny; ++y) {
//         for (int x = 0; x < forces_struct.nx; ++x) {
//             size_t idx = y * forces_struct.nx + x;
        
//             double dx = center_x - x;
//             double dy = center_y - y;
//             // Softening (eps2 = 15^2 = 225.0): evita força infinita em r=1 e cria um funil suave!
//             double r2_soft = dx * dx + dy * dy + 225.0; 
//             double r_soft  = std::sqrt(r2_soft);

//             forces_struct.fx_ptr[idx] = g_mass * dx / (r2_soft * r_soft);
//             forces_struct.fy_ptr[idx] = g_mass * dy / (r2_soft * r_soft);
//         }
        
//     }
// }

void Custom::apply_custom_initialization(Grid& grid, const Simconfig& cfg) {
    const double rho0   = cfg.rho0;
    const double u_base = (std::abs(cfg.ux0) > 1e-9) ? cfg.ux0 : 0.01;
    const int nx        = grid.nx;
    const int ny        = grid.ny;
    const size_t pitch  = grid.plane_pitch;
    constexpr double PI = 3.14159265358979323846;

    const bool pure_analytical_shear = (std::abs(cfg.uy0) < 1e-12);
    const bool black_hole_vortex     = (cfg.uy0 < -0.5); // Ativado quando uy0 = -1.0 no JSON!

    #pragma omp parallel
    {
        std::mt19937 gen(1337 + omp_get_thread_num());
        std::normal_distribution<double> randn(0.0, 1.0);

        #pragma omp for schedule(static)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const size_t i = static_cast<size_t>(y) * nx + x;

                if (grid.flag_field_ptr[i] == CellType::SOLID) {
                    grid.f0_ptr[i] = D2Q9::compute_feq0(rho0, 0.0);
                    for (int d = 0; d < 8; ++d) {
                        const double val = D2Q9::compute_feq(d, rho0, 0.0, 0.0, 0.0);
                        grid.f1_ptr[d * pitch + i] = val;
                        grid.f2_ptr[d * pitch + i] = val;
                    }
                    continue;
                }

                double ux = 0.0, uy = 0.0;

                if (black_hole_vortex && cfg.use_forces) {
                    // Inicialização em Órbita de Acreção Kepleriana + Queda Radial Suave (Espiral!)
                    const double dx = static_cast<double>(x - cfg.center_x);
                    const double dy = static_cast<double>(y - cfg.center_y);
                    const double r2 = dx * dx + dy * dy;
                    const double r  = std::sqrt(r2 + 1e-12);
                    const double r_soft = std::sqrt(r2 + 900.0); // eps = 30.0

                    // Velocidade orbital onde Força Centrífuga = Gravidade: u_theta = sqrt(g_mass * r^2 / (rho0 * r_soft^3))
                    const double u_theta = std::sqrt((cfg.g_mass / rho0) * r2 / (r_soft * r_soft * r_soft));
                    // Pequena velocidade radial de queda (infall) para formar a espiral de acreção
                    const double u_infall = -0.15 * u_theta * (30.0 / r_soft);

                    // Conversão de coordenadas polares (u_infall, u_theta) para cartesianas (ux, uy)
                    ux = (u_infall * dx - u_theta * dy) / r;
                    uy = (u_infall * dy + u_theta * dx) / r;
                } else if (pure_analytical_shear) {
                    ux = u_base * std::sin(2.0 * PI * static_cast<double>(y) / static_cast<double>(ny));
                    uy = 0.0;
                } else {
                    ux = u_base * (1.0 + 0.2 * std::cos(2.0 * PI * static_cast<double>(x) / nx * 4.0))
                       + (0.01 * u_base) * randn(gen);
                    uy = (0.01 * u_base) * randn(gen);
                }

                const double u_sq = ux * ux + uy * uy;
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

void Custom::initialize_custom_forces(Grid& grid, const Simconfig& cfg, Custom::CustomForces& forces_struct) {
    if (!cfg.use_forces) {
        forces_struct.fx_ptr = nullptr;
        forces_struct.fy_ptr = nullptr;
        forces_struct.nx = grid.nx;
        forces_struct.ny = grid.ny;
        return;
    }

    const int center_x  = cfg.center_x;
    const int center_y  = cfg.center_y;
    const double g_mass = cfg.g_mass;
    const bool black_hole_vortex = (cfg.uy0 < -0.5);

    if (!grid.fx_ptr || !grid.fy_ptr) {
        throw std::runtime_error("Forças externas não alocadas na grelha.");
    }

    forces_struct.fx_ptr = grid.fx_ptr;
    forces_struct.fy_ptr = grid.fy_ptr;
    forces_struct.nx = grid.nx;
    forces_struct.ny = grid.ny;

    // Se estiver no modo Vórtice de Buraco Negro, usa eps^2 = 900 (eps=30); nos testes standard mantém 225 (eps=15)
    const double eps2 = black_hole_vortex ? 900.0 : 225.0;

    #pragma omp parallel for schedule(static)
    for (int y = 0; y < forces_struct.ny; ++y) {
        for (int x = 0; x < forces_struct.nx; ++x) {
            size_t idx = static_cast<size_t>(y) * forces_struct.nx + x;

            if (center_x < 0) {
                forces_struct.fx_ptr[idx] = g_mass;
                forces_struct.fy_ptr[idx] = 0.0;
            }else {
                double total_fx = 0.0;
                double total_fy = 0.0;

                // Se houver buracos (SPONGE) na lista de obstacles, soma a força de cada um deles!
                if (!cfg.obstacles.empty() && std::abs(cfg.rot_factor) > 1e-9) {
                    for (const auto& obs : cfg.obstacles) {
                        double dx = obs.offset_x - x;
                        double dy = obs.offset_y - y;
                        double r2 = dx * dx + dy * dy;

                        double r2_soft = r2 + 3600.0;
                        double r_soft  = std::sqrt(r2_soft);
                        double inv_r3  = 1.0 / (r2_soft * r_soft);

                        double f_rad = g_mass * inv_r3;
                        double f_rot = g_mass * cfg.rot_factor * inv_r3 / r_soft;

                        total_fx += f_rad * dx + f_rot * dy;
                        total_fy += f_rad * dy - f_rot * dx;
                    }
                } else {
                    // Caso com 1 único centro (center_x, center_y) ou testes de validação
                    double dx = center_x - x;
                    double dy = center_y - y;
                    double r2 = dx * dx + dy * dy;

                    if (std::abs(cfg.rot_factor) > 1e-9) {
                        double r2_soft = r2 + 3600.0;
                        double r_soft  = std::sqrt(r2_soft);
                        double inv_r3  = 1.0 / (r2_soft * r_soft);

                        double f_rad = g_mass * inv_r3;
                        double f_rot = g_mass * cfg.rot_factor * inv_r3 / r_soft;

                        total_fx = f_rad * dx + f_rot * dy;
                        total_fy = f_rad * dy - f_rot * dx;
                    } else {
                        double r2_soft = r2 + 225.0;
                        double r_soft  = std::sqrt(r2_soft);
                        total_fx = g_mass * dx / (r2_soft * r_soft);
                        total_fy = g_mass * dy / (r2_soft * r_soft);
                    }
                }

                forces_struct.fx_ptr[idx] = total_fx;
                forces_struct.fy_ptr[idx] = total_fy;
            }
        }
    }
}