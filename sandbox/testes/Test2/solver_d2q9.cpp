#include "solver_d2q9.hpp"

namespace SolverD2Q9 {

// ============================================================================
// 1. AS 3 FUNÇÕES DE COLISÃO SEPARADAS (static inline -> Zero overhead)
// ============================================================================

static inline void collide_bgk(double& f0, double f[8], double rho, double ux, double uy, double omega) noexcept {
    const double u_sq = ux * ux + uy * uy;

    f0   -= omega * (f0 - D2Q9::compute_feq0(rho, u_sq));
    f[0] -= omega * (f[0] - D2Q9::compute_feq(0, rho, ux, uy, u_sq));
    f[1] -= omega * (f[1] - D2Q9::compute_feq(1, rho, ux, uy, u_sq));
    f[2] -= omega * (f[2] - D2Q9::compute_feq(2, rho, ux, uy, u_sq));
    f[3] -= omega * (f[3] - D2Q9::compute_feq(3, rho, ux, uy, u_sq));
    f[4] -= omega * (f[4] - D2Q9::compute_feq(4, rho, ux, uy, u_sq));
    f[5] -= omega * (f[5] - D2Q9::compute_feq(5, rho, ux, uy, u_sq));
    f[6] -= omega * (f[6] - D2Q9::compute_feq(6, rho, ux, uy, u_sq));    
    f[7] -= omega * (f[7] - D2Q9::compute_feq(7, rho, ux, uy, u_sq));
}

static inline void collide_trt(double& f0, double f[8], double rho, double ux, double uy,
                               double omega_p, double omega_m) noexcept {
    const double u_sq = ux * ux + uy * uy;
    f0 -= omega_p * (f0 - D2Q9::compute_feq0(rho, u_sq));

    constexpr int pairs[4][2] = {{0, 2}, {1, 3}, {4, 6}, {5, 7}};
    for (int p = 0; p < 4; ++p) {
        const int d     = pairs[p][0];
        const int d_opp = pairs[p][1];

        const double feq_d   = D2Q9::compute_feq(d, rho, ux, uy, u_sq);
        const double feq_opp = D2Q9::compute_feq(d_opp, rho, ux, uy, u_sq);

        const double f_plus    = 0.5 * (f[d] + f[d_opp]);
        const double f_minus   = 0.5 * (f[d] - f[d_opp]);
        const double feq_plus  = 0.5 * (feq_d + feq_opp);
        const double feq_minus = 0.5 * (feq_d - feq_opp);

        const double corr_p = omega_p * (f_plus - feq_plus);
        const double corr_m = omega_m * (f_minus - feq_minus);

        f[d]     -= (corr_p + corr_m);
        f[d_opp] -= (corr_p - corr_m);
    }
}

static inline void collide_mrt(double& f0, double f[8], double rho, double ux, double uy,
                               const double s_mrt[9]) noexcept {
    const double u_sq = ux * ux + uy * uy;
    const double f_card = f[0] + f[1] + f[2] + f[3];
    const double f_diag = f[4] + f[5] + f[6] + f[7];

    // Momentos atuais (m1=e, m2=eps, m4=qx, m6=qy, m7=pxx, m8=pxy)
    const double m1 = -4.0 * f0 - f_card + 2.0 * f_diag;
    const double m2 =  4.0 * f0 - 2.0 * f_card + f_diag;
    const double m4 = -2.0 * (f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7]);
    const double m6 = -2.0 * (f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7]);
    const double m7 = (f[0] - f[1] + f[2] - f[3]);
    const double m8 = (f[4] - f[5] + f[6] - f[7]);

    // Momentos de equilíbrio
    const double m1_eq = rho * (-2.0 + 3.0 * u_sq);
    const double m2_eq = rho * ( 1.0 - 3.0 * u_sq);
    const double m4_eq = -rho * ux;
    const double m6_eq = -rho * uy;
    const double m7_eq = rho * (ux * ux - uy * uy);
    const double m8_eq = rho * (ux * uy);

    const double dm1 = s_mrt[1] * (m1 - m1_eq) * (1.0 / 36.0);
    const double dm2 = s_mrt[2] * (m2 - m2_eq) * (1.0 / 36.0);
    const double dm4 = s_mrt[4] * (m4 - m4_eq) * (1.0 / 12.0);
    const double dm6 = s_mrt[6] * (m6 - m6_eq) * (1.0 / 12.0);
    const double dm7 = s_mrt[7] * (m7 - m7_eq) * (1.0 / 4.0);
    const double dm8 = s_mrt[8] * (m8 - m8_eq) * (1.0 / 4.0);

    f0   -= ( -4.0 * dm1 + 4.0 * dm2 );
    f[0] -= ( -dm1 - 2.0 * dm2 - 2.0 * dm4 + dm7 );
    f[1] -= ( -dm1 - 2.0 * dm2 - 2.0 * dm6 - dm7 );
    f[2] -= ( -dm1 - 2.0 * dm2 + 2.0 * dm4 + dm7 );
    f[3] -= ( -dm1 - 2.0 * dm2 + 2.0 * dm6 - dm7 );
    f[4] -= (  2.0 * dm1 + dm2 + dm4 + dm6 + dm8 );
    f[5] -= (  2.0 * dm1 + dm2 - dm4 + dm6 - dm8 );
    f[6] -= (  2.0 * dm1 + dm2 - dm4 - dm6 + dm8 );
    f[7] -= (  2.0 * dm1 + dm2 + dm4 - dm6 - dm8 );
}

// Disparo em tempo de compilação para a função de colisão correta
template <CollisionType ColType>
static inline void dispatch_collision(double& f0, double f[8], double rho, double ux, double uy,
                                      double omega_bgk, double omega_trt, const double s_mrt[9]) noexcept {
    if constexpr (ColType == CollisionType::BGK) {
        collide_bgk(f0, f, rho, ux, uy, omega_bgk);
    } else if constexpr (ColType == CollisionType::TRT) {
        collide_trt(f0, f, rho, ux, uy, omega_bgk, omega_trt);
    } else if constexpr (ColType == CollisionType::MRT) {
        collide_mrt(f0, f, rho, ux, uy, s_mrt);
    }
}

// ============================================================================
// 2. PULL + COLLIDE POR NÓ DE FLUIDO
// ============================================================================
template <CollisionType ColType, bool IsInterior>
static inline void pull_and_collide_node(int x, int y, int nx, int ny, size_t pitch, Grid& grid,
                                         double omega_bgk, double omega_trt, const double s_mrt[9]) noexcept {
    const size_t idx = static_cast<size_t>(y) * nx + x;
    const CellType cell = grid.flag_field_ptr[idx];
    if (cell != CellType::FLUID && cell != CellType::SPONGE) return;

    double feq;
    double f0 = grid.f0_ptr[idx];
    double f[8];

    /*for (int d = 0; d < 8; ++d) {
        size_t src_idx;
        if constexpr (IsInterior) {
            src_idx = idx - (D2Q9::cx[d] + D2Q9::cy[d] * nx);
        } else {
            int src_x = (x - D2Q9::cx[d] + nx) % nx;
            int src_y = (y - D2Q9::cy[d] + ny) % ny;
            src_idx = static_cast<size_t>(src_y) * nx + src_x;
        }
        f[d] = grid.f1_ptr[d * pitch + src_idx];
    }*/

    if constexpr (IsInterior) {
        
        f[0] = grid.f1_ptr[0 * pitch + D2Q9::src_idx(idx,nx,0)];
        f[1] = grid.f1_ptr[1 * pitch + D2Q9::src_idx(idx,nx,1)];
        f[2] = grid.f1_ptr[2 * pitch + D2Q9::src_idx(idx,nx,2)];
        f[3] = grid.f1_ptr[3 * pitch + D2Q9::src_idx(idx,nx,3)];
        f[4] = grid.f1_ptr[4 * pitch + D2Q9::src_idx(idx,nx,4)];
        f[5] = grid.f1_ptr[5 * pitch + D2Q9::src_idx(idx,nx,5)];
        f[6] = grid.f1_ptr[6 * pitch + D2Q9::src_idx(idx,nx,6)];
        f[7] = grid.f1_ptr[7 * pitch + D2Q9::src_idx(idx,nx,7)];            

    } else {

        f[0] = grid.f1_ptr[0 * pitch + D2Q9::src_idx_p(x,y,nx,ny,0)];
        f[1] = grid.f1_ptr[1 * pitch + D2Q9::src_idx_p(x,y,nx,ny,1)];
        f[2] = grid.f1_ptr[2 * pitch + D2Q9::src_idx_p(x,y,nx,ny,2)];
        f[3] = grid.f1_ptr[3 * pitch + D2Q9::src_idx_p(x,y,nx,ny,3)];
        f[4] = grid.f1_ptr[4 * pitch + D2Q9::src_idx_p(x,y,nx,ny,4)];
        f[5] = grid.f1_ptr[5 * pitch + D2Q9::src_idx_p(x,y,nx,ny,5)];
        f[6] = grid.f1_ptr[6 * pitch + D2Q9::src_idx_p(x,y,nx,ny,6)];
        f[7] = grid.f1_ptr[7 * pitch + D2Q9::src_idx_p(x,y,nx,ny,7)];    
    
    }

    const double rho = f0 + f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
    const double inv_rho = 1.0 / rho;
    const double ux = ((f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7])) * inv_rho;
    const double uy = ((f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7])) * inv_rho;

    dispatch_collision<ColType>(f0, f, rho, ux, uy, omega_bgk, omega_trt, s_mrt);

    grid.f0_ptr[idx] = f0;
    grid.f2_ptr[0 * pitch + idx] = f[0];
    grid.f2_ptr[1 * pitch + idx] = f[1];
    grid.f2_ptr[2 * pitch + idx] = f[2];
    grid.f2_ptr[3 * pitch + idx] = f[3];
    grid.f2_ptr[4 * pitch + idx] = f[4];
    grid.f2_ptr[5 * pitch + idx] = f[5];
    grid.f2_ptr[6 * pitch + idx] = f[6];
    grid.f2_ptr[7 * pitch + idx] = f[7];
}

// ============================================================================
// 3. FRONTEIRAS (Esponja, Zou-He Inlet/Outlet e Bouzidi/Bounce-Back)
// ============================================================================
static inline void apply_sponges(Grid& grid) noexcept {
    const size_t pitch = grid.plane_pitch;
    const size_t n_sponges = grid.sponge_nodes.size();

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n_sponges; ++i) {
        const SpongeNode& sp = grid.sponge_nodes[i];
        const size_t idx = sp.idx;
        const double sigma = sp.sigma;
        const double oms = 1.0 - sigma;
        const double u_sq = static_cast<double>(sp.vx_ref) * sp.vx_ref + static_cast<double>(sp.vy_ref) * sp.vy_ref;

        grid.f0_ptr[idx] = oms * grid.f0_ptr[idx] + sigma * D2Q9::compute_feq0(sp.rho_ref, u_sq);
        for (int d = 0; d < 8; ++d) {
            double feq_ref = D2Q9::compute_feq(d, sp.rho_ref, sp.vx_ref, sp.vy_ref, u_sq);
            grid.f2_ptr[d * pitch + idx] = oms * grid.f2_ptr[d * pitch + idx] + sigma * feq_ref;
        }
    }
}

template <CollisionType ColType, bool IsInlet>
static inline void apply_open_boundaries(Grid& grid, const std::vector<OpenBoundaryNode>& nodes,
                                         double omega_bgk, double omega_trt, const double s_mrt[9]) noexcept {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const size_t pitch = grid.plane_pitch;
    const size_t n_nodes = nodes.size();

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n_nodes; ++i) {
        const OpenBoundaryNode& bnd = nodes[i];
        const size_t idx = bnd.idx;
        const int x = static_cast<int>(idx % nx);
        const int y = static_cast<int>(idx / nx);

        double f0 = grid.f0_ptr[idx];
        double f[8];

        // 1. Puxar populações conhecidas dos vizinhos interiores
        for (int d = 0; d < 8; ++d) {
            int src_x = x - D2Q9::cx[d];
            int src_y = y - D2Q9::cy[d];
            if (src_x >= 0 && src_x < nx && src_y >= 0 && src_y < ny) {
                f[d] = grid.f1_ptr[d * pitch + (static_cast<size_t>(src_y) * nx + src_x)];
            } else {
                f[d] = 0.0;
            }
        }

        double rho = 1.0, ux = 0.0, uy = 0.0;

        // 2. Reconstrução de Zou-He consoante a normal (0:+x, 1:-x, 2:+y, 3:-y)
        if constexpr (IsInlet) {
            ux = bnd.vx;
            uy = bnd.vy;
            if (bnd.normal == 0) { // Face Oeste (x=0, normal +x): desconhecidos 0, 4, 7
                rho = (f0 + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / (1.0 - ux);
                f[0] = f[2] + (2.0 / 3.0) * rho * ux;
                f[4] = f[6] - 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * ux + 0.5 * rho * uy;
                f[7] = f[5] + 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * ux - 0.5 * rho * uy;
            } else if (bnd.normal == 1) { // Face Este (x=nx-1, normal -x): desconhecidos 2, 5, 6
                rho = (f0 + f[1] + f[3] + 2.0 * (f[0] + f[4] + f[7])) / (1.0 + ux);
                f[2] = f[0] - (2.0 / 3.0) * rho * ux;
                f[5] = f[7] - 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * ux + 0.5 * rho * uy;
                f[6] = f[4] + 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * ux - 0.5 * rho * uy;
            } else if (bnd.normal == 2) { // Face Sul (y=0, normal +y): desconhecidos 1, 4, 5
                rho = (f0 + f[0] + f[2] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - uy);
                f[1] = f[3] + (2.0 / 3.0) * rho * uy;
                f[4] = f[6] - 0.5 * (f[0] - f[2]) + 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
                f[5] = f[7] + 0.5 * (f[0] - f[2]) - 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
            } else if (bnd.normal == 3) { // Face Norte (y=ny-1, normal -y): desconhecidos 3, 6, 7
                rho = (f0 + f[0] + f[2] + 2.0 * (f[1] + f[4] + f[5])) / (1.0 + uy);
                f[3] = f[1] - (2.0 / 3.0) * rho * uy;
                f[6] = f[4] + 0.5 * (f[0] - f[2]) - 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
                f[7] = f[5] - 0.5 * (f[0] - f[2]) + 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
            }
        } else {
            // Outlet de pressão/densidade imposta (rho guardado em bnd.vx)
            rho = bnd.vx;
            if (bnd.normal == 0) {
                ux = -1.0 + (f0 + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / rho;
                f[0] = f[2] + (2.0 / 3.0) * rho * ux;
                f[4] = f[6] - 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * ux;
                f[7] = f[5] + 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * ux;
            } else if (bnd.normal == 1) {
                ux = 1.0 - (f0 + f[1] + f[3] + 2.0 * (f[0] + f[4] + f[7])) / rho;
                f[2] = f[0] - (2.0 / 3.0) * rho * ux;
                f[5] = f[7] - 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * ux;
                f[6] = f[4] + 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * ux;
            } else if (bnd.normal == 2) {
                uy = -1.0 + (f0 + f[0] + f[2] + 2.0 * (f[3] + f[6] + f[7])) / rho;
                f[1] = f[3] + (2.0 / 3.0) * rho * uy;
                f[4] = f[6] - 0.5 * (f[0] - f[2]) + (1.0 / 6.0) * rho * uy;
                f[5] = f[7] + 0.5 * (f[0] - f[2]) + (1.0 / 6.0) * rho * uy;
            } else if (bnd.normal == 3) {
                uy = 1.0 - (f0 + f[0] + f[2] + 2.0 * (f[1] + f[4] + f[5])) / rho;
                f[3] = f[1] - (2.0 / 3.0) * rho * uy;
                f[6] = f[4] + 0.5 * (f[0] - f[2]) - (1.0 / 6.0) * rho * uy;
                f[7] = f[5] - 0.5 * (f[0] - f[2]) - (1.0 / 6.0) * rho * uy;
            }
        }

        // 3. Colidir também o nó de fronteira e gravar em f2_ptr
        dispatch_collision<ColType>(f0, f, rho, ux, uy, omega_bgk, omega_trt, s_mrt);

        grid.f0_ptr[idx] = f0;
        for (int d = 0; d < 8; ++d) {
            grid.f2_ptr[d * pitch + idx] = f[d];
        }
    }
}

static inline void apply_wall_links(Grid& grid, double* f_target, bool periodic_x, bool periodic_y) noexcept {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const size_t pitch = grid.plane_pitch;
    const size_t n_links = grid.wall_links.size();

    // Sem nowait aqui: funciona como barreira final antes do swap_populations()!
    #pragma omp for schedule(static)
    for (size_t i = 0; i < n_links; ++i) {
        const WallLink& wl = grid.wall_links[i];
        const size_t xf  = wl.fluid_idx;
        const size_t xff = wl.neighbor_f;
        const uint8_t d     = wl.dir;
        const uint8_t d_opp = wl.opp_dir;
        const double q      = wl.q_frac;

        int x = static_cast<int>(xf % nx) + D2Q9::cx[d];
        int y = static_cast<int>(xf / nx) + D2Q9::cy[d];
        if (periodic_x) x = (x + nx) % nx;
        if (periodic_y) y = (y + ny) % ny;
        const size_t solid_idx = static_cast<size_t>(y) * nx + x;

        const double f_d_xf = f_target[d * pitch + xf];
        double f_reflected;

        if (q == 0.5) {
            f_reflected = f_d_xf;
        } else if (q < 0.5) {
            f_reflected = 2.0 * q * f_d_xf + (1.0 - 2.0 * q) * f_target[d * pitch + xff];
        } else {
            const double inv_2q = 1.0 / (2.0 * q);
            f_reflected = inv_2q * f_d_xf + (1.0 - inv_2q) * f_target[d_opp * pitch + xf];
        }

        f_target[d_opp * pitch + solid_idx] = f_reflected;
    }
}

// ============================================================================
// 4. RUN_SIMULATION (A Função Principal que Funde Tudo)
// ============================================================================
template <CollisionType ColType>
void run_simulation(Grid& grid, const Simconfig& cfg) {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const size_t pitch = grid.plane_pitch;

    const double omega_bgk = (cfg.tau_BGK > 0.0) ? (1.0 / cfg.tau_BGK) : 1.0;
    const double omega_trt = (cfg.tau_TRT > 0.0) ? (1.0 / cfg.tau_TRT) : 1.0;
    const bool per_x = (cfg.boundary_conditions[0] == BoundaryCondition::PERIODIC);
    const bool per_y = (cfg.boundary_conditions[2] == BoundaryCondition::PERIODIC);

    // As threads OpenMP são criadas aqui UMA ÚNICA VEZ
    #pragma omp parallel
    {
        // Inicializar os slots SOLID em f1_ptr antes do passo 0
        apply_wall_links(grid, grid.f1_ptr, per_x, per_y);

        for (int step = 0; step <= cfg.max_steps; ++step) {

            // 1. Exportar VTK na frequência pedida
            if (cfg.export_freq > 0 && (step % cfg.export_freq == 0)) {
                #pragma omp single
                {
                    VTKExporter::export_to_vti(grid, step, cfg.output_prefix);
                }
            }

            if (step == cfg.max_steps) break;

            // 2. Fused Pull + Collide em todo o domínio fluido
                #pragma omp for schedule(static) nowait
                for (int y = 1; y < ny - 1; ++y) {
                    for (int x = 1; x < nx - 1; ++x) {
                        pull_and_collide_node<ColType, true>(x, y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT);
                    }
                }

                // 2. BORDAS PERIÓDICAS (Só corre se a fronteira for realmente PERIODIC)
                if (per_y) {
                    #pragma omp for schedule(static) nowait
                    for (int x = 0; x < nx; ++x) {
                        pull_and_collide_node<ColType, false>(x, 0,      nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT);
                        pull_and_collide_node<ColType, false>(x, ny - 1, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT);
                    }
                }
                if (per_x) {
                    #pragma omp for schedule(static) nowait
                    for (int y = 1; y < ny - 1; ++y) {
                        pull_and_collide_node<ColType, false>(0,      y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT);
                        pull_and_collide_node<ColType, false>(nx - 1, y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT);
                    }
                }
                #pragma omp barrier // Sincroniza tudo antes de aplicar as listas de fronteira

            // 3. Aplicar Fronteiras (Esponjas, Inlets, Outlets correm com 'nowait' e sincronizam no fim de wall_links)
            apply_sponges(grid);
            apply_open_boundaries<ColType, true>(grid, grid.inlet_nodes, omega_bgk, omega_trt, cfg.s_MRT);
            apply_open_boundaries<ColType, false>(grid, grid.outlet_nodes, omega_bgk, omega_trt, cfg.s_MRT);
            apply_wall_links(grid, grid.f2_ptr, per_x, per_y); // Barreira implícita aqui!

            // 4. Trocar ponteiros f1_ptr <-> f2_ptr
            #pragma omp single
            {
                grid.swap_populations();
            }
        }
    }
}

// ============================================================================
// 5. INSTANCIAÇÃO EXPLÍCITA DOS TEMPLATES (Gera as 3 versões no .o)
// ============================================================================
template void run_simulation<CollisionType::BGK>(Grid& grid, const Simconfig& cfg);
template void run_simulation<CollisionType::TRT>(Grid& grid, const Simconfig& cfg);
template void run_simulation<CollisionType::MRT>(Grid& grid, const Simconfig& cfg);

} // namespace SolverD2Q9


