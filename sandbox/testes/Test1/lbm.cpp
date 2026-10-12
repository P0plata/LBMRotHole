#include "lbm.h"
#include <cstdio>
#include <filesystem>

namespace fs = std::filesystem;

void initialize(Grid *g) {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int idx = get_idx(x, y);
            for (int q = 0; q < 9; q++) {
                g->f[q][idx] = w[q] * 1.0; // rho inicial = 1.0
            }
        }
    }
}

// Calcula densidade e velocidades reais do fluido
void compute_macroscopic(const Grid *g, double rho[NUM_CELLS], double ux[NUM_CELLS], double uy[NUM_CELLS]) {
    for (int idx = 0; idx < NUM_CELLS; idx++) {
        double r = 0.0;
        double mx = 0.0;
        double my = 0.0;

        for (int q = 0; q < 9; q++) {
            double val = g->f[q][idx];
            r += val;
            mx += val * cx[q];
            my += val * cy[q];
        }

        rho[idx] = r;
        // Com forçamento de corpo, a velocidade física real inclui metade da força
        ux[idx] = (mx + 0.5 * FORCE_X) / r;
        uy[idx] = my / r;
    }
}

// Equilíbrio D2Q9 puramente algébrico
void compute_equilibrium(const double rho[NUM_CELLS], const double ux[NUM_CELLS], const double uy[NUM_CELLS], Grid *f_eq) {
    for (int idx = 0; idx < NUM_CELLS; idx++) {
        double r = rho[idx];
        double u_x = ux[idx];
        double u_y = uy[idx];
        double u_sq = u_x * u_x + u_y * u_y;

        for (int q = 0; q < 9; q++) {
            double c_dot_u = cx[q] * u_x + cy[q] * u_y;
            f_eq->f[q][idx] = w[q] * r * (1.0 + 3.0 * c_dot_u + 4.5 * (c_dot_u * c_dot_u) - 1.5 * u_sq);
        }
    }
}

// Operador BGK padrão (se quiseres MRT ou TRT, crias uma função collide_mrt aqui)
void collide_bgk(const Grid *f_in, const Grid *f_eq, Grid *f_out) {
    for (int q = 0; q < 9; q++) {
        for (int idx = 0; idx < NUM_CELLS; idx++) {
            f_out->f[q][idx] = f_in->f[q][idx] - OMEGA * (f_in->f[q][idx] - f_eq->f[q][idx]);
        }
    }
}

// Adiciona o efeito da força externa às populações pós-colisão
void apply_forcing(const double ux[NUM_CELLS], const double uy[NUM_CELLS], Grid *f_post) {
    for (int idx = 0; idx < NUM_CELLS; idx++) {
        double u_x = ux[idx];
        double u_y = uy[idx];

        for (int q = 0; q < 9; q++) {
            double c_dot_u = cx[q] * u_x + cy[q] * u_y;
            double term = (1.0 - 0.5 * OMEGA) * w[q] * (3.0 * (cx[q] - u_x) + 9.0 * c_dot_u * cx[q]) * FORCE_X;
            f_post->f[q][idx] += term;
        }
    }
}

// Propagação pura no interior do domínio com periodicidade em X
void stream(const Grid *f_src, Grid *f_dst) {
    for (int y = 1; y < NY - 1; y++) {
        for (int x = 0; x < NX; x++) {
            int src_idx = get_idx(x, y);

            for (int q = 0; q < 9; q++) {
                int next_x = (x + cx[q] + NX) % NX;
                int next_y = y + cy[q];

                // Só propaga se o destino estiver dentro do fluido ativo
                if (next_y > 0 && next_y < NY - 1) {
                    int dst_idx = get_idx(next_x, next_y);
                    f_dst->f[q][dst_idx] = f_src->f[q][src_idx];
                }
            }
        }
    }
}

// Aplicação de bounce-back nas paredes (y = 0 e y = NY - 1)
void apply_boundary_conditions(const Grid *f_post_collision, Grid *f_dst) {
    for (int y = 1; y < NY - 1; y++) {
        for (int x = 0; x < NX; x++) {
            int idx = get_idx(x, y);

            for (int q = 0; q < 9; q++) {
                int target_y = y + cy[q];

                // Se a direção apontar para fora da parede, reflete na direção oposta
                if (target_y == 0 || target_y == NY - 1) {
                    f_dst->f[opp[q]][idx] = f_post_collision->f[q][idx];
                }
            }
        }
    }
}

void export_to_vtk(const double rho[NUM_CELLS], const double ux[NUM_CELLS], const double uy[NUM_CELLS], int step) {
    fs::create_directories("data");
    char filename[128];
    std::snprintf(filename, sizeof(filename), "data/poiseuille_%06d.vtk", step);

    FILE *fp = std::fopen(filename, "w");
    if (!fp) return;

    std::fprintf(fp, "# vtk DataFile Version 3.0\n");
    std::fprintf(fp, "Poiseuille 2D Modular\n");
    std::fprintf(fp, "ASCII\n");
    std::fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    std::fprintf(fp, "DIMENSIONS %d %d 1\n", NX, NY);
    std::fprintf(fp, "ORIGIN 0 0 0\n");
    std::fprintf(fp, "SPACING 1 1 1\n");
    std::fprintf(fp, "POINT_DATA %d\n", NUM_CELLS);

    std::fprintf(fp, "SCALARS density double 1\n");
    std::fprintf(fp, "LOOKUP_TABLE default\n");
    for (int i = 0; i < NUM_CELLS; i++) std::fprintf(fp, "%.6f\n", rho[i]);

    std::fprintf(fp, "VECTORS velocity double\n");
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int idx = get_idx(x, y);
            if (y == 0 || y == NY - 1) {
                std::fprintf(fp, "0.0 0.0 0.0\n"); // Paredes em repouso
            } else {
                std::fprintf(fp, "%.6f %.6f 0.0\n", ux[idx], uy[idx]);
            }
        }
    }
    std::fclose(fp);
}