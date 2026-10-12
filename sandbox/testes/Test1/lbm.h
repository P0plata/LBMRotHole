#ifndef LBM_H
#define LBM_H

#include "constants.h"

struct Grid {
    double f[9][NUM_CELLS];
};

inline int get_idx(int x, int y) {
    return y * NX + x;
}

// 1. Inicialização
void initialize(Grid *g);

// 2. Cálculo dos Momentos Macroscópicos (rho, u)
void compute_macroscopic(const Grid *g, double rho[NUM_CELLS], double ux[NUM_CELLS], double uy[NUM_CELLS]);

// 3. Cálculo da Distribuição de Equilíbrio (f^eq)
void compute_equilibrium(const double rho[NUM_CELLS], const double ux[NUM_CELLS], const double uy[NUM_CELLS], Grid *f_eq);

// 4. Colisão Local (ex: BGK)
void collide_bgk(const Grid *f_in, const Grid *f_eq, Grid *f_out);

// 5. Adição de Termo de Força (Corpo / Gravidade)
void apply_forcing(const double ux[NUM_CELLS], const double uy[NUM_CELLS], Grid *f_post);

// 6. Propagação (Streaming no interior)
void stream(const Grid *f_src, Grid *f_dst);

// 7. Condições de Fronteira (Paredes / Bounce-Back)
void apply_boundary_conditions(const Grid *f_post_collision, Grid *f_dst);

// 8. Escrita de Ficheiro
void export_to_vtk(const double rho[NUM_CELLS], const double ux[NUM_CELLS], const double uy[NUM_CELLS], int step);

#endif