#include "lbm.h"
#include <cstdio>

int main(int argc, char *argv[]) {
    int time_steps = 2000;

    // Buffers das populações
    Grid *f_current = new Grid();
    Grid *f_eq      = new Grid();
    Grid *f_post    = new Grid();
    Grid *f_next    = new Grid();

    // Arrays para guardar os momentos macroscópicos do passo
    double rho[NUM_CELLS];
    double ux[NUM_CELLS];
    double uy[NUM_CELLS];

    initialize(f_current);

    for (int t = 0; t < time_steps; t++) {
        // 1. Momento macroscópico
        compute_macroscopic(f_current, rho, ux, uy);

        // 2. Distribuição de equilíbrio
        compute_equilibrium(rho, ux, uy, f_eq);

        // 3. Operador de colisão (BGK)
        collide_bgk(f_current, f_eq, f_post);

        // 4. Termo de força de corpo
        apply_forcing(ux, uy, f_post);

        // 5. Streaming (interior)
        stream(f_post, f_next);

        // 6. Condições de fronteira (paredes sólidas)
        apply_boundary_conditions(f_post, f_next);

        // 7. Saída de dados periódica
        if (t % 200 == 0 || t == time_steps - 1) {
            export_to_vtk(rho, ux, uy, t);
            std::printf("Passo %d gravado em ficheiro VTK.\n", t);
        }

        // 8. Troca de ponteiros para o passo seguinte
        Grid *temp = f_current;
        f_current = f_next;
        f_next = temp;
    }

    delete f_current;
    delete f_eq;
    delete f_post;
    delete f_next;

    return 0;
}