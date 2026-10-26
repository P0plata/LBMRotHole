#ifndef D2Q9_HPP
#define D2Q9_HPP

namespace D2Q9 {
    // Ordem das 8 direções móveis (f1_ptr / f2_ptr):
    // 0:E, 1:N, 2:W, 3:S, 4:NE, 5:NW, 6:SW, 7:SE
    constexpr int cx[8]  = { 1,  0, -1,  0,  1, -1, -1,  1};
    constexpr int cy[8]  = { 0,  1,  0, -1,  1,  1, -1, -1};
    
    // Direção oposta para o Bounce-Back (ex: oposto de 0 [E] é 2 [W])
    constexpr int opp[8] = { 2,  3,  0,  1,  6,  7,  4,  5};

    constexpr double w0 = 4.0 / 9.0;
    constexpr double w[8] = {
        1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };

    constexpr inline double compute_feq(int d, double rho, double ux, double uy, double u_sq) {
        double cu = cx[d] * ux + cy[d] * uy;
        return w[d] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq);
    }

    constexpr inline double compute_feq0(double rho, double u_sq) {
        return w0 * rho * (1.0 - 1.5 * u_sq);
    }

    constexpr inline size_t src_idx(size_t idx, int nx, int d) {
        return idx - (cx[d] + cy[d] * nx);
    }
    
    /*
        indice periodico
        int src_x = (x - D2Q9::cx[d] + nx) % nx;
        int src_y = (y - D2Q9::cy[d] + ny) % ny;
        src_idx = static_cast<size_t>(src_y) * nx + src_x;
    
    */

    // src_idx_periodico
    constexpr inline size_t src_idx_p(int x, int y, int nx, int ny, int d){
        return static_cast<size_t>((x - cx[d] + nx) % nx) + (static_cast<size_t>((y - cy[d] + ny) % ny) * nx);
    }

}

#endif