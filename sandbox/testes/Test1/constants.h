#ifndef CONSTANTS_H
#define CONSTANTS_H

// Grelha 2D
inline constexpr int NX = 64;
inline constexpr int NY = 32;
inline constexpr int NUM_CELLS = NX * NY;

// Parâmetros físicos
inline constexpr double TAU = 1.0;            // Tempo de relaxamento (viscosidade nu = 1/6)
inline constexpr double OMEGA = 1.0 / TAU;     // Frequência de colisão
inline constexpr double FORCE_X = 1e-4;       // Aceleração externa constante em X

// Conjunto de velocidades D2Q9
inline constexpr double w[9] = {
    4.0 / 9.0,
    1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
    1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0
};

inline constexpr int cx[9] = {0,  1,  0, -1,  0,  1, -1, -1,  1};
inline constexpr int cy[9] = {0,  0,  1,  0, -1,  1,  1, -1, -1};

// Direção inversa para bounce-back
inline constexpr int opp[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

#endif