#ifndef D_Q_HPP
#define D_Q_HPP

#include "types.hpp"
#include <cstddef>



// ============================================================================
// 1. DEFINIÇÃO DAS REDES (Apenas Geometria, Pesos e Velocidade do Som)
//    Convenção obrigatória para Esoteric Twist:
//      - 0 .. NUM_PAIRS-1            : direções positivas (cx >= 0)
//      - NUM_PAIRS .. 2*NUM_PAIRS-1  : direções opostas exatas (c[p + N] = -c[p])
// ============================================================================

struct D3Q19 {
    static constexpr int DIM       = 3;
    static constexpr int Q         = 19;
    static constexpr int NUM_PAIRS = 9;
    static constexpr int Q_MOVING  = 18;

    // Ordem dos 9 pares:
    //  Faces (3):   0:(+1,0,0), 1:(0,+1,0), 2:(0,0,+1)
    //  Arestas XY:  3:(+1,+1,0), 4:(+1,-1,0)
    //  Arestas XZ:  5:(+1,0,+1), 6:(+1,0,-1)
    //  Arestas YZ:  7:(0,+1,+1), 8:(0,+1,-1)
    //  Opostas:     9..17 = - (0..8)
    static constexpr int cx[Q_MOVING] = {
         1,  0,  0,   1,  1,   1,  1,   0,  0,
        -1,  0,  0,  -1, -1,  -1, -1,   0,  0
    };
    static constexpr int cy[Q_MOVING] = {
         0,  1,  0,   1, -1,   0,  0,   1,  1,
         0, -1,  0,  -1,  1,   0,  0,  -1, -1
    };
    static constexpr int cz[Q_MOVING] = {
         0,  0,  1,   0,  0,   1, -1,   1, -1,
         0,  0, -1,   0,  0,  -1,  1,  -1,  1
    };

    // Versões em lbm_real para evitar conversões int->float no loop interno
    static constexpr lbm_real cxr[Q_MOVING] = {
         1,  0,  0,   1,  1,   1,  1,   0,  0,
        -1,  0,  0,  -1, -1,  -1, -1,   0,  0
    };
    static constexpr lbm_real cyr[Q_MOVING] = {
         0,  1,  0,   1, -1,   0,  0,   1,  1,
         0, -1,  0,  -1,  1,   0,  0,  -1, -1
    };
    static constexpr lbm_real czr[Q_MOVING] = {
         0,  0,  1,   0,  0,   1, -1,   1, -1,
         0,  0, -1,   0,  0,  -1,  1,  -1,  1
    };

    static constexpr int opp[Q_MOVING] = {
         9, 10, 11, 12, 13, 14, 15, 16, 17,
         0,  1,  2,  3,  4,  5,  6,  7,  8
    };

    // Pesos exatos D3Q19: w0 = 1/3, faces = 1/18, arestas = 1/36 (Soma = 1.0)
    static constexpr lbm_real wf = static_cast<lbm_real>(1.0 / 18.0);
    static constexpr lbm_real we = static_cast<lbm_real>(1.0 / 36.0);

    static constexpr lbm_real w0 = static_cast<lbm_real>(1.0 / 3.0);
    static constexpr lbm_real w[Q_MOVING] = {
        wf, wf, wf,  we, we,  we, we,  we, we,
        wf, wf, wf,  we, we,  we, we,  we, we
    };

    // Constantes da velocidade do som (cs^2 = 1/3) — todas em lbm_real!
    static constexpr lbm_real cs2          = static_cast<lbm_real>(1.0 / 3.0);
    static constexpr lbm_real inv_cs2      = static_cast<lbm_real>(3.0);  // 1 / cs^2
    static constexpr lbm_real inv_cs4      = static_cast<lbm_real>(9.0);  // 1 / cs^4
    static constexpr lbm_real inv_2cs4     = static_cast<lbm_real>(4.5);  // 1 / (2 * cs^4)
    static constexpr lbm_real half_inv_cs2 = static_cast<lbm_real>(1.5);  // 1 / (2 * cs^2)

    // Normas inversas ao quadrado das 19 linhas ortogonais da matriz M de d'Humières (2002)
    // M * M^T = diag(19, 2394, 252, 10, 40, 10, 40, 10, 40, 36, 72, 12, 24, 4, 4, 4, 8, 8, 8)
    static constexpr lbm_real inv_norm_mrt[Q] = {
        static_cast<lbm_real>(1.0 / 19.0),   // 0:  rho
        static_cast<lbm_real>(1.0 / 2394.0), // 1:  e
        static_cast<lbm_real>(1.0 / 252.0),  // 2:  eps
        static_cast<lbm_real>(1.0 / 10.0),   // 3:  jx
        static_cast<lbm_real>(1.0 / 40.0),   // 4:  qx
        static_cast<lbm_real>(1.0 / 10.0),   // 5:  jy
        static_cast<lbm_real>(1.0 / 40.0),   // 6:  qy
        static_cast<lbm_real>(1.0 / 10.0),   // 7:  jz
        static_cast<lbm_real>(1.0 / 40.0),   // 8:  qz
        static_cast<lbm_real>(1.0 / 36.0),   // 9:  3*pxx
        static_cast<lbm_real>(1.0 / 72.0),   // 10: 3*pi_xx
        static_cast<lbm_real>(1.0 / 12.0),   // 11: pww
        static_cast<lbm_real>(1.0 / 24.0),   // 12: pi_ww
        static_cast<lbm_real>(1.0 / 4.0),    // 13: pxy
        static_cast<lbm_real>(1.0 / 4.0),    // 14: pyz
        static_cast<lbm_real>(1.0 / 4.0),    // 15: pxz
        static_cast<lbm_real>(1.0 / 8.0),    // 16: mx
        static_cast<lbm_real>(1.0 / 8.0),    // 17: my
        static_cast<lbm_real>(1.0 / 8.0)     // 18: mz
    };

    // Taxas de relaxamento padrão de d'Humières (2002)
    // Valor negativo (-1.0) = preencher automaticamente com s_nu = 1 / tau!
    static constexpr lbm_real s_mrt_default[Q] = {
        static_cast<lbm_real>(0.0),  // 0:  rho (conservado)
        static_cast<lbm_real>(1.19), // 1:  e   (viscosidade volumétrica / bulk)
        static_cast<lbm_real>(1.40), // 2:  eps
        static_cast<lbm_real>(0.0),  // 3:  jx  (conservado)
        static_cast<lbm_real>(1.20), // 4:  qx
        static_cast<lbm_real>(0.0),  // 5:  jy  (conservado)
        static_cast<lbm_real>(1.20), // 6:  qy
        static_cast<lbm_real>(0.0),  // 7:  jz  (conservado)
        static_cast<lbm_real>(1.20), // 8:  qz
        static_cast<lbm_real>(-1.0), // 9:  3*pxx -> s_nu = 1/tau
        static_cast<lbm_real>(1.40), // 10: 3*pi_xx
        static_cast<lbm_real>(-1.0), // 11: pww   -> s_nu = 1/tau
        static_cast<lbm_real>(1.40), // 12: pi_ww
        static_cast<lbm_real>(-1.0), // 13: pxy   -> s_nu = 1/tau
        static_cast<lbm_real>(-1.0), // 14: pyz   -> s_nu = 1/tau
        static_cast<lbm_real>(-1.0), // 15: pxz   -> s_nu = 1/tau
        static_cast<lbm_real>(1.98), // 16: mx
        static_cast<lbm_real>(1.98), // 17: my
        static_cast<lbm_real>(1.98)  // 18: mz
    };

};

// Exemplo de como D3Q27 encaixa exatamente na mesma interface (pronto para usar quando quiseres!)
struct D3Q27 {
    static constexpr int DIM       = 3;
    static constexpr int Q         = 27;
    static constexpr int NUM_PAIRS = 13;
    static constexpr int Q_MOVING  = 26;

    // 0..2: Faces | 3..8: Arestas | 9..12: Cantos | 13..25: Opostas exatas
    static constexpr int cx[Q_MOVING] = {
         1,  0,  0,   1,  1,  1,  1,  0,  0,   1,  1,  1,  1,
        -1,  0,  0,  -1, -1, -1, -1,  0,  0,  -1, -1, -1, -1
    };
    static constexpr int cy[Q_MOVING] = {
         0,  1,  0,   1, -1,  0,  0,  1,  1,   1,  1, -1, -1,
         0, -1,  0,  -1,  1,  0,  0, -1, -1,  -1, -1,  1,  1
    };
    static constexpr int cz[Q_MOVING] = {
         0,  0,  1,   0,  0,  1, -1,  1, -1,   1, -1,  1, -1,
         0,  0, -1,   0,  0, -1,  1, -1,  1,  -1,  1, -1,  1
    };

    static constexpr lbm_real cxr[Q_MOVING] = {
         1,  0,  0,   1,  1,  1,  1,  0,  0,   1,  1,  1,  1,
        -1,  0,  0,  -1, -1, -1, -1,  0,  0,  -1, -1, -1, -1
    };
    static constexpr lbm_real cyr[Q_MOVING] = {
         0,  1,  0,   1, -1,  0,  0,  1,  1,   1,  1, -1, -1,
         0, -1,  0,  -1,  1,  0,  0, -1, -1,  -1, -1,  1,  1
    };
    static constexpr lbm_real czr[Q_MOVING] = {
         0,  0,  1,   0,  0,  1, -1,  1, -1,   1, -1,  1, -1,
         0,  0, -1,   0,  0, -1,  1, -1,  1,  -1,  1, -1,  1
    };

    static constexpr int opp[Q_MOVING] = {
        13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
         0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12
    };

    static constexpr lbm_real wf = static_cast<lbm_real>(2.0 / 27.0);
    static constexpr lbm_real we = static_cast<lbm_real>(1.0 / 54.0);
    static constexpr lbm_real wc = static_cast<lbm_real>(1.0 / 216.0);

    static constexpr lbm_real w0 = static_cast<lbm_real>(8.0 / 27.0);
    static constexpr lbm_real w[Q_MOVING] = {
        wf, wf, wf,  we, we, we, we, we, we,  wc, wc, wc, wc,
        wf, wf, wf,  we, we, we, we, we, we,  wc, wc, wc, wc
    };

    static constexpr lbm_real cs2          = static_cast<lbm_real>(1.0 / 3.0);
    static constexpr lbm_real inv_cs2      = static_cast<lbm_real>(3.0);
    static constexpr lbm_real inv_cs4      = static_cast<lbm_real>(9.0);
    static constexpr lbm_real inv_2cs4     = static_cast<lbm_real>(4.5);
    static constexpr lbm_real half_inv_cs2 = static_cast<lbm_real>(1.5);
};

// ============================================================================
// 2. FUNÇÕES GENÉRICAS DE LBM (Funcionam com QUALQUER rede: D3Q19, D3Q27, D3Q7!)
// ============================================================================
namespace LBM {

    // Deslocamento linear COM SINAL (ptrdiff_t!) até ao vizinho na direção d (usado no Esoteric Twist!)
    template <typename Lattice>
    constexpr LBM_INLINE ptrdiff_t neighbor_offset(int d, int stride_x, int stride_y) noexcept {
        return static_cast<ptrdiff_t>(Lattice::cx[d])
             + static_cast<ptrdiff_t>(Lattice::cy[d]) * stride_x
             + static_cast<ptrdiff_t>(Lattice::cz[d]) * stride_x * stride_y;
    }

    // --- EQUILÍBRIO (Maxwell-Boltzmann 2ª ordem, parametrizado por cs^2 da rede!) ---

    template <typename Lattice>
    constexpr LBM_INLINE lbm_real compute_feq0(lbm_real rho, lbm_real u_sq) noexcept {
        constexpr lbm_real ONE = static_cast<lbm_real>(1.0);
        return Lattice::w0 * rho * (ONE - Lattice::half_inv_cs2 * u_sq);
    }

    // Equilíbrio para 1 direção individual d (0 .. Q_MOVING-1)
    template <typename Lattice>
    constexpr LBM_INLINE lbm_real compute_feq(int d, lbm_real rho,
                                              lbm_real ux, lbm_real uy, lbm_real uz,
                                              lbm_real u_sq) noexcept {
        constexpr lbm_real ONE = static_cast<lbm_real>(1.0);
        const lbm_real cu = Lattice::cxr[d] * ux + Lattice::cyr[d] * uy + Lattice::czr[d] * uz;
        return Lattice::w[d] * rho * (ONE + Lattice::inv_cs2 * cu + Lattice::inv_2cs4 * cu * cu - Lattice::half_inv_cs2 * u_sq);
    }

    // Equilíbrio por PAR OPOSTO p (0 .. NUM_PAIRS-1): devolve (feq_pos, feq_neg) com METADE das contas!
    // Perfeito para BGK, TRT e Esoteric Twist!
    template <typename Lattice>
    constexpr LBM_INLINE void compute_feq_pair(int p, lbm_real rho,
                                               lbm_real ux, lbm_real uy, lbm_real uz,
                                               lbm_real u_sq,
                                               lbm_real& feq_pos, lbm_real& feq_neg) noexcept {
        constexpr lbm_real ONE = static_cast<lbm_real>(1.0);
        const lbm_real cu    = Lattice::cxr[p] * ux + Lattice::cyr[p] * uy + Lattice::czr[p] * uz;
        const lbm_real w_rho = Lattice::w[p] * rho;
        const lbm_real even  = ONE - Lattice::half_inv_cs2 * u_sq + Lattice::inv_2cs4 * cu * cu;
        const lbm_real odd   = Lattice::inv_cs2 * cu;

        feq_pos = w_rho * (even + odd);
        feq_neg = w_rho * (even - odd);
    }

    // --- FORÇAS DE GUO (2002) GENÉRICAS (Serve para BGK e TRT!) ---

    template <typename Lattice>
    constexpr LBM_INLINE lbm_real compute_S0_guo(lbm_real u_dot_F, lbm_real force_factor) noexcept {
        return (force_factor * Lattice::w0) * (-Lattice::inv_cs2 * u_dot_F);
    }

    // Termo de Guo para 1 direção individual d
    template <typename Lattice>
    constexpr LBM_INLINE lbm_real compute_S_guo(int d,
                                                lbm_real ux, lbm_real uy, lbm_real uz,
                                                lbm_real fx, lbm_real fy, lbm_real fz,
                                                lbm_real u_dot_F, lbm_real force_factor) noexcept {
        const lbm_real cu = Lattice::cxr[d] * ux + Lattice::cyr[d] * uy + Lattice::czr[d] * uz;
        const lbm_real cF = Lattice::cxr[d] * fx + Lattice::cyr[d] * fy + Lattice::czr[d] * fz;
        return (force_factor * Lattice::w[d]) * (Lattice::inv_cs2 * (cF - u_dot_F) + Lattice::inv_cs4 * cu * cF);
    }

    // Termo de Guo por PAR OPOSTO p (0 .. NUM_PAIRS-1):
    // - No BGK passas: fp = fm = (1 - 0.5 * omega)
    // - No TRT passas: fp = (1 - 0.5 * omega_p), fm = (1 - 0.5 * omega_m)
    template <typename Lattice>
    constexpr LBM_INLINE void compute_S_guo_pair(int p,
                                                 lbm_real ux, lbm_real uy, lbm_real uz,
                                                 lbm_real fx, lbm_real fy, lbm_real fz,
                                                 lbm_real u_dot_F,
                                                 lbm_real fp, lbm_real fm,
                                                 lbm_real& S_pos, lbm_real& S_neg) noexcept {
        const lbm_real cu = Lattice::cxr[p] * ux + Lattice::cyr[p] * uy + Lattice::czr[p] * uz;
        const lbm_real cF = Lattice::cxr[p] * fx + Lattice::cyr[p] * fy + Lattice::czr[p] * fz;

        const lbm_real S_even = (fp * Lattice::w[p]) * (-Lattice::inv_cs2 * u_dot_F + Lattice::inv_cs4 * cu * cF);
        const lbm_real S_odd  = (fm * Lattice::w[p]) * ( Lattice::inv_cs2 * cF);

        S_pos = S_even + S_odd;
        S_neg = S_even - S_odd;
    }

    // Preenche s_mrt[] com os valores padrão da rede (substituindo -1.0 por 1/tau)
    // No initializer.cpp, se houver tau_MRT_filename (ou valores custom no JSON), basta sobrescrever por cima!
    template <typename Lattice>
    constexpr LBM_INLINE void fill_default_s_mrt(lbm_real s_mrt[Lattice::Q], lbm_real tau) noexcept {
        const lbm_real s_nu = static_cast<lbm_real>(1.0) / tau;
        for (int k = 0; k < Lattice::Q; ++k) {
            s_mrt[k] = (Lattice::s_mrt_default[k] < static_cast<lbm_real>(0.0))
                     ? s_nu
                     : Lattice::s_mrt_default[k];
        }
    }

} // namespace LBM

#endif // D_Q_HPP