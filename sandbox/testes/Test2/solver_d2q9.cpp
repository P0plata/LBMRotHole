#include "solver_d2q9.hpp"

namespace SolverD2Q9 {

// ============================================================================
// 1. AS 3 FUNÇÕES DE COLISÃO SEPARADAS (static  -> Zero overhead)
// ============================================================================
template <bool UseForces>
LBM_INLINE  void collide_bgk(double& f0, double f[8], double rho, double ux, double uy, double fx, double fy, double omega) noexcept {
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

    if constexpr (UseForces) {
        const double force_factor = (1.0- 0.5 * omega);
        const double u_dot_F = ux * fx + uy * fy;
        f0   += D2Q9::compute_S0_guo(u_dot_F, force_factor);
        f[0] += D2Q9::compute_S_guo(0, ux, uy, fx, fy, u_dot_F, force_factor);
        f[1] += D2Q9::compute_S_guo(1, ux, uy, fx, fy, u_dot_F, force_factor);
        f[2] += D2Q9::compute_S_guo(2, ux, uy, fx, fy, u_dot_F, force_factor);
        f[3] += D2Q9::compute_S_guo(3, ux, uy, fx, fy, u_dot_F, force_factor);
        f[4] += D2Q9::compute_S_guo(4, ux, uy, fx, fy, u_dot_F, force_factor);
        f[5] += D2Q9::compute_S_guo(5, ux, uy, fx, fy, u_dot_F, force_factor);
        f[6] += D2Q9::compute_S_guo(6, ux, uy, fx, fy, u_dot_F, force_factor);
        f[7] += D2Q9::compute_S_guo(7, ux, uy, fx, fy, u_dot_F, force_factor);
    }

}

template <bool UseForces>
LBM_INLINE void collide_trt(double& f0, double f[8], double rho, double ux, double uy,
                               double fx, double fy, double omega_p, double omega_m) noexcept {
    const double u_sq = ux * ux + uy * uy;
    
    f0 -= omega_p * (f0 - D2Q9::compute_feq0(rho, u_sq));

    // constexpr int pairs[4][2] = {{0, 2}, {1, 3}, {4, 6}, {5, 7}};
    // for (int p = 0; p < 4; ++p) {
    //     const int d     = pairs[p][0];
    //     const int d_opp = pairs[p][1];

    //     const double feq_d   = D2Q9::compute_feq(d, rho, ux, uy, u_sq);
    //     const double feq_opp = D2Q9::compute_feq(d_opp, rho, ux, uy, u_sq);

    //     const double f_plus    = 0.5 * (f[d] + f[d_opp]);
    //     const double f_minus   = 0.5 * (f[d] - f[d_opp]);
    //     const double feq_plus  = 0.5 * (feq_d + feq_opp);
    //     const double feq_minus = 0.5 * (feq_d - feq_opp);

    //     const double corr_p = omega_p * (f_plus - feq_plus);
    //     const double corr_m = omega_m * (f_minus - feq_minus);

    //     f[d]     -= (corr_p + corr_m);
    //     f[d_opp] -= (corr_p - corr_m);
    // }

    // Termos pares (0,2) (1,3) (4,6) (5,7)

    const double f_plus_02 = D2Q9::compute_f_plus_trt(0, f);
    const double f_plus_13 = D2Q9::compute_f_plus_trt(1, f);
    const double f_plus_46 = D2Q9::compute_f_plus_trt(4, f);
    const double f_plus_57 = D2Q9::compute_f_plus_trt(5, f);

    const double f_minus_02 = D2Q9::compute_f_minus_trt(0, f);
    const double f_minus_13 = D2Q9::compute_f_minus_trt(1, f);
    const double f_minus_46 = D2Q9::compute_f_minus_trt(4, f);
    const double f_minus_57 = D2Q9::compute_f_minus_trt(5, f);

    const double feq_plus_02 = D2Q9::compute_feq_plus_trt(0, rho, ux, uy, u_sq);
    const double feq_plus_13 = D2Q9::compute_feq_plus_trt(1, rho, ux, uy, u_sq);
    const double feq_plus_46 = D2Q9::compute_feq_plus_trt(4, rho, ux, uy, u_sq);
    const double feq_plus_57 = D2Q9::compute_feq_plus_trt(5, rho, ux, uy, u_sq);

    const double feq_minus_02 = D2Q9::compute_feq_minus_trt(0, rho, ux, uy, u_sq);
    const double feq_minus_13 = D2Q9::compute_feq_minus_trt(1, rho, ux, uy, u_sq);
    const double feq_minus_46 = D2Q9::compute_feq_minus_trt(4, rho, ux, uy, u_sq);
    const double feq_minus_57 = D2Q9::compute_feq_minus_trt(5, rho, ux, uy, u_sq);

    const double corr_plus_02 = omega_p * (f_plus_02 - feq_plus_02);
    const double corr_plus_13 = omega_p * (f_plus_13 - feq_plus_13);
    const double corr_plus_46 = omega_p * (f_plus_46 - feq_plus_46);
    const double corr_plus_57 = omega_p * (f_plus_57 - feq_plus_57);

    const double corr_minus_02 = omega_m * (f_minus_02 - feq_minus_02);
    const double corr_minus_13 = omega_m * (f_minus_13 - feq_minus_13);
    const double corr_minus_46 = omega_m * (f_minus_46 - feq_minus_46);
    const double corr_minus_57 = omega_m * (f_minus_57 - feq_minus_57);

    f[0] -= corr_plus_02 + corr_minus_02;
    f[1] -= corr_plus_13 + corr_minus_13;
    f[2] -= corr_plus_02 - corr_minus_02;
    f[3] -= corr_plus_13 - corr_minus_13;
    f[4] -= corr_plus_46 + corr_minus_46;
    f[5] -= corr_plus_57 + corr_minus_57;
    f[6] -= corr_plus_46 - corr_minus_46;
    f[7] -= corr_plus_57 - corr_minus_57;


    if constexpr (UseForces) {
        const double force_factor_p = (1.0 - 0.5 * omega_p);
        const double force_factor_m = (1.0 - 0.5 * omega_m);
        const double u_dot_F = ux * fx + uy * fy;

        f0   += D2Q9::compute_S0_guo(u_dot_F, force_factor_p);
        
        const double force_plus_02 = D2Q9::compute_S_guo_plus_trt(0, ux, uy, fx, fy, u_dot_F, force_factor_p);
        const double force_plus_13 = D2Q9::compute_S_guo_plus_trt(1, ux, uy, fx, fy, u_dot_F, force_factor_p);
        const double force_plus_46 = D2Q9::compute_S_guo_plus_trt(4, ux, uy, fx, fy, u_dot_F, force_factor_p);
        const double force_plus_57 = D2Q9::compute_S_guo_plus_trt(5, ux, uy, fx, fy, u_dot_F, force_factor_p);

        const double force_minus_02 = D2Q9::compute_S_guo_minus_trt(0, fx, fy, force_factor_m);
        const double force_minus_13 = D2Q9::compute_S_guo_minus_trt(1, fx, fy, force_factor_m);
        const double force_minus_46 = D2Q9::compute_S_guo_minus_trt(4, fx, fy, force_factor_m);
        const double force_minus_57 = D2Q9::compute_S_guo_minus_trt(5, fx, fy, force_factor_m);
        
        f[0] += force_plus_02 + force_minus_02;
        f[1] += force_plus_13 + force_minus_13;
        f[2] += force_plus_02 - force_minus_02;
        f[3] += force_plus_13 - force_minus_13;
        f[4] += force_plus_46 + force_minus_46;
        f[5] += force_plus_57 + force_minus_57;
        f[6] += force_plus_46 - force_minus_46;
        f[7] += force_plus_57 - force_minus_57;
    }


}



/* =======================================================================================
 * OTIMIZAÇÃO DO OPERADOR MRT COMPLETO ATÉ O(u^4) (Lallemand & Luo, 2000 + Guo, 2002 - D2Q9)
 * =======================================================================================
 * Este solver evita a multiplicação pesada de matrizes 9x9 (M e M^-1) no ciclo principal,
 * aplicando a ortogonalidade de Gram-Schmidt e simetrias de paridade.
 *
 * 1. ORDEM DAS POPULAÇÕES NO SOLVER:
 *    - f0:       Repouso (0, 0)
 *    - Cardinais (f_card): f[0]=E(1,0), f[1]=N(0,1), f[2]=W(-1,0), f[3]=S(0,-1)
 *    - Diagonais (f_diag): f[4]=NE(1,1), f[5]=NW(-1,1), f[6]=SW(-1,-1), f[7]=SE(1,-1)
 *
 * 2. TRATAMENTO DOS 3 MOMENTOS CONSERVADOS (m0, m3, m5):
 *    - m0 (rho): Massa/Densidade      => Norma(M0)^2 = 9. Como m0 = m0_eq e Psi_0 = 0, dm0 = 0 sempre.
 *    - m3 (jx):  Momento linear em X  => Norma(M3)^2 = 6.
 *    - m5 (jy):  Momento linear em Y  => Norma(M5)^2 = 6.
 *    - SEM FORÇAS (UseForces = false): m3 = m3_eq e m5 = m5_eq, logo dm3 = dm5 = 0 (poupam-se FLOPs).
 *    - COM FORÇAS (UseForces = true):  Como a velocidade macroscópica já inclui o shift de Guo
 *      (rho*u = j + 0.5*F), temos (m3 - m3_eq) = -0.5*Fx. Na equação de colisão MRT:
 *      Delta_m3 = s3*(-0.5*Fx) - (1 - 0.5*s3)*Fx = -Fx (independente de s3!).
 *      Dividindo pela norma^2 (6), resulta em dm3 = -Fx / 6 e dm5 = -Fy / 6.
 *
 * 3. CÁLCULO DOS 6 MOMENTOS NÃO-CONSERVADOS ATIVOS (Matriz M Ortogonal):
 *    - m1 (e):   Energia cinética           => Norma(M1)^2 = 36
 *    - m2 (eps): Energia ao quadrado        => Norma(M2)^2 = 36  [Inclui termo 4ª ordem: +9*rho*ux^2*uy^2]
 *    - m4 (qx):  Fluxo de energia em X      => Norma(M4)^2 = 12  [Inclui termo 3ª ordem: +3*rho*ux*uy^2]
 *    - m6 (qy):  Fluxo de energia em Y      => Norma(M6)^2 = 12  [Inclui termo 3ª ordem: +3*rho*uy*ux^2]
 *    - m7 (pxx): Tensão normal (X vs Y)     => Norma(M7)^2 = 4
 *    - m8 (pxy): Tensão de corte (Diagonal) => Norma(M8)^2 = 4
 *
 * 4. CORREÇÃO DE INVARIÂNCIA DE GALILEU (Termos Cúbicos O(u^3) e Quárticos O(u^4)):
 *    - Enquanto c_ix^3 = c_ix na grelha D2Q9 impede termos puros ux^3, os momentos cruzados
 *      c_ix*c_iy^2 e c_ix^2*c_iy^2 são independentes e retêm os termos de Maxwell-Boltzmann
 *      (3*rho*ux*uy^2, 3*rho*uy*ux^2 e 9*rho*ux^2*uy^2) e as respetivas derivadas de força
 *      Psi_k = (F / rho) . grad_u(m_k^eq), aumentando a estabilidade a números de Mach elevados.
 *
 * 5. OS FATORES DE NORMALIZAÇÃO (1/36, 1/12, 1/6, 1/4) E ATUALIZAÇÃO REVERSA (M^-1):
 *    - Como as linhas de M são mutuamente ortogonais, M * M^T = diag(9, 36, 36, 6, 12, 6, 12, 4, 4).
 *    - Assim, a matriz inversa exata é simplesmente M^-1 = M^T * diag(1/Norma^2).
 *    - Fundimos a relaxação s_k*(m - meq), a fonte (1 - 0.5*s_k)*Psi_k e (1/Norma^2) num só "dm_k".
 *    - A projeção final f^* = f - M^T * dm aplica diretamente as colunas de M a cada "dm_k".
 * ======================================================================================= */
template <bool UseForces>
LBM_INLINE void collide_mrt(double& f0, double f[8], double rho, double ux, double uy,
                            double fx, double fy, const double s_mrt[9]) noexcept {
    const double ux2    = ux * ux;
    const double uy2    = uy * uy;
    const double u_sq   = ux2 + uy2;
    const double ux2uy2 = ux2 * uy2;

    const double f_card = f[0] + f[1] + f[2] + f[3];
    const double f_diag = f[4] + f[5] + f[6] + f[7];

    // 1. Momentos atuais não-conservados (m = M * f)
    const double m1 = -4.0 * f0 - f_card + 2.0 * f_diag;
    const double m2 =  4.0 * f0 - 2.0 * f_card + f_diag;
    const double m4 = -2.0 * (f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7]);
    const double m6 = -2.0 * (f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7]);
    const double m7 = (f[0] - f[1] + f[2] - f[3]);
    const double m8 = (f[4] - f[5] + f[6] - f[7]);

    // 2. Momentos de equilíbrio completos até O(u^4) (Maxwell-Boltzmann fatorizado)
    const double m1_eq = rho * (-2.0 + 3.0 * u_sq);
    const double m2_eq = rho * ( 1.0 - 3.0 * u_sq + 9.0 * ux2uy2);
    const double m4_eq = rho * (-ux  + 3.0 * ux * uy2);
    const double m6_eq = rho * (-uy  + 3.0 * uy * ux2);
    const double m7_eq = rho * (ux2 - uy2);
    const double m8_eq = rho * (ux * uy);

    // 3. Relaxação no espaço de momentos: s_k * (m_k - m_k^eq)
    double dm1 = s_mrt[1] * (m1 - m1_eq);
    double dm2 = s_mrt[2] * (m2 - m2_eq);
    double dm4 = s_mrt[4] * (m4 - m4_eq);
    double dm6 = s_mrt[6] * (m6 - m6_eq);
    double dm7 = s_mrt[7] * (m7 - m7_eq);
    double dm8 = s_mrt[8] * (m8 - m8_eq);

    // 4. Inclusão das Forças de Guo completas no espaço de momentos: -(1 - 0.5*s_k) * Psi_k
    double dm3 = 0.0;
    double dm5 = 0.0;

    if constexpr (UseForces) {
        const double u_dot_F = ux * fx + uy * fy;

        const double psi1 =  6.0 * u_dot_F;
        const double psi2 = -6.0 * u_dot_F + 18.0 * (ux * fx * uy2 + ux2 * uy * fy);
        const double psi4 = -fx + 3.0 * (fx * uy2 + 2.0 * ux * uy * fy);
        const double psi6 = -fy + 3.0 * (fy * ux2 + 2.0 * uy * ux * fx);
        const double psi7 =  2.0 * (ux * fx - uy * fy);
        const double psi8 =  ux * fy + uy * fx;

        dm1 -= (1.0 - 0.5 * s_mrt[1]) * psi1;
        dm2 -= (1.0 - 0.5 * s_mrt[2]) * psi2;
        dm3  = (-1.0 / 6.0) * fx; // Já inclui a divisão pela norma^2 (6) da linha m3 (jx)
        dm4 -= (1.0 - 0.5 * s_mrt[4]) * psi4;
        dm5  = (-1.0 / 6.0) * fy; // Já inclui a divisão pela norma^2 (6) da linha m5 (jy)
        dm6 -= (1.0 - 0.5 * s_mrt[6]) * psi6;
        dm7 -= (1.0 - 0.5 * s_mrt[7]) * psi7;
        dm8 -= (1.0 - 0.5 * s_mrt[8]) * psi8;
    }

    // 5. Divisão pela norma ao quadrado de cada linha (M^-1 = M^T * diag(1 / Norma^2))
    dm1 *= (1.0 / 36.0);
    dm2 *= (1.0 / 36.0);
    dm4 *= (1.0 / 12.0);
    dm6 *= (1.0 / 12.0);
    dm7 *= (1.0 / 4.0);
    dm8 *= (1.0 / 4.0);

    // 6. Transformação inversa para o espaço de populações (f^* = f - M^T * dm)
    f0   -= ( -4.0 * dm1 + 4.0 * dm2 );
    f[0] -= ( -dm1 - 2.0 * dm2 + dm3 - 2.0 * dm4 + dm7 );
    f[1] -= ( -dm1 - 2.0 * dm2 + dm5 - 2.0 * dm6 - dm7 );
    f[2] -= ( -dm1 - 2.0 * dm2 - dm3 + 2.0 * dm4 + dm7 );
    f[3] -= ( -dm1 - 2.0 * dm2 - dm5 + 2.0 * dm6 - dm7 );
    f[4] -= (  2.0 * dm1 + dm2 + dm3 + dm4 + dm5 + dm6 + dm8 );
    f[5] -= (  2.0 * dm1 + dm2 - dm3 - dm4 + dm5 + dm6 - dm8 );
    f[6] -= (  2.0 * dm1 + dm2 - dm3 - dm4 - dm5 - dm6 + dm8 );
    f[7] -= (  2.0 * dm1 + dm2 + dm3 + dm4 - dm5 - dm6 - dm8 );
}

// Disparo em tempo de compilação para a função de colisão correta
template <CollisionType ColType, bool UseForces>
LBM_INLINE  void dispatch_collision(double& f0, double f[8], double rho, double ux, double uy,
                                      double fx, double fy, double omega_bgk, double omega_trt, const double s_mrt[9]) noexcept {
                                        
    if constexpr (ColType == CollisionType::BGK) {
        collide_bgk<UseForces>(f0, f, rho, ux, uy, fx, fy, omega_bgk);
    } else if constexpr (ColType == CollisionType::TRT) {
        collide_trt<UseForces>(f0, f, rho, ux, uy, fx, fy, omega_bgk, omega_trt);
    } else if constexpr (ColType == CollisionType::MRT) {
        collide_mrt<UseForces>(f0, f, rho, ux, uy, fx, fy, s_mrt);
    }
}

// ============================================================================
// 2. PULL + COLLIDE POR NÓ DE FLUIDO
// ============================================================================

// template <CollisionType ColType, bool IsInterior, bool UseForces, bool UseNTStore = false>
// static  void pull_and_collide_node(int x, int y, int nx, int ny, size_t pitch, Grid& grid,
//                                          double omega_bgk, double omega_trt, const double s_mrt[9], int step,  Custom::CustomForces& forces_struct) noexcept {
//     const size_t idx = static_cast<size_t>(y) * nx + x;
    
    
//     if constexpr (!IsInterior) {
//         const CellType cell = grid.flag_field_ptr[idx];
//         if (cell != CellType::FLUID && cell != CellType::SPONGE) return;
//     }

//     double* __restrict__ f0_p       = grid.f0_ptr;
//     const double* __restrict__ f1_p = grid.f1_ptr;
//     double* __restrict__ f2_p       = grid.f2_ptr;

//     //double feq;
//     double f0 = f0_p[idx];
//     double f[8];

//     /*for (int d = 0; d < 8; ++d) {
//         size_t src_idx;
//         if constexpr (IsInterior) {
//             src_idx = idx - (D2Q9::cx[d] + D2Q9::cy[d] * nx);
//         } else {
//             int src_x = (x - D2Q9::cx[d] + nx) % nx;
//             int src_y = (y - D2Q9::cy[d] + ny) % ny;
//             src_idx = static_cast<size_t>(src_y) * nx + src_x;
//         }
//         f[d] = grid.f1_ptr[d * pitch + src_idx];
//     }*/

//     if constexpr (IsInterior) {
        
//         f[0] = f1_p[0 * pitch + D2Q9::src_idx(idx,nx,0)];
//         f[1] = f1_p[1 * pitch + D2Q9::src_idx(idx,nx,1)];
//         f[2] = f1_p[2 * pitch + D2Q9::src_idx(idx,nx,2)];
//         f[3] = f1_p[3 * pitch + D2Q9::src_idx(idx,nx,3)];
//         f[4] = f1_p[4 * pitch + D2Q9::src_idx(idx,nx,4)];
//         f[5] = f1_p[5 * pitch + D2Q9::src_idx(idx,nx,5)];
//         f[6] = f1_p[6 * pitch + D2Q9::src_idx(idx,nx,6)];
//         f[7] = f1_p[7 * pitch + D2Q9::src_idx(idx,nx,7)];            

//     } else {

//         f[0] = f1_p[0 * pitch + D2Q9::src_idx_p(x,y,nx,ny,0)];
//         f[1] = f1_p[1 * pitch + D2Q9::src_idx_p(x,y,nx,ny,1)];
//         f[2] = f1_p[2 * pitch + D2Q9::src_idx_p(x,y,nx,ny,2)];
//         f[3] = f1_p[3 * pitch + D2Q9::src_idx_p(x,y,nx,ny,3)];
//         f[4] = f1_p[4 * pitch + D2Q9::src_idx_p(x,y,nx,ny,4)];
//         f[5] = f1_p[5 * pitch + D2Q9::src_idx_p(x,y,nx,ny,5)];
//         f[6] = f1_p[6 * pitch + D2Q9::src_idx_p(x,y,nx,ny,6)];
//         f[7] = f1_p[7 * pitch + D2Q9::src_idx_p(x,y,nx,ny,7)];    
    
//     }

//     const double rho = f0 + f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
//     const double inv_rho = 1.0 / rho;

//     double ux, uy;
//     double fx = 0.0;
//     double fy = 0.0;

//     if constexpr (UseForces) {
//         // Calcula as forças
//         fx = Custom::compute_custom_fx(x, y, step, forces_struct);
//         fy = Custom::compute_custom_fy(x, y, step, forces_struct);

//         // Velocidade com o Shift de Guo
//         ux = (((f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7])) + 0.5 * fx) * inv_rho;
//         uy = (((f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7])) + 0.5 * fy) * inv_rho;
//     } else {
//         // Velocidade normal sem forças
//         ux = ((f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7])) * inv_rho;
//         uy = ((f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7])) * inv_rho;
//     }

//     dispatch_collision<ColType, UseForces>(f0, f, rho, ux, uy, fx, fy, omega_bgk, omega_trt, s_mrt);

//     // Avisar o compilador que os ponteiros base estão alinhados a 64 bytes (Cache Line / AVX-512)
//     double* __restrict__ f0_aligned = static_cast<double*>(__builtin_assume_aligned(f0_p, 64));
//     double* __restrict__ f2_aligned = static_cast<double*>(__builtin_assume_aligned(f2_p, 64));

//     D2Q9::store_population<UseNTStore>(&f0_aligned[idx], f0);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[0 * pitch + idx], f[0]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[1 * pitch + idx], f[1]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[2 * pitch + idx], f[2]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[3 * pitch + idx], f[3]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[4 * pitch + idx], f[4]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[5 * pitch + idx], f[5]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[6 * pitch + idx], f[6]);
//     D2Q9::store_population<UseNTStore>(&f2_aligned[7 * pitch + idx], f[7]);
// }


template <CollisionType ColType, bool IsInterior, bool UseForces>
LBM_INLINE void pull_and_collide_node(int x, int y, int nx, int ny, size_t pitch, Grid& grid,
                                      double omega_bgk, double omega_trt, const double s_mrt[9],
                                      int step, Custom::CustomForces& forces_struct) noexcept {
    const size_t idx = static_cast<size_t>(y) * nx + x;

    if constexpr (!IsInterior) {
        const CellType cell = grid.flag_field_ptr[idx];
        if (cell != CellType::FLUID && cell != CellType::SPONGE) return;
    }

    double* __restrict__ f0_p       = grid.f0_ptr;
    const double* __restrict__ f1_p = grid.f1_ptr;
    double* __restrict__ f2_p       = grid.f2_ptr;

    double f0 = f0_p[idx];
    double f[8];

    if constexpr (IsInterior) {
        f[0] = f1_p[0 * pitch + D2Q9::src_idx(idx, nx, 0)];
        f[1] = f1_p[1 * pitch + D2Q9::src_idx(idx, nx, 1)];
        f[2] = f1_p[2 * pitch + D2Q9::src_idx(idx, nx, 2)];
        f[3] = f1_p[3 * pitch + D2Q9::src_idx(idx, nx, 3)];
        f[4] = f1_p[4 * pitch + D2Q9::src_idx(idx, nx, 4)];
        f[5] = f1_p[5 * pitch + D2Q9::src_idx(idx, nx, 5)];
        f[6] = f1_p[6 * pitch + D2Q9::src_idx(idx, nx, 6)];
        f[7] = f1_p[7 * pitch + D2Q9::src_idx(idx, nx, 7)];
    } else {
        f[0] = f1_p[0 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 0)];
        f[1] = f1_p[1 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 1)];
        f[2] = f1_p[2 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 2)];
        f[3] = f1_p[3 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 3)];
        f[4] = f1_p[4 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 4)];
        f[5] = f1_p[5 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 5)];
        f[6] = f1_p[6 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 6)];
        f[7] = f1_p[7 * pitch + D2Q9::src_idx_p(x, y, nx, ny, 7)];
    }

    const double rho = f0 + f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7];
    const double inv_rho = 1.0 / rho;

    double ux, uy;
    double fx = 0.0;
    double fy = 0.0;

    if constexpr (UseForces) {
        fx = Custom::compute_custom_fx(x, y, step, forces_struct);
        fy = Custom::compute_custom_fy(x, y, step, forces_struct);
        ux = (((f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7])) + 0.5 * fx) * inv_rho;
        uy = (((f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7])) + 0.5 * fy) * inv_rho;
    } else {
        ux = ((f[0] - f[2]) + (f[4] - f[5] - f[6] + f[7])) * inv_rho;
        uy = ((f[1] - f[3]) + (f[4] + f[5] - f[6] - f[7])) * inv_rho;
    }

    dispatch_collision<ColType, UseForces>(f0, f, rho, ux, uy, fx, fy, omega_bgk, omega_trt, s_mrt);

    f0_p[idx] = f0;
    f2_p[0 * pitch + idx] = f[0];
    f2_p[1 * pitch + idx] = f[1];
    f2_p[2 * pitch + idx] = f[2];
    f2_p[3 * pitch + idx] = f[3];
    f2_p[4 * pitch + idx] = f[4];
    f2_p[5 * pitch + idx] = f[5];
    f2_p[6 * pitch + idx] = f[6];
    f2_p[7 * pitch + idx] = f[7];
}


// ============================================================================
// 3. FRONTEIRAS (Esponja, Zou-He Inlet/Outlet e Bouzidi/Bounce-Back)
// ============================================================================
static  void apply_sponges(Grid& grid) noexcept {
    const size_t pitch = grid.plane_pitch;
    const size_t n_sponges = grid.sponge_nodes.size();

    double* __restrict__ f0_p       = grid.f0_ptr;
    double* __restrict__ f2_p       = grid.f2_ptr;

    #pragma omp for schedule(static) nowait
    for (size_t i = 0; i < n_sponges; ++i) {
        const SpongeNode& sp = grid.sponge_nodes[i];
        const size_t idx = sp.idx;
        const double sigma = sp.sigma;
        const double oms = 1.0 - sigma;
        const double u_sq = static_cast<double>(sp.vx_ref) * sp.vx_ref + static_cast<double>(sp.vy_ref) * sp.vy_ref;

        f0_p[idx] = oms * f0_p[idx] + sigma * D2Q9::compute_feq0(sp.rho_ref, u_sq);
        for (int d = 0; d < 8; ++d) {
            double feq_ref = D2Q9::compute_feq(d, sp.rho_ref, sp.vx_ref, sp.vy_ref, u_sq);
            f2_p[d * pitch + idx] = oms * f2_p[d * pitch + idx] + sigma * feq_ref;
        }
    }
}

template <CollisionType ColType, bool IsInlet, bool UseForces>
static  void apply_open_boundaries(Grid& grid, const std::vector<OpenBoundaryNode>& nodes,
                                         double omega_bgk, double omega_trt, const double s_mrt[9],
                                         bool periodic_x, bool periodic_y, int step, Custom::CustomForces& forces_struct) noexcept {
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
            if (periodic_x) src_x = (src_x + nx) % nx;
            if (periodic_y) src_y = (src_y + ny) % ny;
            if (src_x >= 0 && src_x < nx && src_y >= 0 && src_y < ny) {
                f[d] = grid.f1_ptr[d * pitch + (static_cast<size_t>(src_y) * nx + src_x)];
            } else {
                f[d] = grid.f1_ptr[d * pitch + idx];
            }
        }

        double fx = 0.0, fy = 0.0;
        if constexpr (UseForces) {
            fx = Custom::compute_custom_fx(x, y, step, forces_struct);
            fy = Custom::compute_custom_fy(x, y, step, forces_struct);
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

        // Incorpora o shift de Guo na velocidade antes da colisão na fronteira
        if constexpr (UseForces) {
            ux += 0.5 * fx / rho;
            uy += 0.5 * fy / rho;
        }

        dispatch_collision<ColType, UseForces>(f0, f, rho, ux, uy, fx, fy, omega_bgk, omega_trt, s_mrt);

        grid.f0_ptr[idx] = f0;
        for (int d = 0; d < 8; ++d) {
            grid.f2_ptr[d * pitch + idx] = f[d];
        }
    }
}

static  void apply_wall_links(Grid& grid, double* f_target, bool periodic_x, bool periodic_y) noexcept {
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
template <CollisionType ColType, bool UseForces>
void run_simulation(Grid& grid, const Simconfig& cfg, Custom::CustomForces& forces_struct) noexcept {
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
                    double t0 = omp_get_wtime();
                    VTKExporter::export_to_vti<UseForces>(grid, forces_struct, step, cfg.output_prefix);
                    grid.vtk_io_time_sec += (omp_get_wtime() - t0);
                }
            }

            if (step == cfg.max_steps) break;

            // 2. Fused Pull + Collide em todo o domínio fluido

            #pragma omp for schedule(static) nowait
            for (int y = 1; y < ny - 1; ++y) {
                #pragma omp simd
                for (int x = 1; x < nx - 1; ++x) {
                    pull_and_collide_node<ColType, true, UseForces>(x, y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT, step, forces_struct);
                }
            }

            // 2. BORDAS PERIÓDICAS (Só corre se a fronteira for realmente PERIODIC)
            if (per_y) {
                #pragma omp for schedule(static) nowait
                for (int x = 0; x < nx; ++x) {
                    pull_and_collide_node<ColType, false, UseForces>(x, 0,      nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT, step, forces_struct);
                    pull_and_collide_node<ColType, false, UseForces>(x, ny - 1, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT, step, forces_struct);
                }
            }
            if (per_x) {
                #pragma omp for schedule(static) nowait
                for (int y = 1; y < ny - 1; ++y) {
                    pull_and_collide_node<ColType, false, UseForces>(0,      y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT, step, forces_struct);
                    pull_and_collide_node<ColType, false, UseForces>(nx - 1, y, nx, ny, pitch, grid, omega_bgk, omega_trt, cfg.s_MRT, step, forces_struct);
                }
            }

            // 3. Aplicar Fronteiras (Esponjas, Inlets, Outlets correm com 'nowait' e sincronizam no fim de wall_links)
            //#pragma omp barrier // Sincroniza tudo antes de aplicar as listas de fronteira

            apply_open_boundaries<ColType, true,  UseForces>(grid, grid.inlet_nodes,  omega_bgk, omega_trt, cfg.s_MRT, per_x, per_y, step, forces_struct);
            apply_open_boundaries<ColType, false, UseForces>(grid, grid.outlet_nodes, omega_bgk, omega_trt, cfg.s_MRT, per_x, per_y, step, forces_struct);
            
            #pragma omp barrier

            apply_sponges(grid);
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
// BGK (4 combinações)
template void run_simulation<CollisionType::BGK, true>(Grid&, const Simconfig&, Custom::CustomForces&);
template void run_simulation<CollisionType::BGK, false>(Grid&, const Simconfig&, Custom::CustomForces&);

// TRT (4 combinações)
template void run_simulation<CollisionType::TRT, true>(Grid&, const Simconfig&, Custom::CustomForces&);
template void run_simulation<CollisionType::TRT, false>(Grid&, const Simconfig&, Custom::CustomForces&);

// MRT (4 combinações)
template void run_simulation<CollisionType::MRT, true>(Grid&, const Simconfig&, Custom::CustomForces&);
template void run_simulation<CollisionType::MRT, false>(Grid&, const Simconfig&, Custom::CustomForces&);


} // namespace SolverD2Q9


