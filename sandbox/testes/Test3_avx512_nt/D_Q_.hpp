#ifndef D2Q9_HPP
#define D2Q9_HPP

#include <random>
#include <stdexcept>
#include <immintrin.h>
#include <cstdint>
#include <type_traits>

#define LBM_INLINE __attribute__((always_inline)) inline

namespace D2Q9 {
    constexpr int cx[8]     = { 1,  0, -1,  0,  1, -1, -1,  1};
    constexpr int cy[8]     = { 0,  1,  0, -1,  1,  1, -1, -1};
    constexpr double cxd[8] = { 1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0};
    constexpr double cyd[8] = { 0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0};
    constexpr int opp[8]    = { 2,  3,  0,  1,  6,  7,  4,  5};

    constexpr double w0 = 4.0 / 9.0;
    constexpr double w[8] = {
        1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };

    template <typename Real>
    constexpr LBM_INLINE Real compute_feq(int d, Real rho, Real ux, Real uy, Real u_sq) noexcept {
        Real cu = cxd[d] * ux + cyd[d] * uy;
        return w[d] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sq);
    }

    template <typename Real>
    constexpr LBM_INLINE Real compute_feq0(Real rho, Real u_sq) noexcept {
        return w0 * rho * (1.0 - 1.5 * u_sq);
    }

    template <typename Real>
    constexpr LBM_INLINE Real compute_S_guo(int d, Real ux, Real uy, Real Fx, Real Fy, Real u_dot_F, double force_factor) noexcept {
        Real cu = cxd[d] * ux + cyd[d] * uy;
        Real cF = cxd[d] * Fx + cyd[d] * Fy;
        return (force_factor * w[d]) * (3.0 * (cF - u_dot_F) + 9.0 * cu * cF);
    }

    template <typename Real>
    constexpr LBM_INLINE Real compute_S0_guo(Real u_dot_F, double force_factor) noexcept {
        return (force_factor * w0) * (-3.0 * u_dot_F);
    }

    constexpr LBM_INLINE size_t src_idx(size_t idx, int nx, int d) noexcept {
        return idx - (cx[d] + cy[d] * nx);
    }

    constexpr LBM_INLINE size_t src_idx_p(int x, int y, int nx, int ny, int d) noexcept {
        return static_cast<size_t>((x - cx[d] + nx) % nx) + (static_cast<size_t>((y - cy[d] + ny) % ny) * nx);
    }

    template <typename Real>
    LBM_INLINE Real load_population(const double* __restrict__ ptr) noexcept {
        if constexpr (std::is_same_v<Real, double>) {
            return *ptr;
        } else {
            return _mm512_loadu_pd(ptr);
        }
    }

    template <bool UseNTStore, typename Real>
    LBM_INLINE void store_population(double* __restrict__ ptr, Real val) noexcept {
        if constexpr (std::is_same_v<Real, double>) {
            *ptr = val;
        } else if constexpr (!UseNTStore) {
            _mm512_store_pd(ptr, val);  // Escrita AVX-512 alinhada para a Cache L3
        } else {
            _mm512_stream_pd(ptr, val); // Escrita AVX-512 Non-Temporal direta para a RAM DDR5
        }
    }
}

#endif