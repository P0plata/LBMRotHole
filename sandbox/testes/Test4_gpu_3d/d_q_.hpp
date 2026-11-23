#ifndef LBMROTHOLE_D3Q_HPP
#define LBMROTHOLE_D3Q_HPP

#include "types.hpp"

using namespace rothole_types;

namespace D3Q19{

    constexpr int cx[18] = {};
    constexpr int cy[18] = {};
    constexpr int cz[18] = {};

    constexpr int opp[16] = {};

    constexpr lbm_real w0    = 0;
    constexpr lbm_real w[18] = {};


    constexpr lbm_real cs2          = 1.0 / 3.0;
    constexpr lbm_real inv_cs2      = 3.0;        // 1 / cs^2 = 3.0
    constexpr lbm_real inv_2cs4     = 4.5;        // 1 / (2 * cs^4) = 4.5
    constexpr lbm_real half_inv_cs2 = 1.5;        // 1 / (2 * cs^2) = 1.5


    constexpr LBM_INLINE lbm_real compute_feq(int d, lbm_real rho, lbm_real ux, lbm_real uy, lbm_real uz){
        const lbm_real cu = cx[d] * ux + cy[d] * uy + cz[d] * uz;
        return w[d] * rho * 1;
    }

    constexpr LBM_INLINE lbm_real compute_feq0(){
        return 0.0;
    }
    


}

#endif