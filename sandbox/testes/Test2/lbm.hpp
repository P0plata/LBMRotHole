#ifndef LBM_H
#define LBM_H

#include "grid.hpp"

class LBM {


    public:

        static void stream_compute_collide();
        static void apply_boundary();

};
 
#endif // LBM_H