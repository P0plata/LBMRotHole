
#ifndef LBMROTHOLE_TYPES
#define LBMROTHOLE_TYPES

#include <cstdint>
#include <vector>

#define LBM_INLINE __attribute__((always_inline)) static inline

// Trocar entre float (FP32 - GPU/3D grande) e double (FP64 - verificação de precisão)
using lbm_real = float;
using lbm_int  = uint32_t;

enum class CellType : uint8_t {
    FLUID  = 0,
    SOLID  = 1,
    INLET  = 2,
    OUTLET = 3,
    SPONGE = 4
};

enum class CollisionType : uint8_t {
    BGK = 0,
    TRT = 1,
    MRT = 2
};

enum class BoundaryCondition : uint8_t {
    PERIODIC  = 0,
    DIRICHLET = 1,
    NEUMANN   = 2
};

// ============================================================================
// BLOCO 32x32x32 (Metadados para Cache Blocking em CPU e Warps em GPU)
// ============================================================================
enum class BlockZone : uint8_t {
    EMPTY_SOLID = 0, // 100% sólido -> ignorado pelo solver (0 FLOPs!)
    OUTER_FAST  = 1, // Longe do centro -> corre TRT/BGK rápido
    INNER_HIGH  = 2  // Perto do buraco negro -> corre MRT completo
};

// Tamanho: 16 bytes
struct Block32 {
    uint16_t x0; // Coordenada inicial EM MEMÓRIA (já inclui HALO_X, HALO_Y, HALO_Z)
    uint16_t y0;
    uint16_t z0; 
    uint16_t sx; // Dimensões efetivas do bloco (<= 32, suporta malhas pequenas de validação!)
    uint16_t sy;
    uint16_t sz;
    BlockZone zone;      // Tipo de colisor a despachar
    uint8_t   _pad[3];
};

// ============================================================================
// LISTAS DE FRONTEIRA EM SoA PURO (Zero Padding, 100% Coalescido para GPU)
// ============================================================================

// 1. Paredes / Bouzidi (14 Bytes/link)
struct WallLinkList {
    std::vector<lbm_int> fluid_idx;   // Célula fluida x_f (índice no array com halo)
    std::vector<lbm_int> neighbor_f;  // Célula fluida atrás x_ff (para Bouzidi q < 0.5)
    std::vector<float>   q_frac;      // Distância fracionária à parede curva (0.5 = half-way)
    std::vector<uint8_t> dir;         // Direção que sai do fluido para a parede
    std::vector<uint8_t> opp_dir;     // Direção oposta que regressa ao fluido

    void push_back(lbm_int xf, lbm_int xff, uint8_t d, uint8_t d_opp, float q) {
        fluid_idx.push_back(xf);
        neighbor_f.push_back(xff);
        q_frac.push_back(q);
        dir.push_back(d);
        opp_dir.push_back(d_opp);
    }

    size_t size() const noexcept { return fluid_idx.size(); }
    void clear() { fluid_idx.clear(); neighbor_f.clear(); q_frac.clear(); dir.clear(); opp_dir.clear(); }
    void shrink_to_fit() {
        fluid_idx.shrink_to_fit(); neighbor_f.shrink_to_fit(); q_frac.shrink_to_fit();
        dir.shrink_to_fit(); opp_dir.shrink_to_fit();
    }
    size_t memory_bytes() const noexcept {
        return (fluid_idx.capacity() + neighbor_f.capacity()) * sizeof(lbm_int)
             + q_frac.capacity() * sizeof(float)
             + (dir.capacity() + opp_dir.capacity()) * sizeof(uint8_t);
    }
};

// 2. Fronteiras Abertas 3D (Inlet / Outlet em SoA — 23 Bytes/nó)
struct OpenBoundaryList {
    std::vector<lbm_int>  idx;
    std::vector<lbm_real> rho;  // Densidade alvo (Outlet)
    std::vector<lbm_real> vx;   // Velocidade alvo (Inlet)
    std::vector<lbm_real> vy;
    std::vector<lbm_real> vz;
    std::vector<int8_t>   nx;   // Normal COM SINAL (-1, 0, +1)!
    std::vector<int8_t>   ny;
    std::vector<int8_t>   nz; 

    void push_back(lbm_int i, lbm_real r, lbm_real x, lbm_real y, lbm_real z,
                   int8_t norm_x, int8_t norm_y, int8_t norm_z) {
        idx.push_back(i); rho.push_back(r);
        vx.push_back(x); vy.push_back(y); vz.push_back(z);
        nx.push_back(norm_x); ny.push_back(norm_y); nz.push_back(norm_z);
    }

    size_t size() const noexcept { return idx.size(); }
    void clear() { idx.clear(); rho.clear(); vx.clear(); vy.clear(); vz.clear(); nx.clear(); ny.clear(); nz.clear(); }
    void shrink_to_fit() {
        idx.shrink_to_fit(); rho.shrink_to_fit(); vx.shrink_to_fit(); vy.shrink_to_fit(); vz.shrink_to_fit();
        nx.shrink_to_fit(); ny.shrink_to_fit(); nz.shrink_to_fit();
    }
    size_t memory_bytes() const noexcept {
        return idx.capacity() * sizeof(lbm_int)
             + (rho.capacity() + vx.capacity() + vy.capacity() + vz.capacity()) * sizeof(lbm_real)
             + (nx.capacity() + ny.capacity() + nz.capacity()) * sizeof(int8_t);
    }
};

// 3. Camadas Esponja em SoA com Tabela de Alvos (Apenas 9 Bytes/nó!)
struct SpongeTarget {
    float   rho_ref{1.0f};
    float   vx_ref{0.0f};
    float   vy_ref{0.0f};
    float   vz_ref{0.0f};
    uint8_t is_sink{0}; // 1 = Sumidouro central (entra no cálculo de dM/dt de Bondi!), 0 = Parede exterior
};

struct SpongeList {
    std::vector<SpongeTarget> targets; // Pequena tabela de alvos (ex: buraco central + 6 paredes)
    std::vector<lbm_int>      idx;     // 4 bytes
    std::vector<float>        sigma;   // 4 bytes
    std::vector<uint8_t>      zone_id; // 1 byte -> índice em targets[]

    void push_back(lbm_int i, float sig, uint8_t zid) {
        idx.push_back(i); sigma.push_back(sig); zone_id.push_back(zid);
    }

    size_t size() const noexcept { return idx.size(); }
    void clear() { targets.clear(); idx.clear(); sigma.clear(); zone_id.clear(); }
    void shrink_to_fit() {
        targets.shrink_to_fit(); idx.shrink_to_fit(); sigma.shrink_to_fit(); zone_id.shrink_to_fit();
    }
    size_t memory_bytes() const noexcept {
        return targets.capacity() * sizeof(SpongeTarget)
             + idx.capacity() * sizeof(lbm_int)
             + sigma.capacity() * sizeof(float)
             + zone_id.capacity() * sizeof(uint8_t);
    }
};


#endif 