#ifndef GRID3D_HPP
#define GRID3D_HPP

#include "types.hpp"
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <omp.h>

#ifdef __linux__
#include <sys/mman.h>
#endif

// ============================================================================
// CAMPO DE POPULAÇÕES GENÉRICO (Funciona para D3Q19, D3Q27, D3Q7 térmico/MHD!)
// Usa 1 Única Grelha com Esoteric Twist (f0 + NUM_PAIRS positivos + NUM_PAIRS opostos)
// ============================================================================
template <typename Lattice>
struct LatticeField {
    static constexpr size_t ALIGNMENT    = 64;
    static constexpr size_t NUM_PAIRS    = Lattice::NUM_PAIRS;
    static constexpr size_t TOTAL_PLANES = 1 + 2 * NUM_PAIRS; // 19 no D3Q19, 27 no D3Q27, 7 no D3Q7

    size_t plane_pitch{0}; // Em elementos lbm_real (SEMPRE size_t de 64 bits!)

    lbm_real* __restrict__ f_raw_ptr{nullptr};
    lbm_real* __restrict__ f0_ptr{nullptr};    // 1 plano (repouso)
    lbm_real* __restrict__ f_pos_ptr{nullptr}; // NUM_PAIRS planos (direções 0 .. NUM_PAIRS-1)
    lbm_real* __restrict__ f_neg_ptr{nullptr}; // NUM_PAIRS planos (direções opostas)

    LatticeField() = default;
    ~LatticeField() { deallocate(); }

    LatticeField(const LatticeField&) = delete;
    LatticeField& operator=(const LatticeField&) = delete;

    void allocate(size_t pitch_elements) {
        deallocate();
        plane_pitch = pitch_elements;

        // Nota: plane_pitch já é múltiplo de 16 -> total_bytes é sempre múltiplo de 64B (exigência do aligned_alloc)
        const size_t total_bytes = TOTAL_PLANES * plane_pitch * sizeof(lbm_real);
        f_raw_ptr = static_cast<lbm_real*>(std::aligned_alloc(ALIGNMENT, total_bytes));
        if (!f_raw_ptr) {
            throw std::bad_alloc();
        }
#ifdef __linux__
        madvise(f_raw_ptr, total_bytes, MADV_HUGEPAGE);
#endif
        f0_ptr    = f_raw_ptr;
        f_pos_ptr = f_raw_ptr + plane_pitch;
        f_neg_ptr = f_raw_ptr + (1 + NUM_PAIRS) * plane_pitch;
    }

    // Troca instantânea no fim de cada passo do Esoteric Twist (0 bytes movidos, elimina if (step % 2)!)
    inline void swap_halves() noexcept {
        std::swap(f_pos_ptr, f_neg_ptr);
    }

    void deallocate() noexcept {
        if (f_raw_ptr) {
            std::free(f_raw_ptr);
            f_raw_ptr = f0_ptr = f_pos_ptr = f_neg_ptr = nullptr;
        }
        plane_pitch = 0;
    }

    size_t memory_bytes() const noexcept {
        return f_raw_ptr ? (TOTAL_PLANES * plane_pitch * sizeof(lbm_real)) : 0;
    }
};

// ============================================================================
// GRELHA GEOMÉTRICA 3D + HALO ALINHADO + BLOCOS 32^3 + LISTAS SoA
// ============================================================================
template <typename FluidLattice>
struct Grid3D {
    static constexpr size_t ALIGNMENT = 64;
    static constexpr int    BLOCK_DIM = 32;

    // Halo: 16 em X dá alinhamento 64B em float (128B em double) com apenas +32 células/linha!
    // Em Y e Z basta 1 célula de cada lado.
    static constexpr int HALO_X = 16;
    static constexpr int HALO_Y = 1;
    static constexpr int HALO_Z = 1;

    // Dimensões físicas reais do fluido
    int nx{0}, ny{0}, nz{0};

    // Dimensões em memória (incluindo Halo e Padding de alinhamento)
    int stride_x{0}, stride_y{0}, stride_z{0};

    size_t phys_cells{0};   // nx * ny * nz
    size_t alloc_cells{0};  // stride_x * stride_y * stride_z
    size_t plane_pitch{0};  // Alinhado a 64B e protegido contra 4K Aliasing

    // Mapa de geometria (mesmo layout padded das populações -> 1 só índice idx para tudo!)
    CellType* __restrict__ flag_field_ptr{nullptr};

    // Campo de populações do fluido (D3Q19 ou D3Q27)
    // (No futuro: basta adicionar LatticeField<D3Q7> thermal; aqui!)
    LatticeField<FluidLattice> fluid;

    // Listas de blocos 32x32x32 pré-classificados (Zero if-else no interior!)
    std::vector<Block32> outer_blocks; // Blocos exteriores -> TRT / BGK rápido
    std::vector<Block32> inner_blocks; // Blocos perto do buraco negro -> MRT completo

    // Listas esparsas de fronteira em SoA puro ("Scaffold")
    WallLinkList     wall_links;   // Paredes retas (q=0.5) ou curvas (Bouzidi)
    OpenBoundaryList inlet_nodes;  // Entradas de velocidade/pressão
    OpenBoundaryList outlet_nodes; // Saídas / extrapolação
    SpongeList       sponge_nodes; // Zonas de amortecimento

    double vtk_io_time_sec{0.0};

    Grid3D() = default;
    ~Grid3D() { deallocate(); }

    Grid3D(const Grid3D&) = delete;
    Grid3D& operator=(const Grid3D&) = delete;

    // --- Conversores explícitos de coordenadas (evitam bugs off-by-one com o Halo!) ---

    // Coordenada física (0..nx-1, 0..ny-1, 0..nz-1) -> índice linear em memória
    inline lbm_int idx_phys(int px, int py, int pz) const noexcept {
        return static_cast<lbm_int>(
            (static_cast<size_t>(pz + HALO_Z) * stride_y + (py + HALO_Y)) * stride_x + (px + HALO_X)
        );
    }

    // Coordenada já em memória (HALO_X..HALO_X+nx-1, 1..ny, 1..nz) -> índice linear
    inline lbm_int idx_mem(int mx, int my, int mz) const noexcept {
        return static_cast<lbm_int>(
            (static_cast<size_t>(mz) * stride_y + my) * stride_x + mx
        );
    }

    void allocate_geometry(int x, int y, int z) {
        if (x <= 0 || y <= 0 || z <= 0) {
            throw std::invalid_argument("Dimensões 3D da grelha têm de ser > 0.");
        }
        deallocate();

        nx = x; ny = y; nz = z;
        phys_cells = static_cast<size_t>(nx) * ny * nz;

        // stride_x: [HALO_X à esquerda] + [nx reais] + [>=1 halo à direita], arredondado a múltiplo de 16
        stride_x = ((nx + HALO_X + 1 + 15) / 16) * 16;
        stride_y = ny + 2 * HALO_Y;
        stride_z = nz + 2 * HALO_Z;

        alloc_cells = static_cast<size_t>(stride_x) * stride_y * stride_z;

        // Verifica se o índice máximo cabe em lbm_int (uint32_t -> até ~4.29 mil milhões de células)
        if (alloc_cells >= static_cast<size_t>(UINT32_MAX)) {
            throw std::overflow_error("Grelha excede o limite de uint32_t (4.29G células). Altera lbm_int para uint64_t.");
        }

        // Alinha plane_pitch a múltiplo de 16 elementos (>= 64B em float e >= 128B em double)
        plane_pitch = ((alloc_cells + 15) / 16) * 16;

        // Proteção anti-4K Aliasing entre planos na Cache L1
        constexpr size_t page_elems = 4096 / sizeof(lbm_real);
        if ((plane_pitch % page_elems) == 0) {
            plane_pitch += 16;
        }

        // Alocar flag_field com HugePages e NUMA First-Touch paralelo por fatias Z
        const size_t flag_bytes = ((alloc_cells + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT;
        flag_field_ptr = static_cast<CellType*>(std::aligned_alloc(ALIGNMENT, flag_bytes));
        if (!flag_field_ptr) {
            deallocate();
            throw std::bad_alloc();
        }
#ifdef __linux__
        madvise(flag_field_ptr, flag_bytes, MADV_HUGEPAGE);
#endif
        #pragma omp parallel for schedule(static)
        for (int mz = 0; mz < stride_z; ++mz) {
            for (int my = 0; my < stride_y; ++my) {
                const size_t row = (static_cast<size_t>(mz) * stride_y + my) * stride_x;
                for (int mx = 0; mx < stride_x; ++mx) {
                    flag_field_ptr[row + mx] = CellType::FLUID;
                }
            }
        }
        // Limpa padding final para o Valgrind/ASan
        for (size_t i = alloc_cells; i < flag_bytes; ++i) {
            flag_field_ptr[i] = CellType::FLUID;
        }
    }

    void allocate_populations(size_t max_allowed_bytes = 0) {
        if (alloc_cells == 0) {
            throw std::logic_error("allocate_geometry() tem de ser chamado primeiro.");
        }
        const size_t needed_f_bytes = LatticeField<FluidLattice>::TOTAL_PLANES * plane_pitch * sizeof(lbm_real);
        const size_t projected = memory_usage_bytes() + needed_f_bytes;

        if (max_allowed_bytes > 0 && projected > max_allowed_bytes) {
            deallocate();
            throw std::runtime_error("Limite de memória excedido! Necessário: " +
                                     std::to_string(projected / (1024 * 1024)) + " MiB.");
        }
        fluid.allocate(plane_pitch);
    }

    void deallocate() noexcept {
        if (flag_field_ptr) {
            std::free(flag_field_ptr);
            flag_field_ptr = nullptr;
        }
        fluid.deallocate();
        outer_blocks.clear();
        inner_blocks.clear();
        wall_links.clear();
        inlet_nodes.clear();
        outlet_nodes.clear();
        sponge_nodes.clear();
        nx = ny = nz = stride_x = stride_y = stride_z = 0;
        phys_cells = alloc_cells = plane_pitch = 0;
    }

    size_t memory_usage_bytes() const noexcept {
        size_t bytes = 0;
        if (flag_field_ptr) bytes += ((alloc_cells + ALIGNMENT - 1) / ALIGNMENT) * ALIGNMENT;
        bytes += fluid.memory_bytes();
        bytes += (outer_blocks.capacity() + inner_blocks.capacity()) * sizeof(Block32);
        bytes += wall_links.memory_bytes();
        bytes += inlet_nodes.memory_bytes();
        bytes += outlet_nodes.memory_bytes();
        bytes += sponge_nodes.memory_bytes();
        return bytes;
    }
};

#endif // GRID3D_HPP