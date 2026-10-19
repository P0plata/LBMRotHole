#ifndef GRID_H
#define GRID_H

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <new>
#include <utility>
#include <stdexcept>
#include "D_Q_.hpp"
#include <cstring>



constexpr size_t ALIGNMENT = 64; // Cache line / AVX-512

enum class CellType : uint8_t { FLUID = 0, SOLID = 1, INLET = 2, OUTLET = 3, SPONGE = 4};

// Para Paredes e Curvas (Bounce-Back / Bouzidi): Ao nível do LINK (direção q)
// Tamanho exato: 16 bytes (cabem exatamente 4 links por cache line de 64B)
struct WallLink {
    int fluid_idx;   // Índice linear da célula de fluido (x_f)
    int neighbor_f;  // Índice linear da célula de fluido vizinha atrás (x_ff) para interpolação de Bouzidi
    uint8_t dir;     // Direção q (0..7 no array f1/f2) que sai do fluido e bate na parede
    uint8_t opp_dir; // Direção oposta q_opp (0..7) que reflete da parede para o fluido
    uint16_t _pad;   // Padding explícito para alinhamento de 4 bytes
    float q_frac;    // Fração de distância até à parede curva: q in [0, 1] (0.5 = half-way)
};

// Para Entradas/Saídas (Zou-He, Inflow/Outflow): Ao nível do NÓ
// Tamanho exato: 16 bytes
struct OpenBoundaryNode {
    int idx;         // Índice linear do nó de fronteira
    float vx;        // Velocidade imposta X (ou densidade/pressão alvo)
    float vy;        // Velocidade imposta Y
    uint8_t normal;  // 0 = Este (+x), 1 = Oeste (-x), 2 = Norte (+y), 3 = Sul (-y)
    uint8_t _pad[3]; // Padding explícito para fechar 16 bytes
};

// Para Camadas de Esponja (Sponge Layer / Absorbing Zone): Ao nível do NÓ
// Tamanho exato: 20 bytes (cabem exatamente 2 por cache line de 64B)
struct SpongeNode {
    int idx;         // 4 bytes
    float sigma;     // 4 bytes
    float rho_ref;   // 4 bytes
    float vx_ref;    // 4 bytes
    float vy_ref;    // 4 bytes
};

static_assert(sizeof(WallLink) == 16, "WallLink tem de ter exatamente 16 bytes!");
static_assert(sizeof(OpenBoundaryNode) == 16, "OpenBoundaryNode tem de ter 16 bytes!");
static_assert(sizeof(SpongeNode) == 20, "SpongeNode tem de ter exatamente 32 bytes!");

struct Grid {
    int nx{0};
    int ny{0};
    size_t size{0};        // nx * ny
    size_t plane_pitch{0}; // Múltiplo de 8 doubles para garantir 64 bytes por plano

    // --- MACRO-BLOCO DE POPULAÇÕES (17 * plane_pitch) ---
    double *f_raw_ptr{nullptr};
    double *f0_ptr{nullptr}; // Direção 0 (repouso)
    double *f1_ptr{nullptr}; // 8 direções móveis - Leitura (Pull)
    double *f2_ptr{nullptr}; // 8 direções móveis - Escrita

    // --- MAPA DENSO DA GEOMETRIA (1 byte/célula) ---
    CellType *flag_field_ptr{nullptr};


    // --- LISTAS ESPARSAS ESPECIALIZADAS ---
    std::vector<WallLink> wall_links;            // Paredes retas (q=0.5) ou curvas (Bouzidi)
    std::vector<OpenBoundaryNode> inlet_nodes;   // Entradas de velocidade/pressão
    std::vector<OpenBoundaryNode> outlet_nodes;  // Saídas / extrapolação
    std::vector<SpongeNode> sponge_nodes;        // Zonas de amortecimento


    // Construtor vazio
    Grid() = default;

    // Destrutor: Liberta a memória automaticamente (RAII)
    ~Grid() {
        deallocate();
    }

    // Impede a cópia acidental da estrutura (Rule of Three/Five)
    Grid(const Grid&) = delete;
    Grid& operator=(const Grid&) = delete;

    void allocate_geometry(int x, int y){
        
        
        if (x <= 0 || y <= 0) {
            throw std::invalid_argument("As dimensões da grelha têm de ser estritamente positivas (> 0).");
        }

        // Se já houver algo alocado, limpa primeiro
        deallocate();

        nx = x;
        ny = y;
        size = static_cast<size_t>(nx) * ny;
        // Arredonda para o múltiplo de 8 doubles superior (64 bytes)
        plane_pitch = ((size + 7) / 8) * 8;
        
        // Mapa de flags (alinhamento standard de 64 bytes)
        size_t flag_bytes = ((size + (ALIGNMENT - 1)) / ALIGNMENT) * ALIGNMENT;
        flag_field_ptr = static_cast<CellType*>(std::aligned_alloc(ALIGNMENT, flag_bytes));

        if (!flag_field_ptr) {
            deallocate();
            throw std::bad_alloc();
        }

    }

        // Aloca os blocos densos assim que nx e ny forem conhecidos
    void allocate_populations(size_t max_allowed_bytes = 0) {

        if (size == 0) throw std::logic_error("allocate_geometry() needed to be called first.");

        // Aloca 17 planos (1 repouso + 8 leitura + 8 escrita)
        size_t total_f_bytes = 17 * plane_pitch * sizeof(double);
        size_t projected_total = memory_usage_bytes() + total_f_bytes;

        if (max_allowed_bytes > 0 && projected_total > max_allowed_bytes) {
            deallocate();
            throw std::runtime_error("Limite de memória excedido! Necessário: " + 
                                     std::to_string(projected_total / (1024 * 1024)) + " MB.");
        }

        f_raw_ptr = static_cast<double*>(std::aligned_alloc(ALIGNMENT, total_f_bytes));

        if (!f_raw_ptr) {
            deallocate();
            throw std::bad_alloc();
        }

        f0_ptr = f_raw_ptr;
        f1_ptr = f_raw_ptr + plane_pitch;
        f2_ptr = f_raw_ptr + (9 * plane_pitch);
    }


    void deallocate() {
        if (f_raw_ptr) {
            std::free(f_raw_ptr);
            f_raw_ptr = f0_ptr = f1_ptr = f2_ptr = nullptr;
        }
        if (flag_field_ptr) {
            std::free(flag_field_ptr);
            flag_field_ptr = nullptr;
        }

        wall_links.clear();   wall_links.shrink_to_fit();
        inlet_nodes.clear();  inlet_nodes.shrink_to_fit();
        outlet_nodes.clear(); outlet_nodes.shrink_to_fit();
        sponge_nodes.clear(); sponge_nodes.shrink_to_fit();

        nx = ny = 0;
        size = plane_pitch = 0;
    }

    void swap_populations() noexcept { 
        std::swap(f1_ptr, f2_ptr); 
    }

    size_t memory_usage_bytes() const noexcept {
        size_t flag_bytes = ((size + (ALIGNMENT - 1)) / ALIGNMENT) * ALIGNMENT;
        size_t pop_bytes  = f_raw_ptr ? (17 * plane_pitch * sizeof(double)) : 0;
        size_t bnd_bytes  = wall_links.capacity()   * sizeof(WallLink) +
                            inlet_nodes.capacity()  * sizeof(OpenBoundaryNode) +
                            outlet_nodes.capacity() * sizeof(OpenBoundaryNode) +
                            sponge_nodes.capacity() * sizeof(SpongeNode);
        return flag_bytes + pop_bytes + bnd_bytes;
        }

};

#endif // GRID_H