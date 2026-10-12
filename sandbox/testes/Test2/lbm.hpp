#ifndef GRID_H
#define GRID_H

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <filesystem>



constexpr size_t ALIGNMENT = 64; // Cache line / AVX-512

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
    uint8_t *flag_field_ptr{nullptr};

    // --- LISTAS ESPARSAS DE FRONTEIRA ---
    std::vector<int> solid_cells;
    std::vector<int> input_cells;
    std::vector<int> output_cells;

    // Construtor vazio
    Grid() = default;

    // Destrutor: Liberta a memória automaticamente (RAII)
    ~Grid() {
        deallocate();
    }

    // Impede a cópia acidental da estrutura (Rule of Three/Five)
    Grid(const Grid&) = delete;
    Grid& operator=(const Grid&) = delete;

    // Aloca os blocos densos assim que nx e ny forem conhecidos
    void allocate(int x, int y) {
        // Se já houver algo alocado, limpa primeiro
        deallocate();

        nx = x;
        ny = y;
        size = static_cast<size_t>(nx) * ny;

        // Arredonda para o múltiplo de 8 doubles superior (64 bytes)
        plane_pitch = ((size + 7) / 8) * 8;

        // Aloca 17 planos (1 repouso + 8 leitura + 8 escrita)
        size_t total_f_bytes = 17 * plane_pitch * sizeof(double);
        f_raw_ptr = static_cast<double*>(std::aligned_alloc(ALIGNMENT, total_f_bytes));

        if (!f_raw_ptr) {
            throw std::bad_alloc();
        }

        f0_ptr = f_raw_ptr;
        f1_ptr = f_raw_ptr + plane_pitch;
        f2_ptr = f_raw_ptr + (9 * plane_pitch);

        // Mapa de flags (alinhamento standard de 64 bytes)
        size_t flag_bytes = ((size + (ALIGNMENT - 1)) / ALIGNMENT) * ALIGNMENT;
        flag_field_ptr = static_cast<uint8_t*>(std::aligned_alloc(ALIGNMENT, flag_bytes));

        if (!flag_field_ptr) {
            std::free(f_raw_ptr);
            throw std::bad_alloc();
        }
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
        solid_cells.clear();
        input_cells.clear();
        output_cells.clear();
        nx = ny = 0;
        size = plane_pitch = 0;
    }
};



class LBM {

    private:
        Grid grid;

    public:
        LBM(int nx, int ny);
        ~LBM();

        void initialize();
        void stream_compute_collide();
        void apply_boundary();
        void export_to_vti(int step);

};
 
#endif