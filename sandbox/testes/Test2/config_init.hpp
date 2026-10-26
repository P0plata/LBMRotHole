#ifndef CONFIG_INIT_HPP
#define CONFIG_INIT_HPP

#include "grid.hpp"
#include "utils.hpp"
#include <string>
#include <vector>

enum class TypeGen : uint8_t { FILE = 0, PROCEDURAL = 1 };
enum class FormatType : uint8_t { RAW = 0, PGM = 1, CIRCLE = 2, RECTANGLE = 3 };
enum class CollisionType : uint8_t {  BGK = 0, SRT = 0, TRT = 1, MRT = 2};
enum class BoundaryCondition : uint8_t { PERIODIC = 0, DIRICHLET = 1, NEUMANN = 2  }; // periodic, 1: Dirichlet (Inlet/Wall), 2: Neumann (Outlet)
enum class BoundarySolver : uint8_t { BOUNCEBACK = 0, ZOUHE = 1 }; // 0: bounce-back, 1: Zou-He

struct GeometryConfig {
    TypeGen type{TypeGen::PROCEDURAL};
    FormatType format{FormatType::CIRCLE};
    CellType typeboundary{CellType::SOLID}; // SOLID, INLET, OUTLET, SPONGE
    std::string filename{""};

    // Dimensões (para ficheiro ou Rectangle)
    int obj_nx{0}, obj_ny{0}, obj_nz{1};

    // Posição na grelha (para Circle: centro X, Y; ex: meio-círculo na base -> offset_x = nx/2, offset_y = 0)
    int offset_x{0}, offset_y{0}, offset_z{0};

    // Raio (para Circle)
    double radius{0.0};

    // Parâmetros de Esponja e Condições no Objeto
    double sponge_width{0.0};     // Se > 0 num objeto SOLID, cria uma camada de esponja ao redor (de R a R + sponge_width)
    float sponge_sigma_max{0.0f}; // Força máxima de atenuação da esponja
    float vx{0.0f};               // Velocidade/valor alvo se o objeto for INLET/OUTLET/SPONGE
    float vy{0.0f};
    float rho{1.0f};
};

struct Simconfig {
    // Simulation
    bool use_custom_init{false};
    int nx{0};
    int ny{0};
    int nz{1};
    size_t grid_size{0}; // nx * ny * nz
    int grid_type{0};    // 0: 2D, 1: 3D
    int max_steps{0};
    int export_freq{100};
    std::string output_prefix{""};
    size_t max_allowed_bytes{0};

    // Estado inicial do fluido
    double rho0{1.0};
    double ux0{0.0};
    double uy0{0.0};

    // Geometry
    std::vector<GeometryConfig> obstacles;

    // Collision Model
    CollisionType type_collision{CollisionType::BGK};            // 0: BGK, 1: TRT, 2: MRT
    double tau_BGK{-1.0};             // Relaxation time BGK
    double tau_TRT{-1.0};             // Relaxation time TRT (par/ímpar)
    std::string tau_MRT_filename{""}; // Ficheiro com os tempos de relaxamento MRT
    double s_MRT[9]{0.0};             // Coeficientes MRT carregados em memória

    // Boundary Conditions nas 6 faces: x0 (Oeste), x1 (Este), y0 (Sul), y1 (Norte), z0, z1
    BoundaryCondition boundary_conditions[6]{BoundaryCondition::PERIODIC};  
    BoundarySolver type_boundary_solver[6]{BoundarySolver::BOUNCEBACK}; // 0: bounce-back, 1: Zou-He
    CellType face_cell_type[6]{CellType::FLUID, CellType::FLUID, CellType::FLUID,
                               CellType::FLUID, CellType::FLUID, CellType::FLUID};
    float face_vx[6]{0.0f};
    float face_vy[6]{0.0f};
    float face_rho[6]{1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
    int face_sponge_width[6]{0};
    float face_sponge_sigma[6]{0.0f};
};

namespace GeometryBuilder {
    void stamp_geometry(Grid& grid, const Simconfig& config);
    void build_boundary_lists(Grid& grid, const Simconfig& config);
}

class Initializer {
public:
    static void load_config_from_json(const std::string& json_filepath, Simconfig& config);
    static void initialize(Grid& grid, const Simconfig& config);
    static void apply_custom_initialization(Grid& grid, const Simconfig& config);
};

#endif // CONFIG_INIT_HPP