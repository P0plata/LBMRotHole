#include "config_init.hpp"
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

namespace GeometryBuilder {

// Uma célula faz parte do domínio fluido se for FLUID ou SPONGE
static inline bool is_fluid_domain(CellType c) {
    return c == CellType::FLUID || c == CellType::SPONGE;
}

// Calcula a interseção exata raio-círculo q in (0, 1] para o método de Bouzidi
static float compute_circle_q_frac(int xf, int yf, int cx_dir, int cy_dir, const GeometryConfig& circle) {
    const double dx = static_cast<double>(xf - circle.offset_x);
    const double dy = static_cast<double>(yf - circle.offset_y);
    const double r2 = circle.radius * circle.radius;

    const double A = static_cast<double>(cx_dir * cx_dir + cy_dir * cy_dir);
    const double B = dx * cx_dir + dy * cy_dir;
    const double C = dx * dx + dy * dy - r2;

    const double disc = B * B - A * C;
    if (disc < 0.0) return 0.5f;

    // Raiz na direção do fluido (fora, C > 0) para o sólido (dentro, C <= 0)
    double q = (-B - std::sqrt(disc)) / A;
    if (q <= 0.0 || q > 1.0) {
        q = (-B + std::sqrt(disc)) / A;
    }
    if (q > 0.0 && q <= 1.0) {
        return static_cast<float>(q);
    }
    return 0.5f;
}



static void load_raw_file(Grid& grid, const GeometryConfig& obs) {
    std::ifstream file(obs.filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o ficheiro RAW: " + obs.filename);
    }

    std::vector<uint8_t> buffer(static_cast<size_t>(obs.obj_nx) * obs.obj_ny);
    file.read(reinterpret_cast<char*>(buffer.data()), buffer.size());

    for (int y = 0; y < obs.obj_ny; ++y) {
        int gy = obs.offset_y + y;
        if (gy < 0 || gy >= grid.ny) continue;
        for (int x = 0; x < obs.obj_nx; ++x) {
            int gx = obs.offset_x + x;
            if (gx < 0 || gx >= grid.nx) continue;
            if (buffer[y * obs.obj_nx + x] != 0) {
                grid.flag_field_ptr[gy * grid.nx + gx] = obs.typeboundary;
            }
        }
    }
}

static void load_pgm_file(Grid& grid, const GeometryConfig& obs) {
    std::ifstream file(obs.filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Não foi possível abrir o ficheiro PGM: " + obs.filename);
    }

    std::string magic;
    int width = 0, height = 0, max_val = 0;
    file >> magic;

    char c;
    while (file.get(c)) {
        if (c == '#') {
            std::string dummy;
            std::getline(file, dummy);
        } else if (!std::isspace(static_cast<unsigned char>(c))) {
            file.putback(c);
            break;
        }
    }

    file >> width >> height >> max_val;
    file.get();

    for (int y = 0; y < height; ++y) {
        int gy = obs.offset_y + (height - 1 - y);
        for (int x = 0; x < width; ++x) {
            int gx = obs.offset_x + x;
            uint8_t pixel = (magic == "P5") ? static_cast<uint8_t>(file.get()) : [&]() {
                int v; file >> v; return static_cast<uint8_t>(v);
            }();

            if (gx >= 0 && gx < grid.nx && gy >= 0 && gy < grid.ny && pixel < 128) {
                grid.flag_field_ptr[gy * grid.nx + gx] = obs.typeboundary;
            }
        }
    }
}

void stamp_geometry(Grid& grid, const Simconfig& config) {
    const int nx = grid.nx;
    const int ny = grid.ny;

    // 1. Marcar as 4 faces do domínio 2D se não forem periódicas/fluido
    if (config.face_cell_type[0] != CellType::FLUID) {
        for (int y = 0; y < ny; ++y) grid.flag_field_ptr[y * nx + 0] = config.face_cell_type[0];
    }
    if (config.face_cell_type[1] != CellType::FLUID) {
        for (int y = 0; y < ny; ++y) grid.flag_field_ptr[y * nx + (nx - 1)] = config.face_cell_type[1];
    }
    if (config.face_cell_type[2] != CellType::FLUID) {
        for (int x = 0; x < nx; ++x) grid.flag_field_ptr[0 * nx + x] = config.face_cell_type[2];
    }
    if (config.face_cell_type[3] != CellType::FLUID) {
        for (int x = 0; x < nx; ++x) grid.flag_field_ptr[(ny - 1) * nx + x] = config.face_cell_type[3];
    }

    // 2. Desenhar os objetos (incluindo meios-círculos, retângulos, ficheiros PGM/RAW)
    for (const auto& obs : config.obstacles) {
        if (obs.type == TypeGen::PROCEDURAL) {
            if (obs.format == FormatType::CIRCLE) {
                const double r2 = obs.radius * obs.radius;
                const int min_x = std::max(0, static_cast<int>(std::floor(obs.offset_x - obs.radius)));
                const int max_x = std::min(nx - 1, static_cast<int>(std::ceil(obs.offset_x + obs.radius)));
                const int min_y = std::max(0, static_cast<int>(std::floor(obs.offset_y - obs.radius)));
                const int max_y = std::min(ny - 1, static_cast<int>(std::ceil(obs.offset_y + obs.radius)));

                for (int y = min_y; y <= max_y; ++y) {
                    for (int x = min_x; x <= max_x; ++x) {
                        double dx = x - obs.offset_x;
                        double dy = y - obs.offset_y;
                        if (dx * dx + dy * dy <= r2) {
                            int idx = y * nx + x;
                            // Uma esponja não deve apagar uma parede SOLID já existente
                            if (obs.typeboundary == CellType::SPONGE && grid.flag_field_ptr[idx] != CellType::FLUID) {
                                continue;
                            }
                            grid.flag_field_ptr[idx] = obs.typeboundary;
                        }
                    }
                }
            } else if (obs.format == FormatType::RECTANGLE) {
                const int min_x = std::max(0, obs.offset_x);
                const int max_x = std::min(nx - 1, obs.offset_x + obs.obj_nx - 1);
                const int min_y = std::max(0, obs.offset_y);
                const int max_y = std::min(ny - 1, obs.offset_y + obs.obj_ny - 1);

                for (int y = min_y; y <= max_y; ++y) {
                    for (int x = min_x; x <= max_x; ++x) {
                        int idx = y * nx + x;
                        if (obs.typeboundary == CellType::SPONGE && grid.flag_field_ptr[idx] == CellType::SOLID) {
                            continue;
                        }
                        grid.flag_field_ptr[idx] = obs.typeboundary;
                    }
                }
            }
        } else if (obs.type == TypeGen::FILE) {
            if (obs.format == FormatType::RAW) load_raw_file(grid, obs);
            else if (obs.format == FormatType::PGM) load_pgm_file(grid, obs);
        }
    }
}

void build_boundary_lists(Grid& grid, const Simconfig& config) {
    const int nx = grid.nx;
    const int ny = grid.ny;
    const bool periodic_x = (config.boundary_conditions[0] == BoundaryCondition::PERIODIC && config.boundary_conditions[1] == BoundaryCondition::PERIODIC);
    const bool periodic_y = (config.boundary_conditions[2] == BoundaryCondition::PERIODIC && config.boundary_conditions[3] == BoundaryCondition::PERIODIC);

    // Varredura única da grelha: Garante que todas as listas (especialmente esponjas) ficam ordenadas
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const int idx = y * nx + x;
            CellType cell = grid.flag_field_ptr[idx];

            // 1. Zonas de Esponja (SPONGE)
            if (is_fluid_domain(cell)) {
                float max_sigma = 0.0f;
                float t_rho = config.rho0;
                float t_vx  = config.ux0;
                float t_vy  = config.uy0;

                // 1A. Esponjas nas 4 faces retilíneas
                if (config.face_sponge_width[0] > 0 && x < config.face_sponge_width[0]) {
                    float xi = static_cast<float>(config.face_sponge_width[0] - x) / config.face_sponge_width[0];
                    float s = config.face_sponge_sigma[0] * xi * xi;
                    if (s > max_sigma) { max_sigma = s; t_rho = config.face_rho[0]; t_vx = config.face_vx[0]; t_vy = config.face_vy[0]; }
                }
                if (config.face_sponge_width[1] > 0 && x >= (nx - config.face_sponge_width[1])) {
                    float xi = static_cast<float>(x - (nx - 1 - config.face_sponge_width[1])) / config.face_sponge_width[1];
                    float s = config.face_sponge_sigma[1] * xi * xi;
                    if (s > max_sigma) { max_sigma = s; t_rho = config.face_rho[1]; t_vx = config.face_vx[1]; t_vy = config.face_vy[1]; }
                }
                if (config.face_sponge_width[2] > 0 && y < config.face_sponge_width[2]) {
                    float xi = static_cast<float>(config.face_sponge_width[2] - y) / config.face_sponge_width[2];
                    float s = config.face_sponge_sigma[2] * xi * xi;
                    if (s > max_sigma) { max_sigma = s; t_rho = config.face_rho[2]; t_vx = config.face_vx[2]; t_vy = config.face_vy[2]; }
                }
                if (config.face_sponge_width[3] > 0 && y >= (ny - config.face_sponge_width[3])) {
                    float xi = static_cast<float>(y - (ny - 1 - config.face_sponge_width[3])) / config.face_sponge_width[3];
                    float s = config.face_sponge_sigma[3] * xi * xi;
                    if (s > max_sigma) { max_sigma = s; t_rho = config.face_rho[3]; t_vx = config.face_vx[3]; t_vy = config.face_vy[3]; }
                }

                // 1B. Esponjas associadas a Objetos
                for (const auto& obs : config.obstacles) {
                    if (obs.format == FormatType::CIRCLE) {
                        const double dx = x - obs.offset_x;
                        const double dy = y - obs.offset_y;
                        const double r2 = dx * dx + dy * dy;

                        // Caso 1: Disco central absorvente (SPONGE com sponge_width == 0)
                        // Força = 0 em r = R, cresce para dentro até força máxima em r = 0
                        if (obs.typeboundary == CellType::SPONGE && obs.sponge_width <= 0.0 && obs.radius > 0.0) {
                            if (r2 <= obs.radius * obs.radius) {
                                double r = std::sqrt(r2);
                                float xi = static_cast<float>((obs.radius - r) / obs.radius);
                                float s  = obs.sponge_sigma_max * xi * xi;
                                if (s > max_sigma) { max_sigma = s; t_rho = obs.rho; t_vx = obs.vx; t_vy = obs.vy; }
                            }
                        }
                        // Caso 2: Anel absorvente EXTERIOR (SPONGE com sponge_width > 0)
                        // Força = 0 em r = R, cresce para fora até força máxima em r >= R + sponge_width
                        else if (obs.typeboundary == CellType::SPONGE && obs.sponge_width > 0.0) {
                            if (r2 >= obs.radius * obs.radius) {
                                double r = std::sqrt(r2);
                                float xi = static_cast<float>( std::min(1.0, (r - obs.radius) / obs.sponge_width) );
                                float s  = obs.sponge_sigma_max * xi * xi;
                                if (s > max_sigma) { max_sigma = s; t_rho = obs.rho; t_vx = obs.vx; t_vy = obs.vy; }
                            }
                        }
                        // Caso 3: Camada de esponja em redor de um objeto SOLID interno (r entre R e R + sponge_width)
                        // Força máxima junto à parede (r = R) e decai até 0 em r = R + sponge_width
                        else if (obs.typeboundary == CellType::SOLID && obs.sponge_width > 0.0 && obs.sponge_sigma_max > 0.0f) {
                            const double r_outer = obs.radius + obs.sponge_width;
                            if (r2 > obs.radius * obs.radius && r2 <= r_outer * r_outer) {
                                double r = std::sqrt(r2);
                                float xi = static_cast<float>((r_outer - r) / obs.sponge_width);
                                float s  = obs.sponge_sigma_max * xi * xi;
                                if (s > max_sigma) { max_sigma = s; t_rho = obs.rho; t_vx = obs.vx; t_vy = obs.vy; }
                            }
                        }
                    }
                }

                // Se detetou influência de esponja, converte o nó e regista-o sequencialmente
                if (max_sigma > 0.0f) {
                    grid.flag_field_ptr[idx] = CellType::SPONGE;
                    cell = CellType::SPONGE; // Atualiza para as fases seguintes
                    grid.sponge_nodes.push_back({idx, max_sigma, t_rho, t_vx, t_vy});
                }
            }

            // 2. WallLinks (Bouzidi / Bounce-Back)
            if (is_fluid_domain(cell)) {
                for (uint8_t d = 0; d < 8; ++d) {
                    int nx_nb = x + D2Q9::cx[d];
                    int ny_nb = y + D2Q9::cy[d];

                    if (periodic_x) nx_nb = (nx_nb + nx) % nx;
                    if (periodic_y) ny_nb = (ny_nb + ny) % ny;

                    if (nx_nb >= 0 && nx_nb < nx && ny_nb >= 0 && ny_nb < ny) {
                        int nb_idx = ny_nb * nx + nx_nb;
                        
                        if (grid.flag_field_ptr[nb_idx] == CellType::SOLID) {
                            int x_ff = x - D2Q9::cx[d];
                            int y_ff = y - D2Q9::cy[d];
                            if (periodic_x) x_ff = (x_ff + nx) % nx;
                            if (periodic_y) y_ff = (y_ff + ny) % ny;

                            // CORREÇÃO: Verificar se a célula de trás (x_ff, y_ff) é efetivamente FLUID
                            bool valid_ff = (x_ff >= 0 && x_ff < nx && y_ff >= 0 && y_ff < ny) &&
                                            is_fluid_domain(grid.flag_field_ptr[y_ff * nx + x_ff]);
                                            
                            int ff_idx = valid_ff ? (y_ff * nx + x_ff) : idx;

                            float q_val = 0.5f; // Bounce-back padrão (canto cego)
                            
                            // Só tentamos Bouzidi se a célula de trás for fluido válido
                            if (valid_ff) {
                                for (const auto& obs : config.obstacles) {
                                    if (obs.type == TypeGen::PROCEDURAL && obs.format == FormatType::CIRCLE &&
                                        obs.typeboundary == CellType::SOLID) {
                                        double dist_nb = std::hypot(nx_nb - obs.offset_x, ny_nb - obs.offset_y);
                                        if (dist_nb <= obs.radius) {
                                            q_val = compute_circle_q_frac(x, y, D2Q9::cx[d], D2Q9::cy[d], obs);
                                            break;
                                        }
                                    }
                                }
                            }

                            grid.wall_links.push_back({idx, ff_idx, d, static_cast<uint8_t>(D2Q9::opp[d]), 0, q_val});
                        }
                    }
                }
            }
            // 3. INLET / OUTLET (Mantém-se igual)
            else if (cell == CellType::INLET || cell == CellType::OUTLET) {
                uint8_t normal = 0;
                int face_id = 0;

                if (x + 1 < nx && is_fluid_domain(grid.flag_field_ptr[y * nx + (x + 1)])) {
                    normal = 0; face_id = 0; // Face Oeste (x0): normal aponta para Este (+x)
                } else if (x - 1 >= 0 && is_fluid_domain(grid.flag_field_ptr[y * nx + (x - 1)])) {
                    normal = 1; face_id = 1; // Face Este (x1): normal aponta para Oeste (-x)
                } else if (y + 1 < ny && is_fluid_domain(grid.flag_field_ptr[(y + 1) * nx + x])) {
                    normal = 2; face_id = 2; // Face Sul (y0): normal aponta para Norte (+y)
                } else if (y - 1 >= 0 && is_fluid_domain(grid.flag_field_ptr[(y - 1) * nx + x])) {
                    normal = 3; face_id = 3; // Face Norte (y1): normal aponta para Sul (-y)
                } else {
                    continue;
                }

                OpenBoundaryNode node{};
                node.idx    = idx;
                node.normal = normal;
                node.vx     = (cell == CellType::INLET) ? config.face_vx[face_id] : config.face_rho[face_id];
                node.vy     = config.face_vy[face_id];

                if (cell == CellType::INLET) grid.inlet_nodes.push_back(node);
                else grid.outlet_nodes.push_back(node);
            }
        }
    }

    grid.wall_links.shrink_to_fit();
    grid.inlet_nodes.shrink_to_fit();
    grid.outlet_nodes.shrink_to_fit();
    grid.sponge_nodes.shrink_to_fit();
}

} // namespace GeometryBuilder