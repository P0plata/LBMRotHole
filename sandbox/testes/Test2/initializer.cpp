#include "config_init.hpp"
#include "D_Q_.hpp"
#include <fstream>
#include <stdexcept>
#include <string>
#include <cctype>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <omp.h>

using json = nlohmann::json;

static size_t parse_memory_bytes(const nlohmann::json& val) {
    // Se no JSON for passado diretamente um número inteiro em bytes: "max_allowed_bytes": 85899345920
    if (val.is_number_integer() || val.is_number_unsigned()) {
        return val.get<size_t>();
    }

    if (!val.is_string()) {
        throw std::runtime_error("max_allowed_bytes deve ser uma string (ex: '80GB') ou inteiro em bytes.");
    }

    std::string str = val.get<std::string>();
    
    // Remover espaços em branco
    str.erase(std::remove_if(str.begin(), str.end(), ::isspace), str.end());
    if (str.empty()) {
        return 0;
    }

    // Identificar onde terminam os dígitos e começam as unidades
    size_t idx = 0;
    double number = std::stod(str, &idx);

    std::string unit = str.substr(idx);
    for (char& c : unit) {
        c = static_cast<char>(std::toupper(c));
    }

    size_t multiplier = 1ULL;

    // Suporte para GiB/GB, MiB/MB, KiB/KB e TiB/TB (base binária 1024)
    if (unit == "B" || unit.empty()) {
        multiplier = 1ULL;
    } else if (unit == "KB" || unit == "KIB" || unit == "K") {
        multiplier = 1024ULL;
    } else if (unit == "MB" || unit == "MIB" || unit == "M") {
        multiplier = 1024ULL * 1024ULL;
    } else if (unit == "GB" || unit == "GIB" || unit == "G") {
        multiplier = 1024ULL * 1024ULL * 1024ULL;
    } else if (unit == "TB" || unit == "TIB" || unit == "T") {
        multiplier = 1024ULL * 1024ULL * 1024ULL * 1024ULL;
    } else {
        throw std::invalid_argument("Unidade de memória desconhecida em max_allowed_bytes: " + unit);
    }

    return static_cast<size_t>(number * multiplier);
}


static CellType parse_cell_type(const std::string& s) {
    if (s == "SOLID")  return CellType::SOLID;
    if (s == "INLET")  return CellType::INLET;
    if (s == "OUTLET") return CellType::OUTLET;
    if (s == "SPONGE") return CellType::SPONGE;
    return CellType::FLUID;
}

static CollisionType parse_collision_type(const std::string& s){
    if (s == "BGK") return CollisionType::BGK;
    if (s == "SRT") return CollisionType::SRT;
    if (s == "TRT") return CollisionType::TRT;
    if (s == "MRT") return CollisionType::MRT;

    throw std::invalid_argument("Tipo de collisor nao defenido ou invalido.");

}

static BoundaryCondition parse_boundary_type(const std::string& s){
    if (s == "PERIODIC" ) return BoundaryCondition::PERIODIC;
    if (s == "DIRICHLET" ) return BoundaryCondition::DIRICHLET;
    if (s == "NEUMANN" ) return BoundaryCondition::NEUMANN;

    throw std::invalid_argument("Tipo boundary invalido.");
}

static BoundarySolver parse_boundary_solver(const std::string& s){
    if (s == "BOUNCEBACK") return BoundarySolver::BOUNCEBACK;
    if (s == "ZOUHE") return BoundarySolver::ZOUHE;

    throw std::invalid_argument("Tipo boundary solver invalido");

}

static TypeGen parse_type_gen(const std::string& s) {
    return (s == "FILE" || s == "file") ? TypeGen::FILE : TypeGen::PROCEDURAL;
}

static FormatType parse_format_type(const std::string& s) {
    if (s == "RAW" )       return FormatType::RAW;
    if (s == "PGM" )       return FormatType::PGM;
    if (s == "RECTANGLE" ) return FormatType::RECTANGLE;
    return FormatType::CIRCLE;
}

void Initializer::load_config_from_json(const std::string& json_filepath, Simconfig& cfg) {
    std::ifstream file(json_filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Erro ao abrir ficheiro JSON de configuração: " + json_filepath);
    }

    json j;
    file >> j;

    // 1. Parâmetros da Simulação
    const auto& sim = j.at("simulation");
    cfg.use_custom_init = sim.value("use_custom_init", false);
    cfg.nx              = sim.at("nx").get<int>();
    cfg.ny              = sim.at("ny").get<int>();
    cfg.ny              = sim.at("ny").get<int>();
    cfg.nz              = sim.value("nz", 1);
    cfg.grid_size       = static_cast<size_t>(cfg.nx) * cfg.ny * cfg.nz;
    cfg.grid_type       = sim.value("grid_type", 0);
    cfg.max_steps       = sim.at("max_steps").get<int>();
    cfg.export_freq     = sim.at("export_freq").get<int>();
    cfg.output_prefix   = sim.at("output_prefix").get<std::string>();
    
    // Leitura e conversão de max_allowed_bytes
    if (sim.contains("max_allowed_bytes")) {
        cfg.max_allowed_bytes = parse_memory_bytes(sim.at("max_allowed_bytes"));
    } else {
        cfg.max_allowed_bytes = 0; // ou o valor default que pretendes
    }

    cfg.rho0 = sim.value("rho0", 1.0);
    cfg.ux0  = sim.value("ux0", 0.0);
    cfg.uy0  = sim.value("uy0", 0.0);

    // 2. Modelo de Colisão
    if (j.contains("collision")) {
        const auto& col = j.at("collision");
        
        cfg.type_collision   = parse_collision_type(col.value("type", "BGK" ));
        
        if (cfg.type_collision != CollisionType::MRT ){
            cfg.tau_BGK          = col.value("tau_BGK", 0.8);
        }
        
        if (cfg.type_collision == CollisionType::TRT ){
            cfg.tau_TRT          = col.value("tau_TRT", 0.8);
        }

        if (cfg.type_collision == CollisionType::MRT ){
            cfg.tau_MRT_filename = col.value("tau_MRT_filename", "");

            if (cfg.type_collision == CollisionType::MRT && !cfg.tau_MRT_filename.empty()) {
                std::ifstream mrt_file(cfg.tau_MRT_filename);
                if (!mrt_file.is_open()) {
                    throw std::runtime_error("Não foi possível abrir tau_MRT_filename: " + cfg.tau_MRT_filename);
                }
                for (int i = 0; i < 9; ++i) {
                    mrt_file >> cfg.s_MRT[i];
                }
            }
        }
    }

    // 3. Condições de Fronteira nas Faces (x0, x1, y0, y1, z0, z1)
    const std::string face_keys[6] = {"x0", "x1", "y0", "y1", "z0", "z1"};
    if (j.contains("boundaries")) {
        const auto& bnds = j.at("boundaries");
        for (int i = 0; i < 6; ++i) {
            if (bnds.contains(face_keys[i])) {
                const auto& f = bnds.at(face_keys[i]);
                
                cfg.boundary_conditions[i]  = parse_boundary_type(f.value("condition", "PERIODIC"));
                cfg.type_boundary_solver[i] = parse_boundary_solver(f.value("solver", "BOUNCEBACK"));

                cfg.face_cell_type[i]       = parse_cell_type(f.value("cell_type", "FLUID"));
                cfg.face_vx[i]              = f.value("vx", 0.0f);
                cfg.face_vy[i]              = f.value("vy", 0.0f);
                cfg.face_rho[i]             = f.value("rho", 1.0f);
                cfg.face_sponge_width[i]    = f.value("sponge_width", 0);
                cfg.face_sponge_sigma[i]    = f.value("sponge_sigma_max", 0.0f);
            }
        }
    }

    // 4. Lista de Obstáculos / Objetos (Círculos, Meios-Círculos, Esponjas, PGM, RAW)
    cfg.obstacles.clear();
    if (j.contains("obstacles")) {
        for (const auto& item : j.at("obstacles")) {
            GeometryConfig obs{};
            obs.type              = parse_type_gen(item.value("type", "PROCEDURAL"));
            obs.format            = parse_format_type(item.value("format", "Circle"));
            obs.typeboundary      = parse_cell_type(item.value("typeboundary", "SOLID"));
            obs.filename          = item.value("filename", "");
            obs.obj_nx            = item.value("obj_nx", 0);
            obs.obj_ny            = item.value("obj_ny", 0);
            obs.obj_nz            = item.value("obj_nz", 1);
            obs.offset_x          = item.value("offset_x", 0);
            obs.offset_y          = item.value("offset_y", 0);
            obs.offset_z          = item.value("offset_z", 0);
            obs.radius            = item.value("radius", 0.0);
            obs.sponge_width      = item.value("sponge_width", 0.0);
            obs.sponge_sigma_max  = item.value("sponge_sigma_max", 0.0f);
            obs.vx                = item.value("vx", static_cast<float>(cfg.ux0));
            obs.vy                = item.value("vy", static_cast<float>(cfg.uy0));
            obs.rho               = item.value("rho", static_cast<float>(cfg.rho0));

            cfg.obstacles.push_back(obs);
        }
    }
}

void Initializer::initialize(Grid& grid, const Simconfig& cfg) {

    utils::check_memory_availability(cfg.max_allowed_bytes);

    // FASE 1: Alocar apenas o mapa de flags da geometria
    grid.allocate_geometry(cfg.nx, cfg.ny);

    // FASE 2: Estampar geometria e construir listas de fronteira/esponjas
    GeometryBuilder::stamp_geometry(grid, cfg);
    GeometryBuilder::build_boundary_lists(grid, cfg);

    // FASE 3: Alocar macro-bloco de 17 planos de populações e inicializar no equilíbrio
    grid.allocate_populations(cfg.max_allowed_bytes);

    if (cfg.use_custom_init) {
        Initializer::apply_custom_initialization(grid, cfg);
        return;
    }

    const double rho  = cfg.rho0;
    const double ux   = cfg.ux0;
    const double uy   = cfg.uy0;
    const double u_sq = ux * ux + uy * uy;
    const size_t pitch = grid.plane_pitch;

    const double feq0 = D2Q9::compute_feq0(rho, u_sq);
    double feq[8];
    for (int d = 0; d < 8; ++d) {
        feq[d] = D2Q9::compute_feq(d, rho, ux, uy, u_sq);
    }

    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < grid.size; ++i) {
        if (grid.flag_field_ptr[i] == CellType::SOLID) {
            grid.f0_ptr[i] = D2Q9::compute_feq0(rho, 0.0);
            for (int d = 0; d < 8; ++d) {
                grid.f1_ptr[d * pitch + i] = D2Q9::compute_feq(d, rho, 0.0, 0.0, 0.0);
                grid.f2_ptr[d * pitch + i] = grid.f1_ptr[d * pitch + i];
            }
        } else {
            grid.f0_ptr[i] = feq0;
            for (int d = 0; d < 8; ++d) {
                grid.f1_ptr[d * pitch + i] = feq[d];
                grid.f2_ptr[d * pitch + i] = feq[d];
            }
        }
    }
}
