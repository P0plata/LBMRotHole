#include "config_init.hpp"
#include "grid.hpp"
#include "solver_d2q9.hpp"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <omp.h>

int main(int argc, char* argv[]) {
    try {
        const std::string config_path = (argc > 1) ? argv[1] : "config.json";

        std::cout << "========================================================\n";
        std::cout << "  LBM D2Q9 High-Performance Solver\n";
        std::cout << "========================================================\n";
        std::cout << "[1/3] A carregar configuracao: " << config_path << "...\n";

        Simconfig cfg{};
        Initializer::load_config_from_json(config_path, cfg);

        std::cout << "[2/3] A construir geometria e alocar memoria SoA...\n";
        auto t_init_start = std::chrono::high_resolution_clock::now();

        Grid grid{};
        Initializer::initialize(grid, cfg);

        auto t_init_end = std::chrono::high_resolution_clock::now();
        double init_sec = std::chrono::duration<double>(t_init_end - t_init_start).count();

        double mem_mb = static_cast<double>(grid.memory_usage_bytes()) / (1024.0 * 1024.0);
        std::cout << "      -> Grelha:         " << grid.nx << " x " << grid.ny 
                  << " (" << grid.size << " celulas)\n";
        std::cout << "      -> Memoria Total:  " << std::fixed << std::setprecision(2) << mem_mb << " MiB\n";
        std::cout << "      -> Wall Links:     " << grid.wall_links.size() << "\n";
        std::cout << "      -> Sponge Nodes:   " << grid.sponge_nodes.size() << "\n";
        std::cout << "      -> Inlet / Outlet: " << grid.inlet_nodes.size() << " / " << grid.outlet_nodes.size() << "\n";
        std::cout << "      -> Threads OpenMP: " << omp_get_max_threads() << "\n";
        std::cout << "      -> Tempo de Init:  " << std::setprecision(3) << init_sec << " s\n";

        std::cout << "[3/3] A executar simulacao (" << cfg.max_steps << " passos)...\n";
        auto t_sim_start = std::chrono::high_resolution_clock::now();

        // Seleciona a versão pré-compilada do template (zero branches lá dentro!)
        switch (cfg.type_collision) {
            case CollisionType::BGK:
                SolverD2Q9::run_simulation<CollisionType::BGK>(grid, cfg);
                break;
            case CollisionType::TRT:
                SolverD2Q9::run_simulation<CollisionType::TRT>(grid, cfg);
                break;
            case CollisionType::MRT:
                SolverD2Q9::run_simulation<CollisionType::MRT>(grid, cfg);
                break;
        }

        auto t_sim_end = std::chrono::high_resolution_clock::now();
        double sim_sec = std::chrono::duration<double>(t_sim_end - t_sim_start).count();
        double mlups   = (static_cast<double>(grid.size) * cfg.max_steps) / (sim_sec * 1e6);

        std::cout << "========================================================\n";
        std::cout << "  Simulacao Concluida com Sucesso!\n";
        std::cout << "  Tempo Total:  " << std::fixed << std::setprecision(3) << sim_sec << " s\n";
        std::cout << "  Performance:  " << std::fixed << std::setprecision(2) << mlups << " MLUPS\n";
        std::cout << "========================================================\n";

    } catch (const std::exception& e) {
        std::cerr << "\n[ERRO FATAL] " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}