#include "utils.hpp"

// Obtém a memória física disponível real em bytes via /proc/meminfo
size_t  utils::get_system_available_ram() {
    std::ifstream meminfo("/proc/meminfo");
    if (!meminfo.is_open()) {
        throw std::runtime_error("Não foi possível ler /proc/meminfo.");
    }

    std::string key;
    size_t value_kb = 0;
    std::string unit;

    while (meminfo >> key >> value_kb >> unit) {
        if (key == "MemAvailable:") {
            return value_kb * 1024ULL; // Converte kB para bytes
        }
    }

    throw std::runtime_error("Campo MemAvailable não encontrado em /proc/meminfo.");
}

// Valida se a RAM disponível é pelo menos 95% do limite configurado
void utils::check_memory_availability(size_t max_allowed_bytes) {
    if (max_allowed_bytes == 0) return;

    size_t available_bytes = get_system_available_ram();
    size_t required_bytes  = static_cast<size_t>(max_allowed_bytes * 0.95);

    if (available_bytes < required_bytes) {
        double avail_gb = static_cast<double>(available_bytes) / (1024.0 * 1024.0 * 1024.0);
        double req_gb   = static_cast<double>(required_bytes) / (1024.0 * 1024.0 * 1024.0);
        double max_gb   = static_cast<double>(max_allowed_bytes) / (1024.0 * 1024.0 * 1024.0);

        std::ostringstream oss;
        oss << "Memória RAM insuficiente! "
            << "Configurado: " << max_gb << " GB, "
            << "Mínimo exigido (95%): " << req_gb << " GB, "
            << "Disponível no sistema: " << avail_gb << " GB.";
        throw std::runtime_error(oss.str());
    }
}

size_t utils::get_system_l3_cache_bytes() {
    namespace fs = std::filesystem;
    size_t total_l3_bytes = 0;
    std::set<std::string> seen_shared_maps;

    const fs::path cpu_base("/sys/devices/system/cpu");
    if (fs::exists(cpu_base)) {
        for (const auto& entry : fs::directory_iterator(cpu_base)) {
            const std::string name = entry.path().filename().string();
            if (name.rfind("cpu", 0) == 0 && name.size() > 3 && std::isdigit(name[3])) {
                fs::path index3 = entry.path() / "cache" / "index3";
                fs::path size_file = index3 / "size";
                fs::path map_file  = index3 / "shared_cpu_list";

                if (fs::exists(size_file) && fs::exists(map_file)) {
                    std::ifstream mf(map_file);
                    std::string shared_list;
                    std::getline(mf, shared_list);

                    // Conta cada bloco físico de L3 (ex: os 2 CCDs de 32MB do 9950X = 64MB) apenas uma vez
                    if (seen_shared_maps.insert(shared_list).second) {
                        std::ifstream sf(size_file);
                        size_t val = 0;
                        char unit = 'K';
                        if (sf >> val >> unit) {
                            if (unit == 'K' || unit == 'k') total_l3_bytes += val * 1024ULL;
                            else if (unit == 'M' || unit == 'm') total_l3_bytes += val * 1024ULL * 1024ULL;
                        }
                    }
                }
            }
        }
    }
    // Fallback de segurança (32 MB) caso não consiga ler o /sys
    return (total_l3_bytes > 0) ? total_l3_bytes : (32ULL * 1024ULL * 1024ULL);
}



void utils::audit_grid_integrity(const Grid& grid) {
    size_t nan_or_inf_count = 0;
    size_t neg_pop_count    = 0;
    double total_mass       = 0.0;
    const size_t pitch      = grid.plane_pitch;

    #pragma omp parallel for reduction(+:nan_or_inf_count, neg_pop_count, total_mass) schedule(static)
    for (size_t i = 0; i < grid.size; ++i) {
        if (grid.flag_field_ptr[i] == CellType::SOLID) continue;
        double r = grid.f0_ptr[i];
        if (!std::isfinite(r)) nan_or_inf_count++;
        if (r < 0.0) neg_pop_count++;
        for (int d = 0; d < 8; ++d) {
            double fd = grid.f1_ptr[d * pitch + i];
            if (!std::isfinite(fd)) nan_or_inf_count++;
            if (fd < 0.0) neg_pop_count++;
            r += fd;
        }
        total_mass += r;
    }

    std::cout << "  ------------------------------------------------------\n";
    std::cout << "  -> Massa Fluida Total:  " << std::scientific << std::setprecision(6) << total_mass << "\n";
    std::cout << "  -> NaN/Inf | Pop < 0:   " << nan_or_inf_count << " | " << neg_pop_count
              << ((nan_or_inf_count == 0 && neg_pop_count == 0) ? " [100% INTEGRO]\n" : " [ALERTA: INSTABILIDADE!]\n");
}