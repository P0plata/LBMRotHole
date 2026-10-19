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
