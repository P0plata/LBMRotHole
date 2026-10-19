#include "vtk_exporter.hpp"

namespace fs = std::filesystem;

void VTKExporter::export_to_vti(const Grid& grid, const int step, const std::string& prefix) {
    fs::path p(prefix);
    if (p.has_parent_path()) {
        fs::create_directories(p.parent_path());
    }

    char filename[256];
    std::snprintf(filename, sizeof(filename), "%s_%06d.vti", prefix.c_str(), step);

    FILE *fp = std::fopen(filename, "wb");
    if (!fp) return;

    // Buffer de escrita contínua de 1 MB
    std::setvbuf(fp, nullptr, _IOFBF, 1024 * 1024);

    const size_t num_cells = grid.size;
    const size_t pitch     = grid.plane_pitch;

    const uint32_t bytes_rho  = static_cast<uint32_t>(num_cells * sizeof(float));
    const uint32_t bytes_vel  = static_cast<uint32_t>(num_cells * 3 * sizeof(float));
    const uint32_t bytes_flag = static_cast<uint32_t>(num_cells * sizeof(uint8_t));

    const uint32_t offset_rho  = 0;
    const uint32_t offset_vel  = offset_rho + bytes_rho + sizeof(uint32_t);
    const uint32_t offset_flag = offset_vel + bytes_vel + sizeof(uint32_t);

    // CRÍTICO: Sem "<?xml version="1.0"?>" no topo para o ParaView aceitar encoding="raw"!
    std::fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    std::fprintf(fp, "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", grid.nx - 1, grid.ny - 1);
    std::fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n", grid.nx - 1, grid.ny - 1);
    std::fprintf(fp, "      <PointData Scalars=\"density\" Vectors=\"velocity\">\n");
    std::fprintf(fp, "        <DataArray type=\"Float32\" Name=\"density\" format=\"appended\" offset=\"%u\"/>\n", offset_rho);
    std::fprintf(fp, "        <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%u\"/>\n", offset_vel);
    std::fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"flag_field\" format=\"appended\" offset=\"%u\"/>\n", offset_flag);
    std::fprintf(fp, "      </PointData>\n");
    std::fprintf(fp, "    </Piece>\n");
    std::fprintf(fp, "  </ImageData>\n");
    std::fprintf(fp, "  <AppendedData encoding=\"raw\">_");


    constexpr size_t CHUNK_SIZE = 65536;
    std::vector<float> chunk_buf(CHUNK_SIZE * 3);

    // 1. Bloco de Densidade (Escrita sequencial pura, sem fseek)
    std::fwrite(&bytes_rho, sizeof(uint32_t), 1, fp);
    for (size_t base = 0; base < num_cells; base += CHUNK_SIZE) {
        const size_t count = std::min(CHUNK_SIZE, num_cells - base);
        for (size_t k = 0; k < count; ++k) {
            const size_t i = base + k;
            if (grid.flag_field_ptr[i] == CellType::SOLID) {
                chunk_buf[k] = 1.0f;
                continue;
            }
            double r = grid.f0_ptr[i];
            for (int d = 0; d < 8; ++d) r += grid.f1_ptr[d * pitch + i];
            chunk_buf[k] = static_cast<float>(r);
        }
        std::fwrite(chunk_buf.data(), sizeof(float), count, fp);
    }

    // 2. Bloco de Velocidade (Escrita sequencial pura, sem fseek)
    std::fwrite(&bytes_vel, sizeof(uint32_t), 1, fp);
    for (size_t base = 0; base < num_cells; base += CHUNK_SIZE) {
        const size_t count = std::min(CHUNK_SIZE, num_cells - base);
        for (size_t k = 0; k < count; ++k) {
            const size_t i = base + k;
            if (grid.flag_field_ptr[i] == CellType::SOLID) {
                chunk_buf[k * 3 + 0] = 0.0f;
                chunk_buf[k * 3 + 1] = 0.0f;
                chunk_buf[k * 3 + 2] = 0.0f;
                continue;
            }
            double r = grid.f0_ptr[i];
            double ux = 0.0, uy = 0.0;
            for (int d = 0; d < 8; ++d) {
                const double f = grid.f1_ptr[d * pitch + i];
                r  += f;
                ux += f * D2Q9::cx[d];
                uy += f * D2Q9::cy[d];
            }
            const double inv_r = 1.0 / r;
            chunk_buf[k * 3 + 0] = static_cast<float>(ux * inv_r);
            chunk_buf[k * 3 + 1] = static_cast<float>(uy * inv_r);
            chunk_buf[k * 3 + 2] = 0.0f;
        }
        std::fwrite(chunk_buf.data(), sizeof(float), count * 3, fp);
    }

    // 3. Bloco de Flag Field
    std::fwrite(&bytes_flag, sizeof(uint32_t), 1, fp);
    std::fwrite(grid.flag_field_ptr, sizeof(uint8_t), num_cells, fp);

    std::fprintf(fp, "\n  </AppendedData>\n");
    std::fprintf(fp, "</VTKFile>\n");
    std::fclose(fp);
}









