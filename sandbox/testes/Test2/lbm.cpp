#include "lbm.hpp"

namespace fs = std::filesystem;

LBM::LBM(int x, int y) {
    grid.allocate(x, y);
}

LBM::~LBM() {
    grid.deallocate();
}

void LBM::initialize(){
    //initialize
}

void LBM::stream_compute_collide(){
    // stream, compute macroscopic, compute equilibrium, collide, apply forcing

}

void LBM::apply_boundary(){
    // apply boundary conditions
}

void LBM::export_to_vti(int step) {
    fs::create_directories("data");
    char filename[128];
    std::snprintf(filename, sizeof(filename), "data/sim_%06d.vti", step);

    // O 'wb' (write binary) é crítico no Windows/WSL para não corromper bits
    FILE *fp = std::fopen(filename, "wb");
    if (!fp) return;

    size_t num_cells = grid.size;
    
    // 1. Buffers temporários para contiguidade de memória no fwrite
    // São destruídos automaticamente no fim da função.
    std::vector<double> rho(num_cells, 0.0);
    std::vector<double> vel(num_cells * 3, 0.0); // O ParaView prefere vetores 3D (Z = 0)

    const double cx[8] = {1, 0, -1, 0, 1, -1, -1, 1};
    const double cy[8] = {0, 1, 0, -1, 1, 1, -1, -1};

    // Calcular on-the-fly
    for (size_t i = 0; i < num_cells; i++) {
        if (grid.flag_field_ptr[i] == 1) { // Assumindo 1 = Sólido (Bounce-Back)
            rho[i] = 1.0;
            // vel já está inicializado a 0.0
            continue; 
        }

        double r = grid.f0_ptr[i];
        double ux = 0.0, uy = 0.0;
        
        for (int d = 0; d < 8; d++) {
            double f = grid.f1_ptr[i + d * grid.plane_pitch];
            r += f;
            ux += f * cx[d];
            uy += f * cy[d];
        }
        
        rho[i] = r;
        vel[i * 3 + 0] = ux / r;
        vel[i * 3 + 1] = uy / r;
        vel[i * 3 + 2] = 0.0;
    }

    // 2. Calcular o tamanho dos blocos binários para o cabeçalho XML
    uint32_t bytes_rho  = num_cells * sizeof(double);
    uint32_t bytes_vel  = num_cells * 3 * sizeof(double);
    uint32_t bytes_flag = num_cells * sizeof(uint8_t);

    // No modo "appended", o VTK precisa de saber a que distância (offset) começa cada array
    // O sizeof(uint32_t) extra é porque o VTK exige que cada bloco binário comece com o seu próprio tamanho
    uint32_t offset_rho  = 0;
    uint32_t offset_vel  = offset_rho + bytes_rho + sizeof(uint32_t);
    uint32_t offset_flag = offset_vel + bytes_vel + sizeof(uint32_t);

    // 3. Escrever Cabeçalho XML
    std::fprintf(fp, "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
    std::fprintf(fp, "  <ImageData WholeExtent=\"0 %d 0 %d 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">\n", grid.nx - 1, grid.ny - 1);
    std::fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 0\">\n", grid.nx - 1, grid.ny - 1);
    std::fprintf(fp, "      <PointData Scalars=\"density\" Vectors=\"velocity\">\n");
    
    std::fprintf(fp, "        <DataArray type=\"Float64\" Name=\"density\" format=\"appended\" offset=\"%u\"/>\n", offset_rho);
    std::fprintf(fp, "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%u\"/>\n", offset_vel);
    std::fprintf(fp, "        <DataArray type=\"UInt8\" Name=\"flag_field\" format=\"appended\" offset=\"%u\"/>\n", offset_flag);
    
    std::fprintf(fp, "      </PointData>\n");
    std::fprintf(fp, "    </Piece>\n");
    std::fprintf(fp, "  </ImageData>\n");
    
    // 4. Iniciar secção binária (VTK exige o marcador '_')
    std::fprintf(fp, "  <AppendedData encoding=\"raw\">\n");
    std::fprintf(fp, "_"); 

    // 5. Despejar a memória bruta para o disco
    std::fwrite(&bytes_rho, sizeof(uint32_t), 1, fp);
    std::fwrite(rho.data(), 1, bytes_rho, fp);

    std::fwrite(&bytes_vel, sizeof(uint32_t), 1, fp);
    std::fwrite(vel.data(), 1, bytes_vel, fp);

    std::fwrite(&bytes_flag, sizeof(uint32_t), 1, fp);
    std::fwrite(grid.flag_field_ptr, 1, bytes_flag, fp);

    // 6. Fechar as tags (a nova linha é intencional para fechar de forma limpa)
    std::fprintf(fp, "\n  </AppendedData>\n");
    std::fprintf(fp, "</VTKFile>\n");

    std::fclose(fp);
}  


