#!/usr/bin/env python3
import os
import re
import json
import struct
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Caminho para o executável compilado
EXECUTABLE = "./testes/Test2/build/run_test"
OUT_DIR = "testes/Test2/validation_results"
os.makedirs(OUT_DIR, exist_ok=True)


def read_vti_binary(filepath, nx, ny):
    """Lê diretamente os blocos binários raw do ficheiro .vti gerado pelo C++."""
    with open(filepath, "rb") as f:
        content = f.read()

    marker = b'<AppendedData encoding="raw">_'
    pos = content.find(marker)
    if pos == -1:
        raise ValueError(f"Bloco AppendedData não encontrado em {filepath}")

    offset = pos + len(marker)
    num_cells = nx * ny

    # 1. Densidade (float32)
    bytes_rho = struct.unpack_from("<I", content, offset)[0]
    offset += 4
    rho = np.frombuffer(content, dtype=np.float32, count=num_cells, offset=offset).reshape((ny, nx))
    offset += bytes_rho

    # 2. Velocidade (3 x float32: ux, uy, uz)
    bytes_vel = struct.unpack_from("<I", content, offset)[0]
    offset += 4
    vel = np.frombuffer(content, dtype=np.float32, count=num_cells * 3, offset=offset).reshape((ny, nx, 3))
    offset += bytes_vel

    # 3. Flag field (uint8)
    bytes_flag = struct.unpack_from("<I", content, offset)[0]
    offset += 4
    flag = np.frombuffer(content, dtype=np.uint8, count=num_cells, offset=offset).reshape((ny, nx))

    return rho.astype(np.float64), vel[:, :, 0].astype(np.float64), vel[:, :, 1].astype(np.float64), flag


def run_lbm_case(config_dict, name):
    """Gera o JSON, executa o solver C++ e devolve o output do terminal."""
    json_path = os.path.join(OUT_DIR, f"{name}.json")
    with open(json_path, "w") as f:
        json.dump(config_dict, f, indent=2)

    res = subprocess.run([EXECUTABLE, json_path], capture_output=True, text=True, check=True)
    return res.stdout


def compute_l2_error(num, exact):
    """Calcula a norma relativa do erro L2."""
    return np.sqrt(np.sum((num - exact) ** 2) / np.sum(exact ** 2))


# ==============================================================================
# TESTE 1: EQUIVALÊNCIA ALGÉBRICA (BGK vs TRT vs MRT) & CONSERVAÇÃO DE MASSA
# ==============================================================================
def test_1_collider_consistency():
    print("\n[TESTE 1/4] A validar Consistencia BGK vs TRT vs MRT e Conservacao de Massa...")
    nx, ny = 80, 80
    tau = 0.70
    s_val = 1.0 / tau

    # Ficheiro de relaxamento MRT com todos os s_k = 1 / tau
    mrt_file = os.path.join(OUT_DIR, "mrt_equiv.txt")
    with open(mrt_file, "w") as f:
        f.write(" ".join([f"{s_val:.12f}"] * 9) + "\n")

    results = {}
    for col in ["BGK", "TRT", "MRT"]:
        prefix = os.path.join(OUT_DIR, f"t1_{col}")
        cfg = {
            "simulation": {
                "use_custom_init": True, "nx": nx, "ny": ny, "max_steps": 500,
                "export_freq": 500, "output_prefix": prefix, "rho0": 100.0, "ux0": 0.02, "uy0": 0.0
            },
            "collision": {"type": col, "tau_BGK": tau, "tau_TRT": tau, "tau_MRT_filename": mrt_file},
            "forces": {"use_forces": True, "center_x": 40, "center_y": 40, "g_mass": 50.0},
            "boundaries": {
                "x0": {"condition": "PERIODIC"}, "x1": {"condition": "PERIODIC"},
                "y0": {"condition": "PERIODIC"}, "y1": {"condition": "PERIODIC"}
            },
            "obstacles": []
        }
        out = run_lbm_case(cfg, f"cfg_t1_{col}")
        integro = "[100% INTEGRO]" in out
        rho_0, _, _, _ = read_vti_binary(f"{prefix}_000000.vti", nx, ny)
        rho_f, ux_f, uy_f, _ = read_vti_binary(f"{prefix}_000500.vti", nx, ny)
        mass_err = abs(np.sum(rho_f) - np.sum(rho_0)) / np.sum(rho_0)
        results[col] = (rho_f, ux_f, uy_f, mass_err, integro)
        print(f"   -> {col:3s}: Integro={integro} | Erro Conservacao Massa = {mass_err:.2e}")

    diff_trt = compute_l2_error(results["TRT"][1], results["BGK"][1])
    diff_mrt = compute_l2_error(results["MRT"][1], results["BGK"][1])
    print(f"   -> Diferenca L2 (TRT vs BGK): {diff_trt:.2e}")
    print(f"   -> Diferenca L2 (MRT O(u^4) vs BGK O(u^2)): {diff_mrt:.2e}")
    passed = (diff_trt < 1e-6) and (diff_mrt < 1e-3) and all(r[4] for r in results.values())
    return passed


# ==============================================================================
# TESTE 2: DECAIMENTO DE ONDA DE CORTE (VALIDAÇÃO DA VISCOSIDADE FÍSICA)
# ==============================================================================
def test_2_shear_wave_viscosity():
    print("\n[TESTE 2/4] A validar Viscosidade Cinematica (Decaimento de Onda de Corte)...")
    nx, ny = 16, 128
    u0 = 0.01
    steps = 2000
    export_f = 200
    taus = [0.58, 0.70, 0.90]
    k_wave = 2.0 * np.pi / ny

    plt.figure(figsize=(8, 5))
    max_l2_err = 0.0

    for tau in taus:
        nu_theo = (1.0 / 3.0) * (tau - 0.5)
        prefix = os.path.join(OUT_DIR, f"t2_tau_{tau:.2f}")
        cfg = {
            "simulation": {
                "use_custom_init": True, "nx": nx, "ny": ny, "max_steps": steps,
                "export_freq": export_f, "output_prefix": prefix, "rho0": 100.0, "ux0": u0, "uy0": 0.0
            },
            "collision": {"type": "BGK", "tau_BGK": tau},
            "forces": {"use_forces": False},
            "boundaries": {
                "x0": {"condition": "PERIODIC"}, "x1": {"condition": "PERIODIC"},
                "y0": {"condition": "PERIODIC"}, "y1": {"condition": "PERIODIC"}
            },
            "obstacles": []
        }
        run_lbm_case(cfg, f"cfg_t2_{tau:.2f}")

        t_list = list(range(0, steps + 1, export_f))
        u_lbm = []
        for s in t_list:
            _, ux, _, _ = read_vti_binary(f"{prefix}_{s:06d}.vti", nx, ny)
            u_lbm.append(np.max(ux))

        u_lbm = np.array(u_lbm)
        u_exact = u0 * np.exp(-nu_theo * (k_wave ** 2) * np.array(t_list))
        err_l2 = compute_l2_error(u_lbm, u_exact)
        max_l2_err = max(max_l2_err, err_l2)
        print(f"   -> tau = {tau:.2f} (nu = {nu_theo:.4f}): Erro L2 = {err_l2 * 100:.4f}%")

        t_fine = np.linspace(0, steps, 200)
        plt.plot(t_fine, u0 * np.exp(-nu_theo * (k_wave ** 2) * t_fine), "-", label=f"Analítico ($\\tau={tau}$)")
        plt.plot(t_list, u_lbm, "o", markersize=5, label=f"LBM ($\\tau={tau}$, Erro={err_l2*100:.2f}%)")

    plt.title("Validação da Viscosidade: Decaimento de Onda de Corte de Navier-Stokes")
    plt.xlabel("Passos de Tempo ($t$)")
    plt.ylabel("Velocidade Máxima $u_{\\max}(t)$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "1_validacao_viscosidade.png"), dpi=200)
    plt.close()
    return max_l2_err < 0.005


# ==============================================================================
# TESTE 3: CANAL DE POISEUILLE + ESTUDO DE CONVERGÊNCIA DE MALHA (2ª ORDEM)
# ==============================================================================
def test_3_poiseuille_and_mesh_convergence():
    print("\n[TESTE 3/7] A executar Poiseuille 2D e Estudo de Convergencia de Malha...")
    tau = 1.00  # <-- tau = 1.00 revela o erro espacial O(dx^2) e evita o ponto zero de tau = 0.75!
    nu = (1.0 / 3.0) * (tau - 0.5)
    rho0 = 100.0
    u_max_target = 0.005

    H_list = [16, 32, 64, 128]
    errors_l2 = []

    plt.figure(figsize=(8, 5))

    for H in H_list:
        ny = H + 2
        nx = 8
        fx_val = 8.0 * rho0 * nu * u_max_target / (H ** 2)
        # 5x o tempo difusivo H^2 / nu é suficiente para atingir regime estacionário (corre muito mais rápido!)
        steps = int(5.0 * (H ** 2) / nu)
        prefix = os.path.join(OUT_DIR, f"t3_H{H}")

        cfg = {
            "simulation": {
                "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
                "export_freq": steps, "output_prefix": prefix, "rho0": rho0, "ux0": 0.0, "uy0": 0.0
            },
            "collision": {"type": "BGK", "tau_BGK": tau},
            "forces": {"use_forces": True, "center_x": -1, "center_y": -1, "g_mass": fx_val},
            "boundaries": {
                "x0": {"condition": "PERIODIC", "cell_type": "FLUID"},
                "x1": {"condition": "PERIODIC", "cell_type": "FLUID"},
                "y0": {"condition": "DIRICHLET", "cell_type": "SOLID"},
                "y1": {"condition": "DIRICHLET", "cell_type": "SOLID"}
            },
            "obstacles": []
        }
        run_lbm_case(cfg, f"cfg_t3_H{H}")

        _, ux, _, _ = read_vti_binary(f"{prefix}_{steps:06d}.vti", nx, ny)
        ux_profile = ux[1:-1, nx // 2]

        y_phys = np.arange(1, ny - 1) - 0.5
        ux_exact = (fx_val / (2.0 * rho0 * nu)) * y_phys * (H - y_phys)

        err = compute_l2_error(ux_profile, ux_exact)
        errors_l2.append(err)
        print(f"   -> Malha H = {H:3d} (ny={ny:3d}, {steps:6d} passos): Erro L2 = {err * 100:.4f}%")

        if H == H_list[-1]:
            y_norm = y_phys / H
            plt.plot(y_norm, ux_exact / u_max_target, "k-", lw=2, label="Solução Analítica de Poiseuille")
            plt.plot(y_norm[::2], ux_profile[::2] / u_max_target, "ro", label=f"LBM ($H={H}$, Erro $L_2={err*100:.3f}\\%$)")

    plt.title("Perfil de Velocidade no Canal de Poiseuille 2D (Bounce-Back + Força de Guo)")
    plt.xlabel("Posição Normalizada no Canal ($y / H$)")
    plt.ylabel("Velocidade Normalizada ($u_x / U_{\\max}$)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "2_validacao_poiseuille.png"), dpi=200)
    plt.close()

    # Cálculo da ordem de convergência espacial p
    slopes = np.polyfit(np.log(H_list), np.log(errors_l2), 1)
    order_p = -slopes[0]
    print(f"   -> Ordem de Convergencia Espacial Medida: p = {order_p:.2f} (Teorico LBM = 2.00)")

    plt.figure(figsize=(7, 5))
    H_arr = np.array(H_list, dtype=np.float64)
    plt.loglog(H_arr, errors_l2, "bo-", lw=2, markersize=7, label=f"Erro $L_2$ LBM (Ordem $p = {order_p:.2f}$)")
    ref_line = errors_l2[0] * (H_arr[0] / H_arr) ** 2
    plt.loglog(H_arr, ref_line, "k--", label="Declive Teórico de 2.ª Ordem $\\mathcal{O}(\\Delta x^2)$")
    plt.title("Estudo de Convergência de Malha (Ordem de Precisão Espacial)")
    plt.xlabel("Resolução do Canal ($H = N_y - 2$)")
    plt.ylabel("Erro Relativo $L_2$")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "3_convergencia_malha_loglog.png"), dpi=200)
    plt.close()

    return (order_p > 1.85) and (errors_l2[-1] < 0.005)

# ==============================================================================
# TESTE 4: EQUILÍBRIO HIDROSTÁTICO GRAVITACIONAL (FORÇAS DE GUO 1/r^2)
# ==============================================================================
def test_4_hydrostatic_gravity_well():
    print("\n[TESTE 4/4] A validar Equilibrio Hidrostatico no Poco Gravitacional (Guo)...")
    nx, ny = 200, 200
    xc, yc = 100, 100
    rho0 = 100.0
    g_mass = 80.0
    steps = 4000
    prefix = os.path.join(OUT_DIR, "t4_grav")

    cfg = {
        "simulation": {
            "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
            "export_freq": steps, "output_prefix": prefix, "rho0": rho0, "ux0": 0.0, "uy0": 0.0
        },
        "collision": {"type": "BGK", "tau_BGK": 0.80},
        "forces": {"use_forces": True, "center_x": xc, "center_y": yc, "g_mass": g_mass},
        "boundaries": {
            "x0": {"condition": "PERIODIC", "rho": rho0, "sponge_width": 15, "sponge_sigma_max": 0.2},
            "x1": {"condition": "PERIODIC", "rho": rho0, "sponge_width": 15, "sponge_sigma_max": 0.2},
            "y0": {"condition": "PERIODIC", "rho": rho0, "sponge_width": 15, "sponge_sigma_max": 0.2},
            "y1": {"condition": "PERIODIC", "rho": rho0, "sponge_width": 15, "sponge_sigma_max": 0.2}
        },
        "obstacles": []
    }
    run_lbm_case(cfg, "cfg_t4_grav")

    rho_f, ux_f, uy_f, _ = read_vti_binary(f"{prefix}_{steps:06d}.vti", nx, ny)

    # Corte horizontal passando pelo centro y = yc (para r <= 75, antes da esponja das bordas)
    x_idx = np.arange(xc, xc + 75)
    r_vals = x_idx - xc
    rho_lbm = rho_f[yc, x_idx]

    # Solução Hidrostática Analítica Exata: c_s^2 * grad(rho) = F(r) => rho(r) = rho_inf + 3 * g_mass / sqrt(r^2 + 225)
    potential = 3.0 * g_mass / np.sqrt(r_vals ** 2 + 225.0)
    rho_inf_fit = rho_lbm[-1] - potential[-1]
    rho_exact = rho_inf_fit + potential

    err_l2 = compute_l2_error(rho_lbm - rho_inf_fit, rho_exact - rho_inf_fit)
    max_spurious_vel = np.max(np.hypot(ux_f, uy_f))
    print(f"   -> Erro L2 no Perfil de Densidade Gravitacional: {err_l2 * 100:.3f}%")
    print(f"   -> Velocidade Espuria Maxima em Repouso:         {max_spurious_vel:.2e}")

    plt.figure(figsize=(8, 5))
    plt.plot(r_vals, rho_exact, "k-", lw=2, label="Solução Analítica Hidrostática ($\\nabla p = \\mathbf{F}$)")
    plt.plot(r_vals[::2], rho_lbm[::2], "ro", label=f"LBM c/ Força de Guo (Erro $L_2 = {err_l2*100:.2f}\\%$)")
    plt.title("Validação do Poço Gravitacional em Equilíbrio Hidrostático")
    plt.xlabel("Distância Radial ao Centro $r$ (células)")
    plt.ylabel("Densidade Estacionária $\\rho(r)$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "4_validacao_gravidade_hidrostatica.png"), dpi=200)
    plt.close()

    return err_l2 < 0.02

# ==============================================================================
# TESTE 5: FRONTEIRAS ABERTAS DE ZOU-HE (INLET VELOCIDADE + OUTLET PRESSÃO)
# ==============================================================================
def test_5_zou_he_open_boundaries():
    print("\n[TESTE 5/7] A validar Fronteiras Abertas de Zou-He (Inlet + Outlet)...")
    nx, ny = 120, 32  # Comprimento ideal para advecção rápida
    H = ny - 2
    u_inlet = 0.02
    tau = 0.80
    steps = 10000     # Transiente estabiliza por completo
    prefix = os.path.join(OUT_DIR, "t5_zouhe")

    cfg = {
        "simulation": {
            "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
            "export_freq": steps, "output_prefix": prefix, "rho0": 100.0, "ux0": u_inlet, "uy0": 0.0
        },
        "collision": {"type": "BGK", "tau_BGK": tau},
        "forces": {"use_forces": False},
        "boundaries": {
            "x0": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "INLET",  "vx": u_inlet, "vy": 0.0, "rho": 100.0},
            "x1": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "OUTLET", "vx": u_inlet, "vy": 0.0, "rho": 100.0},
            "y0": {"condition": "DIRICHLET", "cell_type": "SOLID"},
            "y1": {"condition": "DIRICHLET", "cell_type": "SOLID"}
        },
        "obstacles": []
    }
    out = run_lbm_case(cfg, "cfg_t5_zouhe")
    integro = "[100% INTEGRO]" in out

    rho_f, ux_f, _, _ = read_vti_binary(f"{prefix}_{steps:06d}.vti", nx, ny)

    x_probe = int(0.85 * nx)
    ux_prof = ux_f[1:-1, x_probe]

    u_mean_measured = np.mean(ux_prof)
    y_phys = np.arange(1, ny - 1) - 0.5
    ux_exact = 1.5 * u_mean_measured * (4.0 * y_phys * (H - y_phys) / (H ** 2))

    err_l2 = compute_l2_error(ux_prof, ux_exact)
    flux_in  = np.sum(rho_f[1:-1, int(0.20 * nx)] * ux_f[1:-1, int(0.20 * nx)])
    flux_out = np.sum(rho_f[1:-1, x_probe]        * ux_f[1:-1, x_probe])
    flux_err = abs(flux_out - flux_in) / abs(flux_in)

    print(f"   -> Integro={integro} | Erro Caudal (Inlet vs Outlet) = {flux_err * 100:.4f}%")
    print(f"   -> Erro L2 Perfil Desenvolvido Zou-He = {err_l2 * 100:.4f}%")
    return integro and (err_l2 < 0.015) and (flux_err < 0.01)



# ==============================================================================
# TESTE 6: FRONTEIRA CURVA DE BOUZIDI (q_frac != 0.5) & ESCOAMENTO NO CILINDRO
# ==============================================================================
def test_6_bouzidi_curved_cylinder():
    print("\n[TESTE 6/7] A validar Interpolacao Curva de Bouzidi (Cilindro Circular Re=30)...")
    nx, ny = 240, 100
    xc, yc, radius = 60, 50, 12.5  # Raio fracionário para forçar q < 0.5 e q > 0.5 em Bouzidi!
    u0 = 0.04
    tau = 0.60  # nu = 1/30 => Re = u0 * (2*R) / nu = 0.04 * 25 / (1/30) = 30 (Estacionário c/ recirculação)
    steps = 4000
    prefix = os.path.join(OUT_DIR, "t6_bouzidi")

    cfg = {
        "simulation": {
            "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
            "export_freq": steps, "output_prefix": prefix, "rho0": 100.0, "ux0": u0, "uy0": 0.0
        },
        "collision": {"type": "TRT", "tau_BGK": tau, "tau_TRT": tau},
        "forces": {"use_forces": False},
        "boundaries": {
            "x0": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "INLET",  "vx": u0, "vy": 0.0, "rho": 100.0},
            "x1": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "OUTLET", "vx": u0, "vy": 0.0, "rho": 100.0, "sponge_width": 20, "sponge_sigma_max": 0.15},
            "y0": {"condition": "PERIODIC", "cell_type": "FLUID"},
            "y1": {"condition": "PERIODIC", "cell_type": "FLUID"}
        },
        "obstacles": [{
            "type": "PROCEDURAL", "format": "CIRCLE", "typeboundary": "SOLID",
            "offset_x": xc, "offset_y": yc, "radius": radius
        }]
    }
    out = run_lbm_case(cfg, "cfg_t6_bouzidi")
    integro = "[100% INTEGRO]" in out

    # Extrair número de Wall Links gerados pelo GeometryBuilder
    m_links = re.search(r"Wall Links:\s+(\d+)", out)
    n_links = int(m_links.group(1)) if m_links else 0

    _, ux_f, uy_f, flag_f = read_vti_binary(f"{prefix}_{steps:06d}.vti", nx, ny)

    # Verificar que dentro do cilindro SOLID a velocidade é estritamente zero
    solid_vel_max = np.max(np.hypot(ux_f[flag_f == 1], uy_f[flag_f == 1]))

    # Verificar simetria superior vs inferior (y vs 2*yc - y)
    y_top = np.arange(yc + 1, ny - 5)
    y_bot = 2 * yc - y_top
    sym_err = compute_l2_error(ux_f[y_top, :], ux_f[y_bot, :])

    print(f"   -> Wall Links Curvos Ativos (Bouzidi): {n_links} links | Integro={integro}")
    print(f"   -> Velocidade no Interior Solido:      {solid_vel_max:.2e}")
    print(f"   -> Erro de Simetria Geometrica L2:     {sym_err:.2e}")

    # Guardar gráfico do campo de velocidade + linhas de corrente em torno do cilindro
    plt.figure(figsize=(9, 4))
    speed = np.hypot(ux_f, uy_f)
    speed[flag_f == 1] = np.nan
    plt.imshow(speed, origin="lower", cmap="turbo", extent=[0, nx, 0, ny])
    plt.colorbar(label="Magnitude da Velocidade $|\\mathbf{u}|$")
    Y, X = np.mgrid[0:ny, 0:nx]
    plt.streamplot(X, Y, ux_f, uy_f, color="white", linewidth=0.7, density=1.2)
    circle = plt.Circle((xc, yc), radius, color="black", zorder=10)
    plt.gca().add_patch(circle)
    plt.title(f"Escoamento em Cilindro Curvo com Método de Bouzidi ($Re=30$, {n_links} Links)")
    plt.xlabel("$x$ (células)")
    plt.ylabel("$y$ (células)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "5_validacao_cilindro_bouzidi.png"), dpi=200)
    plt.close()

    return integro and (n_links > 50) and (solid_vel_max == 0.0) and (sym_err < 1e-5)


# ==============================================================================
# TESTE 7: EFICÁCIA DA CAMADA ESPONJA (ABSORÇÃO ACÚSTICA) & ESTABILIDADE TRT/MRT
# ==============================================================================
def test_7_sponge_absorption_and_stability():
    print("\n[TESTE 7/7] A validar Absorcao Acustica da Esponja e Estabilidade a Baixo tau...")
    nx, ny = 120, 120
    steps = 800
    export_f = 40

    # Comparar flutuação acústica na borda SEM esponja vs COM esponja após choque gravitacional
    fluctuations = {}
    plt.figure(figsize=(8, 5))

    for mode, sp_w, sp_sig in [("Sem Esponja (Reflexão Total)", 0, 0.0), ("Com Esponja Absorvente", 20, 0.25)]:
        tag = "no_sp" if sp_w == 0 else "with_sp"
        prefix = os.path.join(OUT_DIR, f"t7_{tag}")
        cfg = {
            "simulation": {
                "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
                "export_freq": export_f, "output_prefix": prefix, "rho0": 100.0, "ux0": 0.0, "uy0": 0.0
            },
            "collision": {"type": "MRT", "tau_BGK": 0.55, "tau_MRT_filename": os.path.join(OUT_DIR, "mrt_equiv.txt")},
            "forces": {"use_forces": True, "center_x": 60, "center_y": 60, "g_mass": 120.0},
            "boundaries": {
                "x0": {"condition": "PERIODIC", "rho": 100.0, "sponge_width": sp_w, "sponge_sigma_max": sp_sig},
                "x1": {"condition": "PERIODIC", "rho": 100.0, "sponge_width": sp_w, "sponge_sigma_max": sp_sig},
                "y0": {"condition": "PERIODIC", "rho": 100.0, "sponge_width": sp_w, "sponge_sigma_max": sp_sig},
                "y1": {"condition": "PERIODIC", "rho": 100.0, "sponge_width": sp_w, "sponge_sigma_max": sp_sig}
            },
            "obstacles": []
        }
        run_lbm_case(cfg, f"cfg_t7_{tag}")

        t_list = list(range(0, steps + 1, export_f))
        ek_waves = []
        for s in t_list:
            _, ux, uy, _ = read_vti_binary(f"{prefix}_{s:06d}.vti", nx, ny)
            # Energia cinética das ondas acústicas parasitas no domínio
            ek_waves.append(np.mean(ux ** 2 + uy ** 2))

        fluctuations[tag] = ek_waves[-1]
        plt.semilogy(t_list, ek_waves, "o-", lw=2, label=mode)

    attenuation = (1.0 - fluctuations["with_sp"] / fluctuations["no_sp"]) * 100.0
    print(f"   -> Atenuacao de Energia Acustica pela Esponja: {attenuation:.1f}% de reducao!")

    plt.title("Eficácia das Camadas Absorventes (Sponge Layers) no Amortecimento Acústico")
    plt.xlabel("Passos de Tempo ($t$)")
    plt.ylabel("Energia Acústica Média $\\langle |\\mathbf{u}|^2 \\rangle$ (Escala Log)")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "6_validacao_esponja_acustica.png"), dpi=200)
    plt.close()

    return attenuation > 85.0

# ==============================================================================
# TESTE 8: GEOMETRIAS VIA FICHEIRO (.PGM, .RAW) E PROCEDURAL (RECTANGLE)
# ==============================================================================
def test_8_pgm_raw_rectangle():
    print("\n[TESTE 8/8] A validar Geometrias PGM (Asa NACA0012),Ficheiro RAW e RECTANGLE...")
    nx, ny = 360, 140

    # 1. Gerar ficheiro PGM (P5 binário) com um Perfil Alar NACA 0012 inclinado (Ângulo de Ataque = 10 graus)
    pgm_w, pgm_h = 90, 50
    pgm_path = os.path.join(OUT_DIR, "airfoil_naca0012.pgm")
    # No PGM do geometry_builder.cpp: pixel < 128 é SOLID (0=preto), pixel >= 128 é FLUID (255=branco)
    pgm_img = np.full((pgm_h, pgm_w), 255, dtype=np.uint8)
    chord = 70.0
    alpha = np.deg2rad(10.0)
    cos_a, sin_a = np.cos(alpha), np.sin(alpha)

    for py in range(pgm_h):
        for px in range(pgm_w):
            # Coordenadas centradas e rodadas pelo ângulo de ataque
            dx = px - 10.0
            dy = py - (pgm_h / 2.0)
            x_c = (dx * cos_a - dy * sin_a) / chord
            y_c = (dx * sin_a + dy * cos_a) / chord
            if 0.0 <= x_c <= 1.0:
                # Equação clássica de espessura do perfil NACA 0012
                yt = 5.0 * 0.12 * (0.2969 * np.sqrt(x_c) - 0.1260 * x_c - 0.3516 * (x_c**2) + 0.2843 * (x_c**3) - 0.1015 * (x_c**4))
                if abs(y_c) <= yt:
                    pgm_img[py, px] = 0  # Sólido (Preto)

    with open(pgm_path, "wb") as f:
        f.write(f"P5\n# NACA 0012 Airfoil\n{pgm_w} {pgm_h}\n255\n".encode("ascii"))
        # O load_pgm_file inverte o eixo Y ao ler (height - 1 - y), por isso gravamos flipud
        f.write(np.flipud(pgm_img).tobytes())

    # 2. Gerar ficheiro binário RAW (Diamante / Losango Sólido 24x24)
    raw_nx, raw_ny = 24, 24
    raw_path = os.path.join(OUT_DIR, "diamond_obstacle.raw")
    raw_buf = np.zeros((raw_ny, raw_nx), dtype=np.uint8)
    for ry in range(raw_ny):
        for rx in range(raw_nx):
            if abs(rx - raw_nx // 2) + abs(ry - raw_ny // 2) <= 9:
                raw_buf[ry, rx] = 1  # No load_raw_file: != 0 é SOLID
    with open(raw_path, "wb") as f:
        f.write(raw_buf.tobytes())

    # 3. Configurar simulação juntando PGM + RAW + RECTANGLE no mesmo túnel de vento!
    u0 = 0.04
    steps = 3000
    prefix = os.path.join(OUT_DIR, "t8_geometries")

    cfg = {
        "simulation": {
            "use_custom_init": False, "nx": nx, "ny": ny, "max_steps": steps,
            "export_freq": steps, "output_prefix": prefix, "rho0": 100.0, "ux0": u0, "uy0": 0.0
        },
        "collision": {"type": "TRT", "tau_BGK": 0.58, "tau_TRT": 1.0},
        "forces": {"use_forces": False},
        "boundaries": {
            "x0": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "INLET",  "vx": u0, "vy": 0.0, "rho": 100.0},
            "x1": {"condition": "DIRICHLET", "solver": "ZOUHE", "cell_type": "OUTLET", "vx": u0, "vy": 0.0, "rho": 100.0, "sponge_width": 25, "sponge_sigma_max": 0.2},
            "y0": {"condition": "DIRICHLET", "cell_type": "SOLID"},
            "y1": {"condition": "DIRICHLET", "cell_type": "SOLID"}
        },
        "obstacles": [
            {
                "type": "FILE", "format": "PGM", "typeboundary": "SOLID",
                "filename": pgm_path, "offset_x": 45, "offset_y": 45
            },
            {
                "type": "FILE", "format": "RAW", "typeboundary": "SOLID",
                "filename": raw_path, "obj_nx": raw_nx, "obj_ny": raw_ny, "offset_x": 175, "offset_y": 80
            },
            {
                "type": "PROCEDURAL", "format": "RECTANGLE", "typeboundary": "SOLID",
                "obj_nx": 30, "obj_ny": 18, "offset_x": 175, "offset_y": 30
            }
        ]
    }

    out = run_lbm_case(cfg, "cfg_t8_geometries")
    integro = "[100% INTEGRO]" in out
    _, ux_f, uy_f, flag_f = read_vti_binary(f"{prefix}_{steps:06d}.vti", nx, ny)

    # Verificar que os 3 obstáculos foram efetivamente estampados como SOLID (flag == 1)
    pgm_solid_cells  = np.sum(flag_f[45:45+pgm_h, 45:45+pgm_w] == 1)
    raw_solid_cells  = np.sum(flag_f[80:80+raw_ny, 175:175+raw_nx] == 1)
    rect_solid_cells = np.sum(flag_f[30:30+18, 175:175+30] == 1)

    print(f"   -> Celulas SOLID lidas do PGM (Asa NACA): {pgm_solid_cells} celulas")
    print(f"   -> Celulas SOLID lidas do RAW (Losango):  {raw_solid_cells} celulas")
    print(f"   -> Celulas SOLID do RECTANGLE (30x18):    {rect_solid_cells} celulas (Esperado = {30*18})")
    print(f"   -> Simulacao Integra:                     {integro}")

    # Desenhar o Túnel de Vento com os 3 tipos de geometria!
    plt.figure(figsize=(10, 4.2))
    speed = np.hypot(ux_f, uy_f)
    speed_masked = np.where(flag_f == 1, np.nan, speed)
    plt.imshow(speed_masked, origin="lower", cmap="turbo", extent=[0, nx, 0, ny])
    plt.colorbar(label="Magnitude da Velocidade $|\\mathbf{u}|$")
    # Pintar os sólidos a cinzento escuro
    solid_overlay = np.where(flag_f == 1, 1.0, np.nan)
    plt.imshow(solid_overlay, origin="lower", cmap="gray", vmin=0, vmax=2, extent=[0, nx, 0, ny])
    Y, X = np.mgrid[0:ny, 0:nx]
    plt.streamplot(X, Y, ux_f, uy_f, color="white", linewidth=0.6, density=1.3)
    plt.text(50, 102, "Ficheiro .PGM\n(Asa NACA 0012)", color="white", fontweight="bold", fontsize=9)
    plt.text(170, 108, "Ficheiro .RAW\n(Losango)", color="white", fontweight="bold", fontsize=9)
    plt.text(170, 12, "Procedural\n(RECTANGLE)", color="white", fontweight="bold", fontsize=9)
    plt.title("Validação do Construtor de Geometrias: Ficheiro PGM + Ficheiro RAW + RECTANGLE")
    plt.xlabel("$x$ (células)")
    plt.ylabel("$y$ (células)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "7_validacao_pgm_raw_rectangle.png"), dpi=200)
    plt.close()

    return integro and (pgm_solid_cells > 400) and (raw_solid_cells > 100) and (rect_solid_cells == 30 * 18)


if __name__ == "__main__":
    print("==================================================================")
    print("   SUITE AUTOMÁTICA DE VALIDAÇÃO CIENTÍFICA LBM D2Q9 (8 TESTES)")
    print("==================================================================")
    p1 = test_1_collider_consistency()
    p2 = test_2_shear_wave_viscosity()
    p3 = test_3_poiseuille_and_mesh_convergence()
    p4 = test_4_hydrostatic_gravity_well()
    p5 = test_5_zou_he_open_boundaries()
    p6 = test_6_bouzidi_curved_cylinder()
    p7 = test_7_sponge_absorption_and_stability()
    p8 = test_8_pgm_raw_rectangle()

    print("\n==================================================================")
    print("   RESUMO FINAL DOS 8 TESTES DE VALIDAÇÃO:")
    print(f"   [1] Consistência BGK/TRT/MRT & Massa: {'[PASSOU]' if p1 else '[FALHOU]'}")
    print(f"   [2] Viscosidade (Onda de Corte):      {'[PASSOU]' if p2 else '[FALHOU]'}")
    print(f"   [3] Poiseuille & Convergência (p~2):  {'[PASSOU]' if p3 else '[FALHOU]'}")
    print(f"   [4] Poço Gravitacional (Guo):         {'[PASSOU]' if p4 else '[FALHOU]'}")
    print(f"   [5] Fronteiras Abertas (Zou-He):      {'[PASSOU]' if p5 else '[FALHOU]'}")
    print(f"   [6] Paredes Curvas (Bouzidi Re=30):   {'[PASSOU]' if p6 else '[FALHOU]'}")
    print(f"   [7] Atenuação Acústica (Esponjas):    {'[PASSOU]' if p7 else '[FALHOU]'}")
    print(f"   [8] Geometrias (PGM, RAW, RECTANGLE): {'[PASSOU]' if p8 else '[FALHOU]'}")
    print(f"   -> 7 Gráficos guardados em: ./{OUT_DIR}/")
    print("==================================================================")