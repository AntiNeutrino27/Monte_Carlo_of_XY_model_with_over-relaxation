import os
import sys
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20})
sys.path.append(".")
from misc import read_readme, calc_temperatures, q_l_s, calc_chi_uv

def run_U_analysis(temp, run):
    path = os.path.dirname(os.path.realpath(__file__)) + os.sep + ".." + os.sep + "simulations" + os.sep + run
    s = 50
    num_files = len(os.listdir(path + os.sep + "simulation0/temp_0"))-2
    rm  = read_readme(f"{path}/simulation0/README.txt")
    temperatures = calc_temperatures(rm.T_min, rm.T_max, rm.n_T)
    T = temperatures[temp]
    steps = np.arange(s, num_files*s, s)
    U = []
    dU = []
    U_q = []
    Q_l = []
    dQ_l = []
    dU_q = []
    i = 0

    mc_step = []

    bs = 50


    for i in tqdm(range(0, len(steps-bs), bs)):
        u_j = []
        q_l_j = []
        for sim_fldr in os.listdir(path):
            if sim_fldr.endswith("clone"):
                continue
            sim_num = int(sim_fldr.replace("simulation", ""))
            ql, u = q_l_s(path, sim_num, temp, steps[i:i+bs])
            u_j.append(u)
            q_l_j.append(ql)
        u_j = np.array(u_j)
        q_l_j = np.array(q_l_j)
        U.append(np.mean(u_j))
        dU.append(np.std(u_j)/np.sqrt(len(u_j) - 1))
        U_q.append(3/T * (np.mean(q_l_j)-1))
        Q_l.append(np.mean(q_l_j))
        dQ_l.append(np.std(q_l_j)/np.sqrt(len(q_l_j) - 1))
        dU_q.append(3/T * np.std(q_l_j)/np.sqrt(len(q_l_j) - 1))
        mc_step.append(np.mean(steps[i:i+bs]))
    # plt.errorbar(mc_step, U, yerr=dU, capsize=2, label="Up", fmt='o-')
    # plt.errorbar(mc_step, U_q, yerr=dU_q, capsize=2, label="$U(q_l)$", fmt='o-')
    # plt.errorbar(mc_step, Q_l, yerr=dQ_l, capsize=2, label="$q_l$", fmt='o')
    # plt.xlabel("Monte Carlo Steps")
    # plt.ylabel("U")
    # plt.title(f"U vs MC Steps at T={T:.3f}")
    # print(f"{T:.3f}")
    # plt.legend()
    # plt.grid()
    # plt.show()

    # plt.plot(mc_step, np.array(U)/np.array(U_q)*np.sqrt(T))
    # plt.hlines([T], mc_step[0], mc_step[-1], colors='r', linestyles='dashed')
    # plt.show()
    return T, np.mean((np.array(U)/np.array(U_q))[-3:])


def run_corr_len(temp, run):
    path = os.path.dirname(os.path.realpath(__file__)) + os.sep + ".." + os.sep + "simulations" + os.sep + run
    s = 50
    num_files = len(os.listdir(path + os.sep + "simulation0/temp_0"))-2
    rm  = read_readme(f"{path}/simulation0/README.txt")
    L = rm.L
    N = L**3
    temperatures = calc_temperatures(rm.T_min, rm.T_max, rm.n_T)
    T = temperatures[temp]
    steps = np.arange(s, num_files*s, s)

    k_min = np.array([2*np.pi/L, 0, 0])
    k_0 = np.array([0, 0, 0])
    bs  = 50

    XI = []
    mc_step = []
    dXI = []

    for i in tqdm(range(0, len(steps-bs), bs)):
        chi_min = []
        chi_0 = []
        for sim_fldr in os.listdir(path):
            if sim_fldr.endswith("clone"):
                continue
            sim_num = int(sim_fldr.replace("simulation", ""))
            chi_mu_nu_min = calc_chi_uv(path, sim_num, temp, steps[i:i+bs], k_min)
            chi_mu_nu_0 = calc_chi_uv(path, sim_num, temp, steps[i:i+bs], k_0)
            chi_min.append(np.sum(chi_mu_nu_min))
            chi_0.append(np.sum(chi_mu_nu_0))
        chi_min_mean = np.mean(chi_min)
        chi_0_mean = np.mean(chi_0)
        chi_min_std = np.std(chi_min)/np.sqrt(len(chi_min) - 1)
        chi_0_std = np.std(chi_0)/np.sqrt(len(chi_0) - 1)

        XI.append(
            1/(2*np.sin(np.linalg.norm(k_min/2))) * np.sqrt(chi_0_mean/chi_min_mean - 1)
        )
        dXI.append(
            1/(2*np.sin(np.linalg.norm(k_min/2))) * \
            1/(2*np.sqrt(chi_0_mean/chi_min_mean - 1)) * \
            np.sqrt((chi_0_std/chi_min_mean)**2 + (chi_0_mean*chi_min_std/chi_min_mean**2)**2)
        )
        mc_step.append(np.mean(steps[i:i+bs]))

    XI = np.array(XI)
    dXI = np.array(dXI)
    plt.figure()
    plt.errorbar(mc_step, XI/L, yerr=dXI/L, capsize=2, fmt='o-')
    plt.xlabel("Monte Carlo Steps")
    plt.ylabel(r"Correlation Length $\dfrac{\xi}{L}$")
    # plt.title(f"Correlation Length vs MC Steps at T={T:.3f} and L = {L:.3f}")
    plt.grid()

    plt.figure()
    plt.errorbar(mc_step, XI/L, yerr=dXI/L, capsize=2, fmt='o-')
    plt.xlabel("Monte Carlo Steps")
    plt.ylabel(r"Correlation Length $\dfrac{\xi}{L}$")
    plt.title(f"Correlation Length vs MC Steps at T={T:.3f} and L = {L:.3f}")
    plt.xscale("log")
    plt.grid()
    plt.show()

def run_corr_len_vs_T(run):
    path = os.path.dirname(os.path.realpath(__file__)) + os.sep + ".." + os.sep + "simulations" + os.sep + run
    s = 50
    num_files = len(os.listdir(path + os.sep + "simulation0/temp_0"))-2
    rm  = read_readme(f"{path}/simulation0/README.txt")
    L = rm.L
    N = L**3
    temperatures = calc_temperatures(rm.T_min, rm.T_max, rm.n_T)

    k_min = np.array([2*np.pi/L, 0, 0])
    k_0 = np.array([0, 0, 0])
    bs  = 50

    XI_T = []
    dXI_T = []
    
    steps = np.arange(s, num_files*s, s)

    for temp in tqdm(range(rm.n_T)):
        chi_min = []
        chi_0 = []
        for sim_fldr in os.listdir(path):
            if sim_fldr.endswith("clone"):
                continue
            sim_num = int(sim_fldr.replace("simulation", ""))
            chi_mu_nu_min = calc_chi_uv(path, sim_num, temp, steps[-bs:], k_min)
            chi_mu_nu_0 = calc_chi_uv(path, sim_num, temp, steps[-bs:], k_0)
            chi_min.append(np.sum(chi_mu_nu_min))
            chi_0.append(np.sum(chi_mu_nu_0))
        chi_min_mean = np.mean(chi_min)
        chi_0_mean = np.mean(chi_0)
        chi_min_std = np.std(chi_min)/np.sqrt(len(chi_min) - 1)
        chi_0_std = np.std(chi_0)/np.sqrt(len(chi_0) - 1)

        XI_T.append(
            1/(2*np.sin(np.linalg.norm(k_min/2))) * np.sqrt(chi_0_mean/chi_min_mean - 1)
        )
        dXI_T.append(
            1/(2*np.sin(np.linalg.norm(k_min/2))) * \
            1/(2*np.sqrt(chi_0_mean/chi_min_mean - 1)) * \
            np.sqrt((chi_0_std/chi_min_mean)**2 + (chi_0_mean*chi_min_std/chi_min_mean**2)**2)  
        )
    XI_T = np.array(XI_T)
    dXI_T = np.array(dXI_T)
    # plt.figure()
    plt.errorbar(temperatures, XI_T/L, yerr=dXI_T/L,
                    capsize=2, fmt='o-', label=f"L={L}")
    
        
if __name__ == "__main__":
    T, U = [], []
    for i in range(0, 27, 3):
        t, u = run_U_analysis(i, "run3")
        T.append(np.log(t))
        U.append(u)
    plt.figure()
    plt.plot(T, U, 'o-')
    plt.show()
    # run_corr_len(16, 'run3')
    # plt.figure()
    # run_corr_len_vs_T('run2')
    # run_corr_len_vs_T('run3')
    # run_corr_len_vs_T('run4')
    # run_corr_len_vs_T('run5')
    # plt.xlabel("T")
    # plt.ylabel(r"$\dfrac{\xi}{L}$")
    # plt.title(f"Correlation Length vs Temperature")
    # plt.legend()
    # plt.grid()
    # plt.show()
