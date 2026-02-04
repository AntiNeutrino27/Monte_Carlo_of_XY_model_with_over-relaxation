import numpy as np
import math
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

class README:
    def __init__(self, n_T, T_min, T_max, steps, L):
        self.n_T = n_T
        self.T_min = T_min
        self.T_max = T_max
        self.steps = steps
        self.L = L

    def __str__(self):
        return f"n_T={self.n_T}, T_min={self.T_min}, T_max={self.T_max}, steps={self.steps}, L={self.L}"


def read_readme(filepath):
    with open(filepath, encoding='utf-8') as file:
        rm = file.read()
        for line in rm.split("\n"):
            if line.startswith("Number of temperatures:"):
                n_T = int(line.split()[-1])
            if line.startswith("Temperature range:"):
                parts = line.split()
                T_min = float(parts[-3])
                T_max = float(parts[-1])
            if line.startswith("Total steps:"):
                N = int(line.split()[-1])
            if line.startswith("Box size: "):
                L = int(line.split()[-1])
        return README(n_T, T_min, T_max, N, L)
    

def read_J(filepath):
    J = []
    with open(filepath, encoding='utf-8') as file:
        lines = file.readlines()
        for line in lines:
            if line is None:
                continue
            J.append(np.array([float(x) for x in line[1:-2].split(',')]))
    return np.array(J)

def calc_temperatures(T_min, T_max, n_T):
    if n_T == 1:
        return np.array([T_min])
    return np.array([T_min * math.pow(T_max/T_min, i/(n_T-1.0)) for i in range(n_T)])

def d3_d1(ix, iy, iz, L):
    ix%=L
    iy%=L
    iz%=L

    return ix + iy*L + iz*L*L

def d1_d3(i, L):
    N = L**3
    i %= N

    iz = i // (L*L)
    iy = (i - iz*L*L) // L
    ix = i - iz*L*L - iy*L

    return np.array([ix, iy, iz])


def q_uv(path_to_simulation, sim_num, temp_num, step, k):
    """
    Docstring for q_uv
    
    :param path_to_simulation: Folder containing `simulation_i` files
    :param sim_num: index of the simulation
    :param temp_num: index of the temperature
    :param step: Monte Carlo step number
    :param k: wave vector 
    """

    rm = read_readme(f"{path_to_simulation}/simulation{sim_num}/README.txt")
    L = rm.L
    N = L**3
    path1 = f"{path_to_simulation}/simulation{sim_num}/temp_{temp_num}"
    path2 = f"{path_to_simulation}/simulation{sim_num}_clone/temp_{temp_num}"

    data1 = np.genfromtxt(f"{path1}/cnf-{step}.mgl")
    data2 = np.genfromtxt(f"{path2}/cnf-{step}.mgl")

    assert data1.shape[0] == N

    R = d1_d3(np.arange(N), L)

    phase = np.sum(k[:, np.newaxis] * R, axis=0)[:, np.newaxis]

    q_uv = data1.T @ (data2 * np.exp(1j * phase)) / np.sqrt(N)
    return q_uv

def calc_chi_uv(path_to_simulation, sim_num, temp_num, steps, k):
    q2 = []
    for step in steps:
        q = q_uv(path_to_simulation, sim_num, temp_num, step, k)
        q2.append(np.abs(q)**2)
    return np.mean(q2, axis=0)


def q_l_s(path_to_simulation, sim_num, temp_num, steps):
    rm = read_readme(f"{path_to_simulation}/simulation{sim_num}/README.txt")
    J = read_J(f"{path_to_simulation}/simulation{sim_num}/J.txt")
    L = rm.L
    N = L**3
    Nb = 3*N
    
    path = f"{path_to_simulation}/simulation{sim_num}/temp_{temp_num}"
    QL = 0
    QS = 0
    Up = 0
    Ub = 0
    for step in steps:
        data = np.genfromtxt(f"{path}/cnf-{step}.mgl")
        ql = 0
        qs = 0
        u = 0 
        for index in range(N):
            S0 = data[index]
            ix, iy, iz = d1_d3(index, L)

            sxn = data[d3_d1(ix-1, iy, iz, L)]
            sxp = data[d3_d1(ix+1, iy, iz, L)]
            Jxn = J[d3_d1(ix-1, iy, iz, L)][0]
            Jxp = J[index][0]  

            syn = data[d3_d1(ix, iy-1, iz, L)]
            syp = data[d3_d1(ix, iy+1, iz, L)]
            Jyn = J[d3_d1(ix, iy-1, iz, L)][1]
            Jyp = J[index][1]

            szn = data[d3_d1(ix, iy, iz-1, L)]
            szp = data[d3_d1(ix, iy, iz+1, L)]
            Jzn = J[d3_d1(ix, iy, iz-1, L)][2]
            Jzp = J[index][2]

            ql+= (S0@sxn + S0@sxp + S0@syn + S0@syp + S0@szn + S0@szp)
            qs+= ((S0@sxn)**2 + (S0@sxp)**2 + (S0@syn)**2 + (S0@syp)**2 + (S0@szn)**2 + (S0@szp)**2)

            H = Jxn*sxn + Jxp*sxp + Jyn*syn + Jyp*syp + Jzn*szn + Jzp*szp
            u += -S0@H
        QL += (((ql)/len(steps))**2)/Nb/2
        QS += ((qs/2)/len(steps))/Nb
        Up += u/2/len(steps)/N
        Ub += u/2/len(steps)/Nb

    return QL, QS, Up, Ub

if __name__ == "__main__":
    J = read_J(os.path.dirname(os.path.realpath(__file__)) + os.sep + ".." + os.sep + "simulations/run3" + os.sep + "simulation6/J.txt")
    print(np.var(J.flatten()))
