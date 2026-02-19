import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.signal import savgol_filter

import os

def calc_temperatures(T_min, T_max, n_T):
    if n_T == 1:
        return np.array([T_min])
    return np.array([T_min * math.pow(T_max/T_min, i/(n_T-1.0)) for i in range(n_T)])


if __name__ == "__main__":
    path = "../run3/simulation0"
    with open(path + os.sep + "README.txt", encoding='utf-8') as file:
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
        print(f"n_T={n_T}, T_min={T_min}, T_max={T_max}, N={N}")
    temps = calc_temperatures(T_min, T_max, n_T)
    for i in range(len(temps)):
        fldr = "temp_" + str(i)
        temp = temps[i]
        energies = np.genfromtxt(f"{path}/{fldr}/energy.dat")
        # energies = np.genfromtxt(f"{path}/{fldr}/theta.dat") # theta
        size = energies.shape[0]//100
        smoothed = savgol_filter(energies[:, 1], window_length=size, polyorder=2) #* 180/np.pi  # convert to degrees
        plt.plot((energies[:, 0])[size:-size], smoothed[size:-size], label=f"T={temp:.3f}")
        if i == 0:
            plt.title("Energy vs Monte Carlo Steps for Different Temperatures")
            plt.xlabel("Monte Carlo Steps")
            plt.ylabel("Energy")

    
    plt.xticks(np.arange(0, 1.1*N, 0.1*N))
    plt.grid()
    plt.legend()
    plt.show()
        