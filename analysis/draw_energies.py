import numpy as np
import matplotlib.pyplot as plt
import math

from scipy.signal import savgol_filter

import os

def calc_temperatures(T_min, T_max, n_T):
    return np.array([T_min * math.pow(T_max/T_min, i/(n_T-1.0)) for i in range(n_T)])


if __name__ == "__main__":
    path = "../simulations/simulation1"
    temps = calc_temperatures(0.265, 0.45, 35)
    for i in range(len(temps)):
        fldr = "temp_" + str(i)
        temp = temps[i]
        energies = np.genfromtxt(f"{path}/{fldr}/energy.dat")
        smoothed = savgol_filter(energies[:, 1], window_length=20, polyorder=3)
        plt.plot(energies[:, 0], smoothed, label=f"T={temp:.3f}")

    plt.legend()
    plt.show()
        