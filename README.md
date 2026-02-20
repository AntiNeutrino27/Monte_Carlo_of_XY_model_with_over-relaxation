# Monte Carlo Simulation Framework

## Overview

This project implements a Monte Carlo (MC) simulation framework for particle-based systems with support for:

- Single-replica simulations
- Parallel tempering (replica exchange)
- Parallel execution over disorder realizations

The codebase separates system representation, simulation logic, and execution strategies.

There are two branches:

- `main` – single-configuration execution
- `parallel` – parallel execution over different disorder realizations of the same system configuration

---

## Project Structure

### `particle.hpp`

Represents a single particle.

Responsibilities:
- Stores particle state (position and other model-dependent properties)
- Provides access/modification of state
- Contains no simulation logic

Design principle: lightweight state container only.

---

### `environment.hpp`

Core simulation engine.

Responsibilities:
- Stores the full system (collection of particles + parameters)
- Computes total system energy
- Performs Monte Carlo updates
- Applies Metropolis acceptance rule
- Handles thermalization and measurement phases

This class performs the complete simulation for a given system.

---

### `parallel_tempering.hpp`

Implements replica exchange Monte Carlo.

Responsibilities:
- Manages `N_T` replicas at different temperatures
- Runs independent MC updates per replica
- Attempts replica exchange between neighboring temperatures
- Applies detailed-balance-consistent exchange probability

Exchange probability:

P = min(1, exp[(β_i − β_j)(E_j − E_i)])

Correct temperature spacing is required to maintain reasonable swap acceptance.

---

### `analysis/`

Contains Python scripts used for:

- Generating plots
- Post-processing simulation data
- Statistical analysis
- Disorder averaging 

All visualization and final data analysis is performed here, not in the C++ simulation code.

---

## Branches

### `main`

- Runs a single simulation for a given configuration
- No parallelization
- Intended for baseline validation and testing

---

### `parallel`

- Executes multiple simulations in parallel
- Same system configuration
- Different disorder realizations per thread/process
- Used for disorder averaging

Parallelization is performed over disorder realizations, not over replicas.

---

## Build Instructions

The simulation parameters are given in `params.hpp` file. 

To build the project in the main branch, one uses

```bash
g++ ./mcsim.cpp -std=c+17 -o mcsim_run.out
```

Consequently, in the parallel branch, one should add the `-pthread` flag. I.e.
```bash
g++ ./mcsim.cpp -std=c+17 -pthread -o mcsim_run.out
```
