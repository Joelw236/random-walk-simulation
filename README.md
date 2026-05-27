# Random Walk Simulation

A small project for simulating and analyzing a one-dimensional random walk.

The main focus is the connection between theory and numerical simulation: drift, variance, quadratic variation, and the distribution of terminal values are studied experimentally and then visualized in Python.

## Overview

The project consists of two parts:

* **C++** for simulating many random walk trajectories and exporting the results to CSV files
* **Python** for data analysis and visualization

The basic idea is simple: simulate many trajectories, compute empirical quantities, and compare them with their theoretical predictions.

## Topics Covered

* individual random walk trajectories
* empirical drift over time
* discrete quadratic variation
* empirical variance over time
* distribution of terminal values
* normalized terminal values compared to the standard normal distribution
* influence of different step probabilities `p`

## Files

### `main.cpp`

Generates the simulation data and writes the results to the following CSV files:

* `data.csv`
* `random_walk_stats.csv`
* `terminal_distribution.csv`
* `variance_convergence.csv`
* `parameter_study.csv`

### `plots.py`

Reads the CSV files and creates the corresponding plots.

## Compilation and Execution

### 1. Compile the C++ program

```bash
g++ main.cpp -O2 -o main.exe
```

### 2. Run the simulation

```bash
./main.exe
```

On Windows:

```bash
.\main.exe
```

### 3. Run the Python analysis

```bash
python plots.py
```

## Mathematical Background

The random walk starts at

$$
S_0 = 0
$$

and at each step moves either one unit upward or downward.

For the symmetric case:

$$
\mathbb{P}(\text{step upward}) = \mathbb{P}(\text{step downward}) = \tfrac{1}{2}
$$

the following classical results hold:

* $\mathbb{E}[S_t] = 0$
* $\mathrm{Var}(S_t) = t$
* the normalized terminal values $S_n / \sqrt{n}$ converge to a normal distribution as $n$ becomes large

For the more general case with downward step probability `p`, the expectation of a single increment is

$$
\mathbb{E}[X_{t+1} - X_t] = 1 - 2p
$$

which implies the expected terminal value

$$
\mathbb{E}[S_n] = n(1 - 2p)
$$

## Output

The Python scripts generate plots for:

* multiple sample trajectories
* empirical and theoretical drift
* quadratic variation
* variance growth
* distribution of terminal values
* normalized terminal values
* parameter comparison for different values of `p`

## Notes

* The simulation is intentionally kept simple and is well suited for illustrating basic stochastic concepts.
* A fixed random seed can be used in the C++ code for reproducible results.
* All plots and CSV files are generated in the current working directory.

## Project Structure

```text
.
├── main.cpp
├── plots.py
├── data.csv
├── random_walk_stats.csv
├── terminal_distribution.csv
├── variance_convergence.csv
├── parameter_study.csv
└── README.md
```

## Summary

A compact project for exploring random walks through simulation, statistics, and visualization.