import numpy as np
import matplotlib.pyplot as plt


#  Trajectories + quantile band + theory band


data = np.loadtxt("data.csv", delimiter=",")
m, T = data.shape
n = T - 1

t = np.arange(T)

q_low = np.quantile(data, 0.025, axis=0)
q_high = np.quantile(data, 0.975, axis=0)

z = 1.96
theo_low = -z * np.sqrt(t)
theo_high = z * np.sqrt(t)
theo_low[0] = 0
theo_high[0] = 0

fig, ax = plt.subplots(figsize=(12, 7))

# Walks
for walk in data:
    ax.plot(t, walk, alpha=0.05, linewidth=1)

# quantile band
ax.fill_between(t, q_low, q_high, color="gray", alpha=0.2)
ax.plot(t, q_low, color="gray", linewidth=3)
ax.plot(t, q_high, color="gray", linewidth=3)

# Theorie
ax.plot(t, theo_low, color="black", linewidth=2, label=r"Theorie: $\pm 1.96\sqrt{t}$")
ax.plot(t, theo_high, color="black", linewidth=2)

ax.axhline(0, linewidth=1)
ax.set_xlabel("t")
ax.set_ylabel("S_t")
ax.grid(True, alpha=0.8)
ax.legend()
plt.title("Random Walk Trajectories")
plt.show()


# Drift + quadratic variation check

stats = np.loadtxt("random_walk_stats.csv", delimiter=",", skiprows=1)

t2 = stats[:, 0]
drift = stats[:, 1]
qv = stats[:, 2]

fig2, ax2 = plt.subplots(figsize=(12, 6))

ax2.plot(t2, drift, label=r"Drift $E[X_{t+1}-X_t]$")
ax2.plot(t2, qv, label=r"Quadratic Variation $E[X_{t+1}^2 - X_t^2]$")

ax2.axhline(0, color="black", linewidth=1)
ax2.axhline(1, color="red", linestyle="--", linewidth=2, label="Theorie QV = 1")

ax2.set_xlabel("t")
ax2.set_ylabel("Wert")
ax2.grid(True, alpha=0.8)
ax2.legend()
plt.title("Martingale & Quadratic Variation Check")
plt.show()


# Empirical variance over time

var_data = np.loadtxt("variance_convergence.csv", delimiter=",", skiprows=1)

t_var = var_data[:, 0]
variance = var_data[:, 1]
variance_over_t = var_data[:, 2]
theory_variance = var_data[:, 3]

fig3, ax3 = plt.subplots(figsize=(12, 6))

ax3.plot(t_var, variance, label="Empirical variance")
ax3.plot(t_var, theory_variance, label="Theory: Var(S_t)=t")

ax3.set_xlabel("t")
ax3.set_ylabel("Variance")
ax3.grid(True, alpha=0.8)
ax3.legend()
plt.title("Variance Convergence")
plt.show()

fig4, ax4 = plt.subplots(figsize=(12, 6))
ax4.plot(t_var[1:], variance_over_t[1:], label=r"$\mathrm{Var}(S_t)/t$")
ax4.axhline(1, color="red", linestyle="--", linewidth=2, label="Theorie = 1")
ax4.set_xlabel("t")
ax4.set_ylabel("Ratio")
ax4.grid(True, alpha=0.8)
ax4.legend()
plt.title("Variance normalized by time")
plt.show()


# Terminal distribution + Gaussian fit

terminal = np.loadtxt("terminal_distribution.csv", delimiter=",", skiprows=1)
terminal_pos = terminal[:, 1]
scaled_terminal = terminal[:, 2]

mu = np.mean(terminal_pos)
sigma = np.std(terminal_pos, ddof=1)

fig5, ax5 = plt.subplots(figsize=(12, 6))

ax5.hist(terminal_pos, bins=50, density=True, alpha=0.6, label="Terminal positions")

x = np.linspace(terminal_pos.min(), terminal_pos.max(), 400)
if sigma > 0:
    normal = (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    ax5.plot(x, normal, linewidth=2, label="Gaussian fit")

ax5.set_xlabel("S_n")
ax5.set_ylabel("Density")
ax5.grid(True, alpha=0.8)
ax5.legend()
plt.title("Terminal Distribution")
plt.show()

fig6, ax6 = plt.subplots(figsize=(12, 6))
ax6.hist(scaled_terminal, bins=50, density=True, alpha=0.6, label=r"$S_n/\sqrt{n}$")

x2 = np.linspace(scaled_terminal.min(), scaled_terminal.max(), 400)
normal2 = (1.0 / np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * x2 ** 2)
ax6.plot(x2, normal2, linewidth=2, label="Standard normal")

ax6.set_xlabel(r"$S_n/\sqrt{n}$")
ax6.set_ylabel("Density")
ax6.grid(True, alpha=0.8)
ax6.legend()
plt.title("Scaled Terminal Distribution")
plt.show()


# Parameter study (different p)


param = np.loadtxt("parameter_study.csv", delimiter=",", skiprows=1)

p_vals = param[:, 0]
emp_step = param[:, 1]
th_step = param[:, 2]
emp_term = param[:, 3]
th_term = param[:, 4]

fig7, ax7 = plt.subplots(figsize=(12, 6))
ax7.plot(p_vals, emp_step, marker="o", label="Empirical step mean")
ax7.plot(p_vals, th_step, marker="o", label="Theoretical step mean = 1 - 2p")
ax7.set_xlabel("p")
ax7.set_ylabel("Mean increment")
ax7.grid(True, alpha=0.8)
ax7.legend()
plt.title("Drift versus p")
plt.show()

fig8, ax8 = plt.subplots(figsize=(12, 6))
ax8.plot(p_vals, emp_term, marker="o", label="Empirical terminal mean")
ax8.plot(p_vals, th_term, marker="o", label="Theoretical terminal mean")
ax8.set_xlabel("p")
ax8.set_ylabel("E[S_n]")
ax8.grid(True, alpha=0.8)
ax8.legend()
plt.title("Terminal Mean versus p")
plt.show()
