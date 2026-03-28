import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ------------------------------
# DHS model
# ------------------------------
def dhs_ln_tau(F, tau0, x_ddag, delG_ddag):
    alpha = 0.5
    arg = 1 - (alpha * F * x_ddag / delG_ddag)
    arg = np.clip(arg, 1e-12, None)
    term1 = np.log(tau0)
    term2 = (1 - 1/alpha) * np.log(arg)
    term3 = 1.69 * delG_ddag * (1 - arg**(1/alpha))
    return term1 + term2 - term3

# ------------------------------
# Bell model
# ------------------------------
def bell_ln_tau(F, tau0, x_ddag):
    return np.log(tau0) - (1.69 * F * x_ddag)

# ------------------------------
# Experimental data
# ------------------------------
F = np.array([25.6, 25.24, 23.07, 21.64, 20.37])
tau = np.array([12.22,12.34,12.31,13.68,22.31])
y_obs = np.log(tau)

# Uncertainties
Sigma_F = np.array([1.13,1.43,1.27,1.67,1.54])
Sigma_tau = np.array([0.21,0.38,0.29,0.57,1.51])
Sigma_ln_tau = Sigma_tau / tau   # error propagation

# ------------------------------
# FIT DHS (all points)
# ------------------------------
p0 = [10.0, 0.3, 15.0]

popt_dhs, pcov_dhs = curve_fit(
    dhs_ln_tau, F, y_obs,
    p0=p0,
    bounds=([0,0,0],[np.inf,np.inf,np.inf]),
    max_nfev=10000
)

tau0_dhs, x_ddag_dhs, delG_ddag_dhs = popt_dhs
err_tau0_dhs, err_x_dhs, err_G_dhs = np.sqrt(np.diag(pcov_dhs))

# R²
y_fit_dhs = dhs_ln_tau(F, *popt_dhs)
ss_res = np.sum((y_obs - y_fit_dhs)**2)
ss_tot = np.sum((y_obs - np.mean(y_obs))**2)
r2_dhs = 1 - ss_res/ss_tot

# koff + error
tau0_dhs_s = tau0_dhs * 1e-9
koff_dhs = 1 / tau0_dhs_s
koff_dhs_err = koff_dhs * (err_tau0_dhs / tau0_dhs)


# ------------------------------
# FIT BELL (ONLY FIRST FOUR POINTS)
# ------------------------------
F_bell = F[:4]
y_bell = y_obs[:4]
Sigma_ln_bell = Sigma_ln_tau[:4]

p0_bell = [10.0, 0.2]
popt_bell, pcov_bell = curve_fit(
    bell_ln_tau, F_bell, y_bell,
    p0=p0_bell
)

tau0_bell, x_ddag_bell = popt_bell
err_tau0_bell, err_x_bell = np.sqrt(np.diag(pcov_bell))

# R²
y_fit_bell = bell_ln_tau(F_bell, *popt_bell)
ss_res_bell = np.sum((y_bell - y_fit_bell)**2)
ss_tot_bell = np.sum((y_bell - np.mean(y_bell))**2)
r2_bell = 1 - ss_res_bell/ss_tot_bell

# koff + error
tau0_bell_s = tau0_bell * 1e-9
koff_bell = 1 / tau0_bell_s
koff_bell_err = koff_bell * (err_tau0_bell / tau0_bell)

# ------------------------------
# OUTPUT RESULTS
# ------------------------------
print("\n==== DHS MODEL RESULTS ====")
print(f"tau0       = {tau0_dhs:.3f} ± {err_tau0_dhs:.3f} ns")
print(f"x_ddag     = {x_ddag_dhs:.3f} ± {err_x_dhs:.3f} Å")
print(f"ΔG_ddag    = {delG_ddag_dhs:.3f} ± {err_G_dhs:.3f} kcal/mol")
print(f"DHS R²     = {r2_dhs:.4f}")
print(f"koff       = {koff_dhs:.3e} ± {koff_dhs_err:.3e} s^-1")

print("\n==== BELL MODEL RESULTS (first 4 points) ====")
print(f"tau0       = {tau0_bell:.3f} ± {err_tau0_bell:.3f} ns")
print(f"x_ddag     = {x_ddag_bell:.3f} ± {err_x_bell:.3f} Å")
print(f"Bell R²    = {r2_bell:.4f}")
print(f"koff       = {koff_bell:.3e} ± {koff_bell_err:.3e} s^-1")

# ------------------------------
# FINAL PLOT WITH ERROR BARS + SHADING
# ------------------------------

F_fit = np.linspace(min(F), max(F), 200)

plt.figure(figsize=(8,6))

# Data with error bars
plt.errorbar(
    F, y_obs,
    xerr=Sigma_F,
    yerr=Sigma_ln_tau,
    fmt='o', color='blue',
    markersize=8, capsize=3
)

# DHS curve
plt.plot(F_fit, dhs_ln_tau(F_fit, *popt_dhs),
         color='red', linewidth=3,
         label=f"DHS Fit (R²={r2_dhs:.3f})")

# Bell curve (first 4 pts)
F_fit_bell = np.linspace(min(F_bell), max(F_bell), 200)
plt.plot(F_fit_bell, bell_ln_tau(F_fit_bell, *popt_bell),
         color='green', linewidth=3,
         label=f"Bell Fit (R²={r2_bell:.3f})")

# Formatting
plt.xlabel("Force (F) (kcal/mol·Å)", fontsize=24, fontweight='bold')
plt.ylabel("ln(τ(F))", fontsize=24, fontweight='bold')

plt.xticks(fontsize=22, fontweight='bold')
plt.yticks(fontsize=22, fontweight='bold')

plt.legend(fontsize=22, frameon=False)
#plt.title("DHS vs Bell Model Fits", fontsize=24, fontweight='bold')

plt.tight_layout()
plt.savefig('kinetics.png')
plt.show()
