"""Figure for Model 2: thermal stress vs confining stress vs depth, with the
mode-I feasible window shaded, for several in-situ stress ratios K0."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geo_constants as C
from model2_spallation import (sigma_thermal_surface, confining_horizontal,
                               feasible_window, T_BDT)

G = C.GEOTHERM_GRADIENT
T_cold, t, h = 60.0, 1.0, 5e4
z = np.linspace(50, 16000, 300)
T_rock = C.SURFACE_TEMP + G * z

# thermal-stress demand line (depends on local hot-rock T via E(T),alpha(T))
sigT = np.array([sigma_thermal_surface(Tr, T_cold, t, h)[0] for Tr in T_rock])

fig, ax = plt.subplots(figsize=(8.5, 7))
ax.plot(sigT / 1e6, z / 1000, color="tab:red", lw=2.5,
        label=r"available quench tension $\sigma_T$ (hot-rock E,$\alpha$)")

for K0, c in [(0.5, "tab:green"), (0.7, "tab:orange"), (0.9, "tab:blue"), (1.0, "tab:purple")]:
    req = confining_horizontal(z, K0) + C.TENSILE_STRENGTH   # demand to open mode-I
    ax.plot(req / 1e6, z / 1000, color=c, lw=1.6, ls="--",
            label=fr"required: $K_0\rho g z + T_0$  (K0={K0})")
    w = feasible_window(G, T_cold, t, h, K0)
    if w:
        ax.axhspan(w[0] / 1000, w[1] / 1000, xmin=0, xmax=0.06, color=c, alpha=0.5)

# brittle-ductile cap
zbdt = (T_BDT - C.SURFACE_TEMP) / G / 1000
ax.axhline(zbdt, color="k", ls=":", lw=1.2)
ax.text(5, zbdt - 0.4, f"brittle-ductile ~{T_BDT:.0f}C ({zbdt:.1f} km)", fontsize=8)

ax.invert_yaxis()
ax.set_xlabel("Stress [MPa]")
ax.set_ylabel("Depth [km]")
ax.set_xlim(0, 450)
ax.set_title("Model 2: mode-I quench spallation feasibility vs depth\n"
             "feasible only where red (available) exceeds dashed (required) AND brittle\n"
             "(shaded bars at left = feasible window per K0)", fontsize=10)
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("model2_feasibility.png", dpi=130)
print(f"saved model2_feasibility.png  (BDT at {zbdt:.1f} km)")
