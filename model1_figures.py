"""Depth-profile figures for Model 1 (coupled counterflow borehole)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geo_constants as C
from model1_coupled import solve

cases = [
    dict(label="m=2 kg/s, k_ins=0.1 (baseline)", m_dot=2.0, k_ins=0.10, color="tab:red"),
    dict(label="m=2 kg/s, k_ins=0.02 (vacuum tubing)", m_dot=2.0, k_ins=0.02, color="tab:orange"),
    dict(label="m=10 kg/s, k_ins=0.02", m_dot=10.0, k_ins=0.02, color="tab:green"),
]

fig, axes = plt.subplots(1, len(cases), figsize=(15, 5.2), sharey=True)
for ax, cs in zip(axes, cases):
    r = solve(m_dot=cs["m_dot"], k_ins=cs["k_ins"], verbose=False)
    zkm = r["z"] / 1000
    ax.plot(r["Trock"], zkm, "k--", lw=1.3, label="rock T_rock(z)")
    ax.plot(r["Td"], zkm, color="tab:blue", lw=2, label="downcomer T_d (to bit)")
    ax.plot(r["Tu"], zkm, color="tab:red", lw=2, label="annulus T_u (return)")
    ax.axvline(C.BHA_SURVIVAL_TEMP, color="gray", ls=":", lw=1)
    ax.text(C.BHA_SURVIVAL_TEMP + 4, 1.0, "200C\nceiling", fontsize=8, color="gray")
    ax.invert_yaxis()
    ax.set_xlabel("Temperature [C]")
    ax.set_title(cs["label"] +
                 f"\nbit {r['T_bottom_delivered']:.0f}C | "
                 f"surf {r['T_return_surface']:.0f}C | "
                 f"{r['Q_product']/1e6:.2f} MW", fontsize=9)
    ax.grid(alpha=0.3)
axes[0].set_ylabel("Depth [km]")
axes[0].legend(fontsize=8, loc="lower left")
fig.suptitle("Model 1: coupled counterflow borehole temperature profiles "
             "(G=35 K/km, 450 C at 12.4 km)", fontsize=11)
fig.tight_layout()
fig.savefig("model1_profiles.png", dpi=130)
print("saved model1_profiles.png")
