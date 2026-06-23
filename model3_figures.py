"""
Siting / operating map from the coupled optimiser.
Axes: target rock temperature (-> depth via G) x in-situ stress ratio K0.
Colour: ROP gain from quench vs no-quench. Markers: spallation regime.
Hatching: design does NOT survive (bit > ceiling).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geo_constants as C
from model3_optimiser import evaluate

G = 0.040                      # 40 K/km province
m_dot, k_ins = 4.0, 0.02
Ttargets = np.array([275, 325, 375, 425, 475])     # degC rock target
K0s = np.array([0.4, 0.6, 0.8, 1.0])

gain = np.full((len(Ttargets), len(K0s)), np.nan)
spall = np.zeros_like(gain, dtype=bool)
survive = np.zeros_like(gain, dtype=bool)

# no-quench ROP depends only on target T (not K0) -> cache per row
for i, Tt in enumerate(Ttargets):
    nq = evaluate(G, 0.7, m_dot, k_ins=k_ins, P_mech=8000.0, quench=False,
                  target_rock_T=Tt, n_iter=3)
    for j, K0 in enumerate(K0s):
        r = evaluate(G, K0, m_dot, k_ins=k_ins, P_mech=8000.0, quench=True,
                     target_rock_T=Tt, n_iter=3)
        gain[i, j] = r["ROP"] / nq["ROP"] if nq["ROP"] > 0 else np.nan
        spall[i, j] = (r["regime"] == "spall")
        survive[i, j] = r["survive"]
        print(f"T={Tt} K0={K0} depth={r['L']/1000:.1f}km regime={r['regime']:>9} "
              f"gain={gain[i,j]:.2f} survive={r['survive']}")

fig, ax = plt.subplots(figsize=(8.5, 6))
im = ax.imshow(gain, origin="lower", aspect="auto", cmap="viridis",
               extent=[K0s[0]-0.1, K0s[-1]+0.1, Ttargets[0]-25, Ttargets[-1]+25])
cb = fig.colorbar(im, ax=ax); cb.set_label("ROP gain from quench  (x)")

for i, Tt in enumerate(Ttargets):
    for j, K0 in enumerate(K0s):
        txt = f"{gain[i,j]:.1f}x"
        if spall[i, j]:
            txt += "\nSPALL"
        ax.text(K0, Tt, txt, ha="center", va="center", fontsize=8,
                color="white" if gain[i, j] < gain[np.isfinite(gain)].mean() else "black")
        if not survive[i, j]:
            ax.add_patch(plt.Rectangle((K0-0.1, Tt-25), 0.2, 50, fill=False,
                                       hatch="xxx", edgecolor="red", lw=0))

# brittle-ductile reference
ax.axhline(400, color="white", ls=":", lw=1.3)
ax.text(0.42, 405, "brittle-ductile ~400C (no spallation above)", color="white", fontsize=8)
ax.set_xlabel("in-situ stress ratio  K0  (low = extensional)")
ax.set_ylabel("target rock temperature [C]  (depth = (T-15)/G)")
ax.set_title(f"Coupled optimiser siting map (G={G*1000:.0f} K/km, m={m_dot} kg/s, "
             f"vacuum tubing)\nSPALL = pure quench spallation available; "
             f"red hatch = tool does not survive", fontsize=9)
fig.tight_layout()
fig.savefig("model3_siting_map.png", dpi=130)
print("saved model3_siting_map.png")
