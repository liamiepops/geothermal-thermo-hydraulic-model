"""Comparative siting figure: the two governing axes (depth-to-superhot set by
gradient; stress anisotropy) and where each European site lands."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from comparative_sites import SITES, classify
from site_evaluation import evaluate

results = [evaluate(s) for s in SITES]
colors = {"GO": "tab:green", "CONDITIONAL": "tab:orange",
          "NO-GO (stability)": "tab:red", "GO* (ductile-drill)": "gold"}

fig, ax = plt.subplots(figsize=(10, 7))
for s, o in zip(SITES, results):
    surv, stab, verdict = classify(o)
    x, y = s.anisotropy, o["z"] / 1000
    c = colors.get(verdict, "gray")
    ax.scatter(x, y, s=420, c=c, edgecolors="black", zorder=3)
    label = (f"{s.name.split('(')[0].split('/')[0].strip()}\n"
             f"{o['MW_prod']:.1f} MW | {o['drill']['ROP']*3600:.1f} m/hr\n{verdict}")
    dx = 0.06 if x < 2.3 else -0.06
    ha = "left" if x < 2.3 else "right"
    ax.annotate(label, (x, y), xytext=(x + dx, y), ha=ha, va="center", fontsize=8.5)

# guide regions
ax.axvspan(1.6, 3.0, alpha=0.06, color="red")
ax.text(2.3, 1.3, "high anisotropy\n(breakout-prone)", color="tab:red",
        fontsize=8, ha="center")
ax.axhspan(11, 14, alpha=0.06, color="gray")
ax.text(1.15, 12.5, "deep (>11 km):\nstress concentration\n+ campaign cost",
        color="dimgray", fontsize=8, va="center")

ax.invert_yaxis()
ax.set_xlabel("in-situ stress anisotropy  SHmax / Shmin")
ax.set_ylabel("depth to 400 C [km]  (set by geothermal gradient)")
ax.set_xlim(1.1, 2.9); ax.set_ylim(14.5, 0.5)
ax.set_title("European superhot candidates: the two axes that decide feasibility\n"
             "(green=GO, orange=conditional, red=no-go). Sweet spot = shallow-ish + low anisotropy.",
             fontsize=10)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("comparative_sites.png", dpi=130)
print("saved comparative_sites.png")
