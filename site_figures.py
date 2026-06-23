"""Dashboard figure for the site evaluation (Soultz / Upper Rhine Graben)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geo_constants as C
import model1_coupled as m1
from site_evaluation import SOULTZ, evaluate

o = evaluate(SOULTZ)
s = SOULTZ
z = s.target_depth

# Model 1 profiles with the layered geotherm at a production flow
r = m1.solve(m_dot=10.0, T_inj=s.T_inj, k_ins=0.02, k_rock=s.k_rock,
             geotherm=s.geotherm, target_depth=z, Q_face=30000.0, verbose=False)
zkm = r["z"] / 1000

fig, ax = plt.subplots(1, 3, figsize=(15, 6))

# panel 1: layered geotherm + fluid profiles
ax[0].plot(r["Trock"], zkm, "k--", lw=1.5, label="rock (layered geotherm)")
ax[0].plot(r["Td"], zkm, color="tab:blue", lw=2, label="downcomer (to bit)")
ax[0].plot(r["Tu"], zkm, color="tab:red", lw=2, label="annulus (return)")
ax[0].axvline(200, color="gray", ls=":", lw=1); ax[0].text(205, 1, "200C\nceiling", fontsize=7, color="gray")
ax[0].axvline(374, color="purple", ls=":", lw=1); ax[0].text(330, 9.5, "supercrit\n374C", fontsize=7, color="purple")
ax[0].invert_yaxis(); ax[0].set_xlabel("Temperature [C]"); ax[0].set_ylabel("Depth [km]")
ax[0].set_title(f"Thermal: bit {r['T_bottom_delivered']:.0f}C, "
                f"return {r['T_return_surface']:.0f}C, {r['Q_product']/1e6:.1f} MW", fontsize=9)
ax[0].legend(fontsize=7, loc="lower left"); ax[0].grid(alpha=0.3)

# panel 2: stress profile + breakout
zz = np.linspace(0, z, 50)
ax[1].plot(s.Sv_grad*zz/1e6, zz/1000, "k-", lw=2, label="Sv (vertical)")
ax[1].plot(s.SHmax_over_Sv*s.Sv_grad*zz/1e6, zz/1000, color="tab:red", lw=1.6, label="SHmax")
ax[1].plot(s.K0_min*s.Sv_grad*zz/1e6, zz/1000, color="tab:green", lw=1.6, label="Shmin (K0=0.54)")
ax[1].plot(s.rho_fluid_grad*zz/1e6, zz/1000, color="tab:blue", lw=1.6, ls="--", label="fluid (hydrostatic)")
ax[1].invert_yaxis(); ax[1].set_xlabel("Stress [MPa]"); ax[1].set_ylabel("Depth [km]")
ax[1].set_title(f"Stress: anisotropy {s.anisotropy:.2f}\nbreakout needs +{o['breakout']['overbalance_MPa']:.0f} MPa mud", fontsize=9)
ax[1].legend(fontsize=7, loc="lower left"); ax[1].grid(alpha=0.3)

# panel 3: scorecard
ax[2].axis("off")
b = o["breakout"]; d = o["drill"]
rows = [
    ("SITE", SOULTZ.name.split("/")[1].strip()[:22], "k"),
    ("Target", f"{o['T_rock']:.0f} C @ {o['z']/1000:.1f} km", "k"),
    ("", "", "k"),
    ("Tool survival", f"OK ({o['m1']['T_bottom_delivered']:.0f} C bit)", "g"),
    ("Energy (10 kg/s)", f"{o['MW_prod']:.1f} MW_th", "k"),
    ("Drillability", f"{d['regime']}, {d['ROP']*3600:.1f} m/hr, {o['rop_gain']:.1f}x", "g"),
    ("Creep closure", f"controlled ({o['creep_hot']/max(o['creep_cold'],1e-30):.0e}x)", "g"),
    ("Isotropic stab.", f"{o['grc']['reg_c']}, {o['grc']['u_cold']*1000:.1f} mm", "g"),
    ("Breakout (aniso)", f"+{b['overbalance_MPa']:.0f} MPa mud weight", "orange"),
    ("", "", "k"),
    ("VERDICT", "CONDITIONAL GO", "g"),
]
y = 0.95
for k, v, c in rows:
    ax[2].text(0.02, y, k, fontsize=9, fontweight="bold" if k in ("SITE","VERDICT") else "normal")
    ax[2].text(0.45, y, v, fontsize=9, color={"g":"green","orange":"darkorange","k":"black"}[c])
    y -= 0.085
ax[2].set_title("Scorecard", fontsize=10)

fig.suptitle("Site evaluation dashboard: Upper Rhine Graben / Soultz-sous-Forets", fontsize=12)
fig.tight_layout()
fig.savefig("site_dashboard.png", dpi=130)
print("saved site_dashboard.png")
