"""
comparative_sites.py
===================
Run the site-evaluation tool across four European candidate provinces and
tabulate. The point: the SAME models score very different verdicts, driven by
real geotherm + stress data -- which is the pitch-ready output.

Sites & sources (gradients/stress from cited literature; deep-stress ratios are
literature-typical estimates where not directly measured at superhot depth --
flagged as the key uncertainty):
  - Upper Rhine Graben / Soultz (FR): 200C@5km, Shmin~0.54Sv, SHmax~Sv. Genter;
    Cornet/Valley "Stress State at Soultz".
  - Larderello (IT): vapour-dominated, 350C@2.2km, supercritical K-horizon 3-7km
    (= brittle-ductile transition); post-collisional EXTENSION, normal faulting.
    Bertini/Gianelli; DESCRAMBLE Venelle-2.
  - United Downs / Carnmenellis (UK): 190C@5km, ~33-35C/km radiogenic granite;
    strike-slip stress (Cornubian), high horizontal anisotropy (Pine&Batchelor,
    Rosemanowes). Reinecker/Ledingham UDDGP.
  - Pannonian Basin (HU): ~45-50C/km, heat flow 90-100 mW/m2, Miocene back-arc
    EXTENSION, thin crust. Lenkey/Horvath; Toth geothermal atlas.
"""
import numpy as np
import geo_constants as C
from site_evaluation import SiteProfile, SOULTZ, evaluate


def make_geotherm(segments):
    """segments: list of (z_top, T_top, gradient) cumulative; returns vectorized fn."""
    def gt(z):
        z = np.asarray(z, dtype=float)
        conds, vals = [], []
        for (z0, T0, g), (z1, *_ ) in zip(segments, segments[1:] + [(np.inf, 0, 0)]):
            conds.append((z >= z0) & (z < z1))
            vals.append(T0 + g * (z - z0))
        return np.select(conds, vals, default=vals[-1])
    return gt


LARDERELLO = SiteProfile(
    name="Larderello (Italy)",
    geotherm=make_geotherm([(0, 15, 0.152), (2200, 350, 0.077)]),
    target_depth=2850.0, target_T=400.0,            # just below the ~450C K-horizon
    Sv_grad=2600 * 9.81, K0_min=0.55, SHmax_over_Sv=0.75,   # normal-fault, low aniso
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.6, E_rock=45e9, UCS=140e6)

CORNWALL = SiteProfile(
    name="United Downs / Carnmenellis (UK)",
    geotherm=make_geotherm([(0, 15, 0.035), (5000, 190, 0.028)]),
    target_depth=12500.0, target_T=400.0,
    Sv_grad=2630 * 9.81, K0_min=0.55, SHmax_over_Sv=1.40,   # strike-slip, HIGH aniso
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=3.3, E_rock=60e9, UCS=180e6)

PANNONIAN = SiteProfile(
    name="Pannonian Basin (Hungary)",
    geotherm=make_geotherm([(0, 15, 0.04625), (4000, 200, 0.035)]),
    target_depth=9714.0, target_T=400.0,
    Sv_grad=2550 * 9.81, K0_min=0.60, SHmax_over_Sv=0.90,   # back-arc, moderate
    rho_fluid_grad=C.HYDROSTATIC_GRAD, k_rock=2.7, E_rock=50e9, UCS=160e6)

SITES = [SOULTZ, LARDERELLO, CORNWALL, PANNONIAN]


def classify(o):
    b = o["breakout"]
    surv = o["m1"]["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP
    if b["frac_limited"]:
        stab = "NO (frac-limited)"
    elif not b["breaks"]:
        stab = "OK"
    else:
        stab = f"+{b['overbalance_MPa']:.0f}MPa mud"
    bdt = o["T_rock"] >= 400  # at/over brittle-ductile edge
    if b["frac_limited"]:
        verdict = "NO-GO (stability)"
    elif bdt and o["T_rock"] >= 420:
        verdict = "GO* (ductile-drill)"
    elif b["breaks"]:
        verdict = "CONDITIONAL"
    else:
        verdict = "GO"
    return surv, stab, verdict


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    results = [evaluate(s) for s in SITES]

    hdr = (f"{'Site':<26}{'depth':>7}{'T_rock':>7}{'K0':>5}{'aniso':>6}"
           f"{'bitC':>6}{'MW':>5}{'ROP':>5}{'gain':>5}{'stability':>18}{'verdict':>20}")
    print("=" * len(hdr))
    print("COMPARATIVE SITING TABLE  (target 400 C, vacuum tubing, quench-assist)")
    print("=" * len(hdr))
    print(hdr)
    print("-" * len(hdr))
    for s, o in zip(SITES, results):
        surv, stab, verdict = classify(o)
        d = o["drill"]
        short = s.name.split("(")[0].strip()[:25]
        print(f"{short:<26}{o['z']/1000:6.1f}k{o['T_rock']:7.0f}{s.K0_min:5.2f}"
              f"{s.anisotropy:6.2f}{o['m1']['T_bottom_delivered']:6.0f}"
              f"{o['MW_prod']:5.1f}{d['ROP']*3600:5.1f}{o['rop_gain']:5.1f}"
              f"{stab:>18}{verdict:>20}")
    print("-" * len(hdr))
    print("notes: MW at 10 kg/s early-life; ROP m/hr; gain = quench ROP multiplier;")
    print("       'frac-limited' = mud weight to stop breakout would exceed Shmin")
    print("       (hydraulic-fracture the wall) -> not mud-controllable.")
    print()
    # one-line readout per site
    for s, o in zip(SITES, results):
        b = o["breakout"]
        print(f"* {s.name}:")
        print(f"    depth {o['z']/1000:.1f} km to {o['T_rock']:.0f} C; "
              f"Sv {o['Sv']/1e6:.0f} / SHmax {o['SHmax']/1e6:.0f} / Shmin {o['Shmin']/1e6:.0f} MPa; "
              f"breakout sig_th {b['sigma_theta']/1e6:.0f} vs MC {b['mc_cold']/1e6:.0f} MPa, "
              f"P_need {b['P_need']/1e6:.0f} vs Shmin {o['Shmin']/1e6:.0f} MPa")
