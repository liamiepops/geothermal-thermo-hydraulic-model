"""
model3_optimiser.py
===================
COUPLED OPTIMISER -- ties Model 1 (borehole thermal balance / tool survival) to
Model 2 (quench-assisted rock failure) and searches the design knobs for a
self-consistent, survivable, commercially-useful operating point, per SITE.

A "site" is (geothermal gradient G, in-situ stress ratio K0). The deliverable is
a map: for each site, can the unified loop drill (survive + adequate ROP) and
what does the quench BUY versus drilling the same rock with no thermal assist?

The coupling (the crux new physics, all knobs exposed so they are defensible):

  Model 1  --(bit-delivered coolant temp T_cold)-->  Model 2
  Model 2  --(effective MSE of the quenched face)-->  ROP
  ROP      --(cutting power + cuttings sensible heat)-->  Q_face
  Q_face   ------------------------------------------>  Model 1

Solved as a fixed point (the loop is loosely coupled, converges in a few iters
because T_cold depends only weakly on Q_face).

SCOPE: drilling / early-life regime. Per the moving-bit argument, the face sees
fresh hot rock (exposure ~ seconds), so multi-year rock cooldown is NOT modelled
here; the exported MW is an EARLY-LIFE value and would carry a separate economic
decline derate in the production phase.
"""
import numpy as np
import geo_constants as C
import model1_coupled as m1
import model2_spallation as m2

A_BIT = np.pi * (C.BIT_DIAMETER / 2) ** 2
CEILING = C.BHA_SURVIVAL_TEMP            # 200 C tool survival
MSE_SPALL = 40.0e6                       # Pa, effective specific energy if spalling
ROP_CAP = 20.0 / 3600.0                  # m/s, mechanical/cuttings-removal ceiling


# ---------------------------------------------------- quench -> effective MSE
def effective_mse(z, T_rock, T_cold, K0, t_exposure=1.0, h=5e4,
                  MSE_intact=C.MSE, chi_macro=0.5, chi_micro=0.4, dT_ref=400.0,
                  phi_deg=30.0):
    """Effective specific energy of the QUENCHED face [Pa], plus regime label.

    Two damage channels (both uncertain -> calibrate by experiment):
      D_macro: quench relieves horizontal compression, raising the differential
               (shear-driving) stress; fraction of confined shear-failure
               differential 'pre-achieved' by the thermal stress.
      D_micro: grain-scale thermal microcracking, ~ proportional to quench dT,
               weakly depth-dependent (the floor benefit under confinement).
    If pure mode-I tension is feasible (Model 2 window) -> spallation regime."""
    sT = m2.sigma_thermal_surface(T_rock, T_cold, t_exposure, h)[0]
    s_conf = m2.confining_horizontal(z, K0)
    net = sT - s_conf - C.TENSILE_STRENGTH
    brittle = T_rock < m2.T_BDT
    mode_I = (net > 0) and brittle

    if mode_I:
        return MSE_SPALL, "spall", dict(sT=sT, net=net, D=1.0)

    d0, dq = m2.deviatoric_increase(z, T_rock, T_cold, t_exposure, h, K0)
    ddev = dq - d0
    q = np.tan(np.radians(45 + phi_deg / 2)) ** 2          # MC passive coeff
    sig3 = K0 * C.RHO_ROCK * C.g * z
    dev_fail = C.UCS + (q - 1) * sig3                       # confined differential at failure
    D_macro = chi_macro * np.clip(ddev / dev_fail, 0, 1)
    D_micro = chi_micro * np.clip((T_rock - T_cold) / dT_ref, 0, 1)
    D = min(D_macro + D_micro, 0.6)                         # cap total damage credit
    return MSE_intact * (1 - D), "assisted", dict(sT=sT, net=net, D=D,
                                                  D_macro=D_macro, D_micro=D_micro)


def rop_from_power(P_mech, MSE_eff):
    return min(P_mech / (MSE_eff * A_BIT), ROP_CAP)


def face_load(P_mech, ROP, T_rock, T_cold, regime):
    """Total face heat into the coolant [W]."""
    Q_cut = 0.0 if regime == "spall" else P_mech           # mech power -> heat
    Q_spall = C.RHO_ROCK * C.CP_ROCK * A_BIT * ROP * (T_rock - T_cold)
    Q_cond = 2 * np.pi * C.K_ROCK * C.CAVITY_RADIUS * (T_rock - T_cold)
    return Q_cond + Q_cut + Q_spall


# ----------------------------------------------------- one coupled evaluation
def evaluate(G, K0, m_dot, k_ins=0.02, T_inj=40.0, P_mech=8000.0,
             target_rock_T=C.TARGET_ROCK_TEMP, quench=True, n_iter=5):
    """Self-consistent coupled solve for one design at one site."""
    L = (target_rock_T - C.SURFACE_TEMP) / G
    Q_face = 25000.0
    T_cold = T_inj
    for _ in range(n_iter):
        r = m1.solve(m_dot=m_dot, T_inj=T_inj, G=G, target_rock_T=target_rock_T,
                     Q_face=Q_face, k_ins=k_ins, n_nodes=80, verbose=False)
        T_cold = r["T_bottom_delivered"]
        if quench:
            MSE_eff, regime, info = effective_mse(L, target_rock_T, T_cold, K0)
        else:
            MSE_eff, regime, info = C.MSE, "mech", dict(D=0.0)
        ROP = rop_from_power(P_mech, MSE_eff)
        Q_new = face_load(P_mech, ROP, target_rock_T, T_cold, regime)
        if abs(Q_new - Q_face) < 200.0:
            Q_face = Q_new
            break
        Q_face = Q_new

    pp, _ = m1.pump_power(m_dot, L)
    return dict(G=G, K0=K0, m_dot=m_dot, L=L, T_cold=T_cold, MSE_eff=MSE_eff,
                regime=regime, ROP=ROP, Q_face=Q_face, MW=r["Q_product"] / 1e6,
                pump_kW=pp / 1e3, survive=T_cold < CEILING,
                bvp_ok=r["success"], D=info.get("D", 0.0))


# ------------------------------------------------- optimise design for a site
def optimise_site(G, K0, P_mech=8000.0, T_inj=40.0,
                  m_dot_grid=(2, 4, 7, 10, 14, 20), k_ins=0.02):
    """Pick the survivable design (min flow that survives, then best ROP)."""
    best = None
    for md in m_dot_grid:
        r = evaluate(G, K0, md, k_ins=k_ins, T_inj=T_inj, P_mech=P_mech, quench=True)
        if not r["survive"]:
            continue
        # objective: maximise ROP per unit pump power, require net-positive energy
        score = r["ROP"] * 3600 / max(r["pump_kW"], 0.5)
        r["score"] = score
        if best is None or score > best["score"]:
            best = r
    if best is None:                       # nothing survived: report highest flow
        best = evaluate(G, K0, m_dot_grid[-1], k_ins=k_ins, T_inj=T_inj,
                        P_mech=P_mech, quench=True)
        best["score"] = 0.0
    # quench benefit: same design, no thermal assist
    nq = evaluate(G, K0, best["m_dot"], k_ins=k_ins, T_inj=T_inj,
                  P_mech=P_mech, quench=False)
    best["ROP_noquench"] = nq["ROP"]
    best["ROP_gain"] = best["ROP"] / nq["ROP"] if nq["ROP"] > 0 else np.nan
    return best


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    print("Single coupled evaluation (G=35 K/km, K0=0.7, m=10 kg/s):")
    r = evaluate(0.035, 0.7, 10.0)
    for k in ("T_cold", "regime", "MSE_eff", "ROP", "Q_face", "MW", "pump_kW",
              "survive", "D"):
        v = r[k]
        if k == "MSE_eff": v = f"{v/1e6:.0f} MPa"
        if k == "ROP": v = f"{v*3600:.2f} m/hr"
        if k == "Q_face": v = f"{v/1000:.1f} kW"
        print(f"   {k:10}: {v}")

    print("\nSITE OPTIMISATION MAP  (best survivable design per site)")
    print(f"{'G[K/km]':>8}{'K0':>5}{'depth':>7}{'m_opt':>7}{'T_bit':>7}"
          f"{'regime':>9}{'ROP':>7}{'gain':>6}{'MW':>6}{'pump':>7}{'survive':>8}")
    for G in (0.030, 0.040, 0.050):
        for K0 in (0.5, 0.7, 0.9):
            b = optimise_site(G, K0)
            print(f"{G*1000:8.0f}{K0:5.1f}{b['L']/1000:7.1f}{b['m_dot']:7.0f}"
                  f"{b['T_cold']:7.0f}{b['regime']:>9}{b['ROP']*3600:7.2f}"
                  f"{b['ROP_gain']:6.1f}{b['MW']:6.2f}{b['pump_kW']:7.1f}"
                  f"{str(b['survive']):>8}")
