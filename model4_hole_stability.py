"""
model4_hole_stability.py
=======================
COOLED-SKIN BOREHOLE STABILITY -- the inversion: cooling the hot ductile rock
both extracts its energy AND hardens it (raises viscosity) so the hole stays
open. "Go deeper -> dwell, cooling/extracting until a stiff shell forms ->
go deeper." This model asks: how thick a cooled shell do you need, how long does
it take to grow, and how fast can the resulting 'ratchet' advance?

Physics
-------
1. CONVECTION CHECK (Peclet number). Test the claim that the surrounding rock
   convects heat back to the cooled nodule. Pe = v*L/alpha. If Pe << 1,
   conduction dominates and convective resupply is negligible.

2. STRESS around the hole (thick-wall, Lame, far-field isotropic S, internal
   fluid pressure P_i):   shear stress  tau(r) = (S - P_i)*(a/r)^2
   Driving deviatoric stress at the wall = S - P_i = (rho_rock - rho_fluid)*g*z
   (reducible by pressurising/weighting the fluid -- a knob).

3. CREEP closure. Power-law (dislocation) creep, Arrhenius in T:
       eps_dot(r) = A * tau(r)^n * exp(-Q/(R*T(r)))
   Wall closure rate ~ (1/a) * integral_a^inf eps_dot(r) dr.
   Baseline flow law: wet quartzite (Gleason & Tullis 1995), A=1.1e-4 MPa^-n/s,
   n=4, Q=223 kJ/mol -- REPRESENTATIVE & UNCERTAIN; Q is swept.

4. COOLED TEMPERATURE PROFILE after dwell time t_d (quasi-steady log shell out
   to the conduction front r_inf = a + 2*sqrt(alpha*t_d)):
       T(r) = T_cold + (T_rock - T_cold)*ln(r/a)/ln(r_inf/a),  clipped to T_rock.
   Longer dwell -> thicker cold shell -> the hot creeping rock is pushed out to
   where stress (tau ~ 1/r^2) is weak -> closure suppressed.

All flow-law parameters are exposed; the ROBUST result is the Arrhenius
sensitivity (cooling buys orders of magnitude), not the absolute rate.
"""
import numpy as np
from scipy.integrate import quad
import geo_constants as C

R_GAS = 8.314
YEAR = 3.156e7
ALPHA = C.K_ROCK / (C.RHO_ROCK * C.CP_ROCK)

# --- flow law (wet quartzite, Gleason & Tullis 1995) -- representative ---
A_FL = 1.1e-4      # MPa^-n s^-1
N_FL = 4.0
Q_FL = 223e3       # J/mol


def driving_stress(z, fluid_grad=C.HYDROSTATIC_GRAD):
    """Far-field minus internal fluid pressure [Pa] -> deviatoric driver."""
    return (C.LITHOSTATIC_GRAD - fluid_grad) * z


def T_profile(r, a, t_dwell, T_rock, T_cold):
    """Radial temperature [degC] after dwell time; log shell to conduction front."""
    r_inf = a + 2.0 * np.sqrt(ALPHA * max(t_dwell, 1.0))
    if r >= r_inf:
        return T_rock
    frac = np.log(r / a) / np.log(r_inf / a)
    return T_cold + (T_rock - T_cold) * frac


def closure_rate(z, T_rock, T_cold, t_dwell, a=C.CAVITY_RADIUS,
                 A=A_FL, n=N_FL, Q=Q_FL, cooled=True, r_out=50.0):
    """Fractional wall closure rate [1/s]. cooled=False -> wall stays at T_rock."""
    dsig = driving_stress(z) / 1e6   # MPa
    def integrand(r):
        tau = dsig * (a / r) ** 2                      # MPa
        T = T_profile(r, a, t_dwell, T_rock, T_cold) if cooled else T_rock
        return A * tau ** n * np.exp(-Q / (R_GAS * (T + 273.15)))
    val, _ = quad(integrand, a, a + r_out, limit=100)
    return val / a


def strength(T_C, UCS0=C.UCS):
    """Unconfined compressive strength [Pa], declining with temperature
    (~ -40% by 500 C). Cold rock is stronger -> cooling helps the brittle side."""
    return UCS0 * max(1.0 - 0.0009 * max(T_C - 25.0, 0.0), 0.4)


def wall_tangential_stress(z, fluid_grad=C.HYDROSTATIC_GRAD, K0=1.0):
    """Kirsch tangential stress at a VERTICAL borehole wall. The driver is the
    HORIZONTAL far-field stress S_h = K0 * sigma_v (NOT the vertical/overburden).
    Isotropic-horizontal: sigma_theta = 2*S_h - P_i  [Pa].
    Extensional crust (low K0) => much lower breakout stress -- same settings that
    favoured spallation in Models 2-3."""
    S_h = K0 * C.LITHOSTATIC_GRAD * z
    P_i = fluid_grad * z
    return 2 * S_h - P_i


def fmt_time(seconds):
    if seconds < 3600: return f"{seconds:.0f} s"
    if seconds < 86400: return f"{seconds/3600:.1f} h"
    if seconds < YEAR: return f"{seconds/86400:.1f} d"
    return f"{seconds/YEAR:.1f} yr"


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 78)
    print("(A) CONVECTION CHECK -- does the rock flow heat back to the nodule?")
    print("=" * 78)
    for v_cm_yr, label in [(1.0, "fast lower-crust"), (10.0, "vigorous"), (0.1, "typical")]:
        v = v_cm_yr / 100 / YEAR        # m/s
        for L in (1.0, 10.0):
            Pe = v * L / ALPHA
            print(f"   v={v_cm_yr:4.1f} cm/yr ({label:16}) L={L:4.0f} m -> Pe={Pe:.1e}")
    print("   Conductive cold-front speed at 10 m shell ~ alpha/L = "
          f"{ALPHA/10*YEAR*100:.1e} cm/yr (>> convection).")
    print("   => Pe << 1: convection is thermally NEGLIGIBLE. Fresh heat comes from")
    print("      ADVANCING THE BIT into new hot rock, not from rock convecting in.")

    print("\n" + "=" * 78)
    print("(B) Creep closure: HOT wall (uncooled) vs COOLED wall (200 C, dwell)")
    print("    time-to-close 10% of radius; G=35 K/km, fluid hydrostatic")
    print("=" * 78)
    G = 0.035
    print(f"{'z[km]':>6}{'T_rock':>8}{'drive[MPa]':>11}"
          f"{'hot: t_close':>16}{'cooled200: t_close':>20}")
    for zkm in (8, 10, 11, 12, 14):
        z = zkm * 1000.0; T_rock = C.SURFACE_TEMP + G * z
        cr_hot = closure_rate(z, T_rock, 200.0, t_dwell=7*86400, cooled=False)
        cr_cold = closure_rate(z, T_rock, 200.0, t_dwell=7*86400, cooled=True)
        th = 0.1 / cr_hot if cr_hot > 0 else np.inf
        tc = 0.1 / cr_cold if cr_cold > 0 else np.inf
        print(f"{zkm:6.0f}{T_rock:8.0f}{driving_stress(z)/1e6:11.0f}"
              f"{fmt_time(th):>16}{fmt_time(tc):>20}")

    print("\n" + "=" * 78)
    print("(C) The DWELL lever: closure rate vs cooling dwell time (z=12 km, 435 C)")
    print("=" * 78)
    z = 12000.0; T_rock = C.SURFACE_TEMP + G * z
    print(f"{'dwell':>8}{'shell~[m]':>11}{'T_cold':>8}{'closure[%/yr]':>15}{'t_close10%':>14}")
    for t_d, lbl in [(3600, "1 h"), (86400, "1 d"), (7*86400, "1 wk"),
                     (30*86400, "1 mo"), (365*86400, "1 yr")]:
        shell = 2 * np.sqrt(ALPHA * t_d)
        cr = closure_rate(z, T_rock, 200.0, t_dwell=t_d, cooled=True)
        print(f"{lbl:>8}{shell:11.2f}{200.0:8.0f}{cr*YEAR*100:15.3f}"
              f"{fmt_time(0.1/cr):>14}")

    print("\n" + "=" * 78)
    print("(D) Robustness: sweep activation energy Q (the big uncertainty)")
    print("    cooled wall 200 C, z=12 km, 1-week dwell -> closure [%/yr]")
    print("=" * 78)
    for Q in (180e3, 223e3, 270e3):
        cr_hot = closure_rate(z, T_rock, 200.0, 7*86400, Q=Q, cooled=False)
        cr_cold = closure_rate(z, T_rock, 200.0, 7*86400, Q=Q, cooled=True)
        ratio = cr_hot / cr_cold if cr_cold > 0 else np.inf
        print(f"   Q={Q/1e3:.0f} kJ/mol:  hot {cr_hot*YEAR*100:10.1f} %/yr | "
              f"cooled {cr_cold*YEAR*100:.2e} %/yr | cooling factor {ratio:.1e}x")

    print("\n" + "=" * 78)
    print("(E) THE REAL DEEP PROBLEM: brittle OVER-STRESS (instantaneous), not creep")
    print("    wall tangential stress sigma_theta=2S-P_i  vs  rock strength.")
    print("    If sigma_theta > strength -> the wall YIELDS (breakouts), regardless of")
    print("    creep rate. This is what cooling only PARTLY fixes.")
    print("=" * 78)
    print("    sigma_theta with pressurised fluid (P_i=0.85*lithostatic), strength@200C cold")
    print(f"{'z[km]':>6}{'T_rock':>8}{'K0=1.0':>9}{'K0=0.7':>9}{'K0=0.5':>9}"
          f"{'str_cold':>10}{'worst-case verdict':>20}")
    s_cold = strength(200.0) / 1e6
    for zkm in (8, 10, 12, 14):
        z = zkm * 1000.0; T_rock = C.SURFACE_TEMP + G * z
        Pi = 0.85 * C.LITHOSTATIC_GRAD
        st = {k: wall_tangential_stress(z, Pi, K0=k) / 1e6 for k in (1.0, 0.7, 0.5)}
        verdict = "yields even cold" if st[1.0] > s_cold else "stable"
        print(f"{zkm:6.0f}{T_rock:8.0f}{st[1.0]:9.0f}{st[0.7]:9.0f}{st[0.5]:9.0f}"
              f"{s_cold:10.0f}{verdict:>20}")
    print("  -> Breakout severity scales with HORIZONTAL stress (K0). In compressional")
    print("     crust (K0>=1) sigma_theta >> strength even cold: fluid pressure dominates,")
    print("     cooling only helps ~1.5x. In EXTENSIONAL crust (K0~0.5) it drops near the")
    print("     strength line -- the SAME low-K0 settings that favoured spallation. Cooling")
    print("     kills time-dependent convergence; a cold load-bearing shell + fluid support")
    print("     is the credible (open) stability strategy -- not cooling alone.")
