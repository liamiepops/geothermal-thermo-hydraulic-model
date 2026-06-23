"""
model2_spallation.py
====================
MODEL 2 -- cutting-face / coolant-assisted rock failure, done HONESTLY with the
in-situ confining stress in the same frame as the thermal stress.

Sign convention: TENSION POSITIVE throughout.

Physics
-------
1. Transient thermal field. A hot rock face (T_rock) is quenched by coolant at
   T_cold via convective coefficient h. Semi-infinite solid with surface
   convection (Carslaw & Jaeger):

       theta(x,t) = erfc(eta) - exp(2*b*eta) * erfcx(eta + b)
       eta = x / (2*sqrt(alpha*t)),   b = h*sqrt(alpha*t)/k_rock   (Biot-like)
       T(x,t) = T_rock + (T_cold - T_rock) * theta(x,t)

   (erfcx = exp(z^2)*erfc(z), used for numerical stability at large argument.)

2. Thermoelastic stress. A thin quenched skin on a massive body is fully
   laterally constrained by the (unchanged) bulk, so the in-plane stress is the
   constrained thermoelastic value, tensile where cooled:

       sigma_T(x,t) = E*alpha/(1-nu) * (T_rock - T(x,t))      (>=0, tensile)

   Peak at the surface: sigma_T_surf = E*alpha/(1-nu) * (T_rock - T_cold) * theta_s.

3. In-situ stress. Vertical (overburden) and horizontal:
       sigma_v  = -rho_rock*g*z            (compressive, negative)
       sigma_h0 = -K0 * rho_rock*g*z       (compressive, negative)
   The thermal stress adds IN-PLANE (horizontal), so near the face:
       sigma_h(x) = sigma_h0 + sigma_T(x)

Failure checks (the honest part)
--------------------------------
* MODE I (tensile spallation): needs net in-plane tension above tensile
  strength T0:   sigma_h0 + sigma_T_surf > T0
  => sigma_T_surf > T0 + K0*rho*g*z. Defines a CROSSOVER DEPTH below which pure
  quench spallation cannot make net tension. This is the bench->depth gap.

* SHEAR (Mohr-Coulomb): cooling REDUCES horizontal compression while vertical
  stays, raising the differential stress (sigma1 - sigma3). We report whether
  quench pushes the stress state toward shear failure even when mode I is denied.

* GRAIN-SCALE microcracking from thermal-expansion mismatch occurs regardless of
  the macroscopic state and lowers effective strength / MSE; we treat that as a
  parameterized MSE-reduction back-coupled to Model 1 rather than a continuum
  stress, and flag it as the likely real payoff at depth.
"""
import numpy as np
from scipy.special import erfc, erfcx
from scipy.optimize import brentq
import geo_constants as C

# ---- optional temperature dependence (granite softens & expands more when hot)
def E_of_T(T_C):
    """Young's modulus declines with T. ~ -0.04 GPa/K from 25C baseline (lab data
    on granite gives ~30-50% drop by 500-600C). Returns Pa."""
    E = C.E_ROCK * (1.0 - 0.0007 * max(T_C - 25.0, 0.0))   # ~ -40% by 600C
    return max(E, 0.3 * C.E_ROCK)

def alpha_of_T(T_C):
    """Linear thermal expansion rises with T; ~8e-6 at 25C to ~11e-6 at 500C."""
    return C.ALPHA_THERMAL * (1.0 + 0.0008 * max(T_C - 25.0, 0.0))


# ---------------------------------------------------------------- thermal field
def theta_profile(x, t, alpha, h, k_rock):
    # theta = erfc(eta) - exp(2 b eta + b^2) erfc(eta+b).
    # Stable form: the exp() prefactor times erfc(eta+b) = exp(-eta^2)*erfcx(eta+b).
    eta = x / (2.0 * np.sqrt(alpha * t))
    b = h * np.sqrt(alpha * t) / k_rock
    return erfc(eta) - np.exp(-eta * eta) * erfcx(eta + b)

def theta_surface(t, alpha, h, k_rock):
    b = h * np.sqrt(alpha * t) / k_rock
    return 1.0 - erfcx(b)            # 1 - exp(b^2) erfc(b)


# --------------------------------------------------------------- thermal stress
def sigma_thermal_surface(T_rock, T_cold, t, h, T_for_props=None,
                          alpha_rock=C.ALPHA_ROCK, k_rock=C.K_ROCK):
    """Peak (surface) constrained thermoelastic tensile stress [Pa]."""
    Tp = T_rock if T_for_props is None else T_for_props
    E = E_of_T(Tp); al = alpha_of_T(Tp); nu = C.NU_ROCK
    th = theta_surface(t, alpha_rock, h, k_rock)
    return E * al / (1 - nu) * (T_rock - T_cold) * th, E, al, th


# ----------------------------------------------------------------- in-situ load
def confining_horizontal(z, K0=0.9):
    """Magnitude of in-situ horizontal compressive stress [Pa] at depth z."""
    return K0 * C.RHO_ROCK * C.g * z


# --------------------------------------------------------- failure / crossover
def net_inplane_stress(z, T_rock, T_cold, t, h, K0=0.9):
    s_T, E, al, th = sigma_thermal_surface(T_rock, T_cold, t, h)
    s_conf = confining_horizontal(z, K0)        # compressive magnitude
    return s_T - s_conf, s_T, s_conf            # net (tension +) , thermal, confining


def margin(z, G, T_cold, t, h, K0, T_surf=C.SURFACE_TEMP, T0=C.TENSILE_STRENGTH):
    """Net mode-I margin [Pa] at depth z: sigma_T - K0 rho g z - T0.
    >0 means quench can open net tensile (mode-I) cracks."""
    T_rock = T_surf + G * z
    s_T, *_ = sigma_thermal_surface(T_rock, T_cold, t, h)
    return s_T - confining_horizontal(z, K0) - T0


T_BDT = 400.0   # degC, approx brittle-ductile transition for granite; above this
                # rock flows rather than cracks, so spallation cannot work there
                # (raised somewhat by fast drilling strain rates -- a real caveat).

def feasible_window(G, T_cold, t, h, K0, zmax=20000.0, n=400, T_surf=C.SURFACE_TEMP):
    """Return (z_lo, z_hi) depth window [m] where mode-I quench spallation is
    feasible: net tension > T0 AND rock still brittle (T_rock < T_BDT)."""
    zs = np.linspace(1.0, zmax, n)
    m = np.array([margin(z, G, T_cold, t, h, K0) for z in zs])
    brittle = (T_surf + G * zs) < T_BDT
    pos = zs[(m > 0) & brittle]
    return (pos.min(), pos.max()) if pos.size else None


def competing_gradients(G, K0, T_rock_ref=450.0):
    """The headline numbers: depth-rate of thermal-stress demand vs confining
    supply [Pa/m]. If thermal slope < confining slope, depth always loses."""
    E = E_of_T(T_rock_ref); al = alpha_of_T(T_rock_ref)
    thermal_slope = E * al / (1 - C.NU_ROCK) * G          # per metre of depth
    confining_slope = K0 * C.RHO_ROCK * C.g
    return thermal_slope, confining_slope


def deviatoric_increase(z, T_rock, T_cold, t, h, K0):
    """Quench reduces horizontal compression (adds sigma_T tensile) while
    vertical stress is unchanged -> raises the differential stress that drives
    SHEAR failure. Returns (differential_no_quench, differential_with_quench) [Pa].
    Tension-positive: sigma_v=-rho g z, sigma_h0=-K0 rho g z, sigma_h=sigma_h0+sigma_T."""
    s_T, *_ = sigma_thermal_surface(T_rock, T_cold, t, h)
    sig_v = -C.RHO_ROCK * C.g * z
    sig_h0 = -K0 * C.RHO_ROCK * C.g * z
    diff0 = abs(sig_v - sig_h0)
    diffq = abs(sig_v - (sig_h0 + s_T))
    return diff0, diffq


def spallation_rate_estimate(t_cycle, alpha=C.ALPHA_ROCK):
    """Crude flake thickness ~ thermal penetration 2*sqrt(alpha*t); implied
    spallation ROP if a flake pops each cycle of duration t_cycle."""
    delta = 2.0 * np.sqrt(alpha * t_cycle)     # m
    return delta, delta / t_cycle              # m, m/s


# ------------------------------------------------------------------- reporting
def report(G=C.GEOTHERM_GRADIENT, T_cold=60.0, t=1.0, h=5e4):
    print("=" * 76)
    print("MODEL 2 -- coolant-assisted rock failure WITH in-situ confining stress")
    print("=" * 76)
    print(f"  coolant T_cold={T_cold} C | exposure t={t}s | h={h:.0e} W/m2K | "
          f"G={G*1000:.0f} K/km | T0={C.TENSILE_STRENGTH/1e6:.0f} MPa")
    print("-" * 76)
    print("  (a) Surface thermal stress saturates fast with jet-class h")
    print("      (=> heat transfer is NOT the limiter; confinement is):")
    for tt in (0.01, 0.1, 1.0):
        s, E, al, th = sigma_thermal_surface(450.0, T_cold, tt, h)
        print(f"      t={tt:5.2f}s  theta_s={th:5.3f}  sigma_T(450C)={s/1e6:5.0f} MPa "
              f"(E={E/1e9:.0f}GPa, alpha={al*1e6:.1f}e-6)")
    print("-" * 76)
    print("  (b) HEADLINE: competing depth-gradients [kPa per metre of depth]")
    print("      thermal-stress demand grows as Eα/(1-ν)·G ; confinement as K0·ρg")
    th_slope, _ = competing_gradients(G, 1.0)
    print(f"      thermal-stress slope (hot-rock E) : {th_slope/1e3:5.1f} kPa/m")
    for K0 in (0.5, 0.7, 0.9, 1.0):
        _, cs = competing_gradients(G, K0)
        verdict = "thermal keeps up" if th_slope > cs else "confinement WINS w/ depth"
        print(f"      confining slope  K0={K0:.1f}        : {cs/1e3:5.1f} kPa/m   -> {verdict}")
    print("-" * 76)
    print("  (c) Mode-I feasible depth WINDOW vs in-situ stress ratio K0:")
    print(f"      (rock must be hot enough AND not too confined; T_rock=15+G·z)")
    for K0 in (0.5, 0.7, 0.9, 1.0):
        w = feasible_window(G, T_cold, t, h, K0)
        if w is None:
            print(f"      K0={K0:.1f}:  NO window -- pure mode-I quench spallation never feasible")
        else:
            print(f"      K0={K0:.1f}:  feasible {w[0]/1000:4.1f}-{w[1]/1000:4.1f} km "
                  f"(rock {15+G*w[0]:.0f}-{15+G*w[1]:.0f} C)")
    print("-" * 76)
    print("  (d) Even where mode-I is DENIED, quench raises the differential")
    print("      (shear-driving) stress -- the likely real mechanism at depth:")
    print(f"      {'z[km]':>6}{'T_rock':>8}{'diff_noquench':>15}{'diff_quench':>13}{'increase':>10}")
    for zkm in (8, 10, 12):
        z = zkm * 1000.0; T_rock = C.SURFACE_TEMP + G * z
        d0, dq = deviatoric_increase(z, T_rock, T_cold, t, h, 0.7)
        print(f"      {zkm:6.0f}{T_rock:8.0f}{d0/1e6:12.0f} MPa{dq/1e6:10.0f} MPa{(dq-d0)/1e6:8.0f} MPa")
    print("-" * 76)
    print("  (e) Flake thickness & implied spallation ROP (if spallation occurs):")
    for tc in (0.1, 1.0, 10.0):
        d, rop = spallation_rate_estimate(tc)
        print(f"      t_cycle={tc:5.1f}s  flake~{d*1000:4.1f} mm  ROP~{rop*1000:5.2f} mm/s = {rop*3600:5.2f} m/hr")
    print("      (Holtzman 2023 bench: ~30 mm/min = 0.50 mm/s at 480C, UNCONFINED)")
    print("=" * 76)


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    report(T_cold=60.0)
