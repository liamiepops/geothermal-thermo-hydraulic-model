"""
model5_convergence_confinement.py
================================
Does a COLD LOAD-BEARING SHELL + fluid support keep a deep hole open even where
the elastic wall stress exceeds rock strength? Breakout != collapse: the rock
yields into a plastic annulus that converges to a finite displacement, and a
support pressure stabilises it. This is the tunnelling 'convergence-confinement'
method (Ground Reaction Curve), applied here with TEMPERATURE-DEPENDENT strength
and stiffness so the cooling enters physically.

Model (Duncan Fama 1993 / Hoek, circular opening, hydrostatic far field p0,
Mohr-Coulomb elastic-perfectly-plastic):

  k     = (1+sin phi)/(1-sin phi)
  p_cr  = (2 p0 - sigma_cm)/(1 + k)            critical support pressure
  if p_i >= p_cr   (ELASTIC):
        u_i = a (p0 - p_i)(1+nu)/E
        r_p = a
  if p_i <  p_cr   (PLASTIC annulus radius r_p, wall convergence u_i):
        r_p = a [ 2(p0(k-1)+sigma_cm) / ((1+k)((k-1)p_i+sigma_cm)) ]^(1/(k-1))
        u_i = a(1+nu)/E [ 2(1-nu)(p0-p_cr)(r_p/a)^2 - (1-2nu)(p0-p_i) ]

Cooling enters via sigma_cm(T) (strength rises as rock cools) and E(T). The
cold-rock GRC is valid out to the cooled-shell radius; we check self-consistency
(is the plastic zone r_p contained within the shell the dwell has cooled?).

KEY CORRECTION over Model 4(E): proper Mohr-Coulomb confines the wall via the
RADIAL stress (the fluid pressure), through the (1+k) factor -- so the simple
'sigma_theta > UCS' test in 4(E) was too pessimistic. Fluid support does a lot
of the work; cooling EXTENDS the stable envelope and cuts convergence.
"""
import numpy as np
import geo_constants as C

ALPHA = C.K_ROCK / (C.RHO_ROCK * C.CP_ROCK)
PHI = np.radians(35.0)
KMC = (1 + np.sin(PHI)) / (1 - np.sin(PHI))


def sigma_cm(T_C):
    """Rock-mass compressive strength [Pa], declining with T (cooling -> stronger)."""
    return C.UCS * max(1.0 - 0.0009 * max(T_C - 25.0, 0.0), 0.4)


def E_of_T(T_C):
    return max(C.E_ROCK * (1.0 - 0.0007 * max(T_C - 25.0, 0.0)), 0.3 * C.E_ROCK)


def grc(p_i, p0, T_wall, a=C.CAVITY_RADIUS, nu=C.NU_ROCK):
    """Ground reaction: returns (u_i [m] wall convergence, r_p/a, p_cr [Pa], regime)."""
    scm = sigma_cm(T_wall)
    E = E_of_T(T_wall)
    p_cr = (2 * p0 - scm) / (1 + KMC)
    if p_i >= p_cr:
        u_i = a * (p0 - p_i) * (1 + nu) / E
        return u_i, 1.0, p_cr, "elastic"
    rp_a = (2 * (p0 * (KMC - 1) + scm) /
            ((1 + KMC) * ((KMC - 1) * p_i + scm))) ** (1.0 / (KMC - 1))
    u_i = a * (1 + nu) / E * (2 * (1 - nu) * (p0 - p_cr) * rp_a ** 2
                              - (1 - 2 * nu) * (p0 - p_i))
    return u_i, rp_a, p_cr, "PLASTIC"


def dwell_to_cool(radius_m, a=C.CAVITY_RADIUS):
    """Time [s] for the conduction cold front (2 sqrt(alpha t)) to reach a radius."""
    delta = max(radius_m - a, 0.0)
    return (delta / 2.0) ** 2 / ALPHA


def fmt_t(s):
    if s < 86400: return f"{s/3600:.1f} h"
    if s < 3.156e7: return f"{s/86400:.1f} d"
    return f"{s/3.156e7:.1f} yr"


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    G = 0.035
    a = C.CAVITY_RADIUS

    print("=" * 80)
    print("(A) Is the hole even plastic? Fluid support vs critical pressure, by K0")
    print("    p_i = hydrostatic fluid column; p0 = K0 * lithostatic")
    print("=" * 80)
    print(f"{'z[km]':>6}{'T_rock':>7}{'p_i[MPa]':>9}{'K0':>5}{'p0[MPa]':>9}"
          f"{'p_cr_hot':>9}{'p_cr_cold':>10}{'regime(hot/cold)':>18}")
    for zkm in (10, 12, 14):
        z = zkm * 1000.0; T_rock = C.SURFACE_TEMP + G * z
        p_i = C.HYDROSTATIC_GRAD * z
        for K0 in (0.7, 1.0, 1.3):
            p0 = K0 * C.LITHOSTATIC_GRAD * z
            _, _, pcr_h, rh = grc(p_i, p0, T_rock, a)
            _, _, pcr_c, rc = grc(p_i, p0, 200.0, a)
            print(f"{zkm:6.0f}{T_rock:7.0f}{p_i/1e6:9.0f}{K0:5.1f}{p0/1e6:9.0f}"
                  f"{pcr_h/1e6:9.0f}{pcr_c/1e6:10.0f}{rh+'/'+rc:>18}")
    print("  -> hydrostatic fluid alone keeps K0<~0.9 ELASTIC even hot; compressional")
    print("     crust (K0>=1) goes plastic -> that is where cooling/shell/mud-weight earn")
    print("     their keep. (Corrects Model 4E, which ignored fluid confinement.)")

    print("\n" + "=" * 80)
    print("(B) Compressional case (K0=1.3, z=12 km): convergence HOT vs COOLED")
    print("    vs support pressure p_i  (the Ground Reaction Curve)")
    print("=" * 80)
    z = 12000.0; T_rock = C.SURFACE_TEMP + G * z
    p0 = 1.3 * C.LITHOSTATIC_GRAD * z
    print(f"  p0={p0/1e6:.0f} MPa, T_rock={T_rock:.0f} C")
    print(f"{'p_i[MPa]':>9}{'  | HOT: u[mm] r_p/a regime':<30}{'| COOLED200: u[mm] r_p/a regime':<32}")
    for pi_frac in (0.30, 0.40, 0.50, 0.60, 0.70):
        p_i = pi_frac * C.LITHOSTATIC_GRAD * z
        uh, rph, _, regh = grc(p_i, p0, T_rock, a)
        uc, rpc, _, regc = grc(p_i, p0, 200.0, a)
        print(f"{p_i/1e6:9.0f}  | {uh*1000:7.1f} {rph:5.2f} {regh:<8}"
              f"| {uc*1000:7.1f} {rpc:5.2f} {regc:<8}")
    print("  -> cooling raises strength+stiffness: smaller plastic zone, less convergence,")
    print("     and pushes the elastic threshold to lower support pressure.")

    print("\n" + "=" * 80)
    print("(C) The cold-shell ratchet: dwell to cool out to the plastic radius r_p")
    print("    (so the strong cold rock contains the yielded zone) -- K0=1.3, z=12km")
    print("=" * 80)
    print(f"{'p_i[MPa]':>9}{'r_p (hot)[m]':>13}{'shell needed[m]':>16}{'dwell':>10}"
          f"{'u_cooled[mm]':>13}")
    for pi_frac in (0.40, 0.50, 0.60):
        p_i = pi_frac * C.LITHOSTATIC_GRAD * z
        uh, rph, _, _ = grc(p_i, p0, T_rock, a)
        uc, rpc, _, regc = grc(p_i, p0, 200.0, a)
        r_p_m = rph * a
        shell = (rpc * a) if regc == "PLASTIC" else (rph * a)  # cool out to contain yield
        t = dwell_to_cool(shell)
        print(f"{p_i/1e6:9.0f}{r_p_m:13.2f}{shell:16.2f}{fmt_t(t):>10}{uc*1000:13.1f}")
    print("  -> plastic zones are sub-metre; cooling out to contain them takes hours-days,")
    print("     well within a drilling dwell. The cold shell as a contained, self-")
    print("     supporting yielded annulus is GEOMECHANICALLY CREDIBLE here.")

    print("\n" + "=" * 80)
    print("(D) Verdict band: deepest stable depth vs stress regime (hydrostatic fluid)")
    print("=" * 80)
    print(f"{'K0':>5}{'cooled? ':>9}{'max stable depth (u<1% radius)':>34}")
    for K0 in (0.7, 1.0, 1.3):
        for Tlab, Tw in (("hot", None), ("cooled200", 200.0)):
            zmax = 0
            for zkm in np.arange(6, 25, 0.5):
                z = zkm * 1000.0
                Tr = C.SURFACE_TEMP + G * z
                Tw_use = Tr if Tw is None else Tw
                p_i = C.HYDROSTATIC_GRAD * z
                p0 = K0 * C.LITHOSTATIC_GRAD * z
                u, rp, _, _ = grc(p_i, p0, Tw_use, a)
                if u < 0.01 * a:
                    zmax = zkm
                else:
                    break
            print(f"{K0:5.1f}{Tlab:>9}{zmax:>30.1f} km")

    print("\n" + "=" * 80)
    print("(E) The honest hard part: STRESS ANISOTROPY drives localized breakout")
    print("    Vertical hole: max hoop stress at the min-horizontal azimuth:")
    print("       sigma_theta = 3*S_H - S_h - P_i   (S_H=max, S_h=min horizontal)")
    print("    Mohr-Coulomb breakout when sigma_theta > k*P_i + sigma_cm.")
    print("    Aniso = S_H/S_h. z=12 km, S_h = K0*sigma_v, hydrostatic P_i.")
    print("=" * 80)
    z = 12000.0; T_rock = C.SURFACE_TEMP + G * z
    P_i = C.HYDROSTATIC_GRAD * z
    sig_v = C.LITHOSTATIC_GRAD * z
    print(f"{'K0':>5}{'Aniso':>7}{'sig_theta[MPa]':>15}{'MC limit hot':>14}"
          f"{'MC limit cold':>15}{'breakout?':>20}")
    for K0 in (0.7, 1.0):
        S_h = K0 * sig_v
        for aniso in (1.0, 1.3, 1.6):
            S_H = aniso * S_h
            sth = 3 * S_H - S_h - P_i
            lim_hot = KMC * P_i + sigma_cm(T_rock)
            lim_cold = KMC * P_i + sigma_cm(200.0)
            v = ("breaks hot&cold" if sth > lim_cold else
                 "breaks hot only" if sth > lim_hot else "stable")
            print(f"{K0:5.1f}{aniso:7.1f}{sth/1e6:15.0f}{lim_hot/1e6:14.0f}"
                  f"{lim_cold/1e6:15.0f}{v:>20}")
    print("  -> anisotropy is the real breakout driver: at Aniso~1.6 the wall breaks out")
    print("     even cold -> needs higher mud weight or accepts a breakout-widened hole.")
    print("     Cooling buys one strength tier (a 'breaks hot only' band becomes stable);")
    print("     it helps but does NOT license arbitrary anisotropy. The defensible claim:")
    print("     in moderate-anisotropy, low-K0 crust, cooling + fluid support keep a")
    print("     stable hole several km below the conventional ductile limit.")
