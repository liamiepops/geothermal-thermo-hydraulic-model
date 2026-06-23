"""
model1_steadystate.py
=====================
MODEL 1, STAGE 1 -- single-depth, steady-state bottom-hole heat balance.

Question: at a chosen depth/rock temperature, what coolant mass flow keeps the
bottom-hole assembly (BHA) below its survival temperature, and which heat load
actually dominates?

We decompose the heat the coolant must carry at the FACE into three terms that
are usually conflated:

  Q_cond   conduction from hot rock into the cooled cavity wall.
           Steady hemispherical sink in a semi-infinite medium:
               Q_cond = 2*pi*k_rock*a*(T_rock - T_coolant)
           (This is your ~470 W number. It is a LOWER bound: a freshly exposed
           face conducts transiently at higher flux, ~1/sqrt(t), decaying to
           this steady value. We report both.)

  Q_cut    mechanical work of cutting, ~all of which becomes heat at the face:
               Q_cut = MSE * A_bit * ROP
           (Zero for a pure spallation head -- replaced by jet hydraulic power.)

  Q_spall  THE TERM THE ENVELOPE CALC MISSED: to advance, you remove rock and it
           leaves at ~rock temperature; cooling that removed volume from T_rock
           to the coolant temperature dumps its SENSIBLE heat into the loop:
               Q_spall = rho_rock * cp_rock * A_bit * ROP * (T_rock - T_coolant)
           This is the "drilling heat = product heat" term made concrete: the
           cuttings/spall ARE harvested geothermal heat. It scales with ROP, so
           fast drilling and big face-heat-load are the same thing.

The coolant carries Q_total = Q_cond + Q_cut + Q_spall with a temperature rise
    dT_coolant = Q_total / (m_dot * cp_water).
We also estimate the metal-to-fluid film drop to check the BHA truly tracks the
coolant temperature, not the rock.
"""
import numpy as np
import geo_constants as C
from water_props import water

KELVIN = 273.15


def depth_for_rock_temp(T_rock_C, G=C.GEOTHERM_GRADIENT, T_surf=C.SURFACE_TEMP):
    return (T_rock_C - T_surf) / G  # metres


def heat_loads(T_rock_C, T_coolant_C, a=C.CAVITY_RADIUS,
               k_rock=C.K_ROCK, A_bit=None, MSE=C.MSE, ROP=C.ROP,
               rho_rock=C.RHO_ROCK, cp_rock=C.CP_ROCK):
    if A_bit is None:
        A_bit = np.pi * (C.BIT_DIAMETER / 2) ** 2
    dT = T_rock_C - T_coolant_C

    Q_cond = 2 * np.pi * k_rock * a * dT                       # W, steady hemispherical
    Q_cut = MSE * A_bit * ROP                                   # W, cutting work
    Q_spall = rho_rock * cp_rock * A_bit * ROP * dT             # W, sensible heat of removed rock
    return dict(Q_cond=Q_cond, Q_cut=Q_cut, Q_spall=Q_spall,
                Q_total=Q_cond + Q_cut + Q_spall, A_bit=A_bit, dT=dT)


def required_flow(Q_total, dT_allow, cp_water):
    """Mass flow [kg/s] to absorb Q_total with allowed coolant temp rise dT_allow."""
    return Q_total / (cp_water * dT_allow)


def transient_conduction_flux(T_rock_C, T_coolant_C, A_bit, t_exposure,
                              k_rock=C.K_ROCK, alpha=C.ALPHA_ROCK):
    """Heat into a suddenly-cooled semi-infinite face after exposure time t.
    q'' = k*dT / sqrt(pi*alpha*t) ; integrate over face area. Upper-bound-ish
    on the conductive load for FRESH rock (vs the steady hemispherical value)."""
    dT = T_rock_C - T_coolant_C
    qpp = k_rock * dT / np.sqrt(np.pi * alpha * t_exposure)
    return qpp * A_bit


def report(T_rock_C=C.TARGET_ROCK_TEMP, T_coolant_C=60.0, dT_allow=20.0,
           ROP=C.ROP, MSE=C.MSE):
    G = C.GEOTHERM_GRADIENT
    depth = depth_for_rock_temp(T_rock_C, G)
    P = C.HYDROSTATIC_GRAD * depth                              # Pa, hydrostatic
    w = water(T_coolant_C, P)
    L = heat_loads(T_rock_C, T_coolant_C, MSE=MSE, ROP=ROP)
    m_dot = required_flow(L["Q_total"], dT_allow, w["cp"])
    # face residence/exposure time ~ time for bit to advance one cavity radius
    t_exp = C.CAVITY_RADIUS / ROP
    Q_trans = transient_conduction_flux(T_rock_C, T_coolant_C, L["A_bit"], t_exp)

    print("=" * 72)
    print(f"SINGLE-DEPTH STEADY-STATE BOTTOM-HOLE BALANCE")
    print("=" * 72)
    print(f"  Geotherm gradient G       : {G*1000:.0f} K/km")
    print(f"  Rock temperature T_rock   : {T_rock_C:.0f} C  ->  depth {depth/1000:.1f} km")
    print(f"  Hydrostatic pressure      : {P/1e6:.0f} MPa")
    print(f"  Coolant temp at face      : {T_coolant_C:.0f} C   (phase: {w['phase']})")
    print(f"  Water props @face: rho={w['rho']:.0f} cp={w['cp']:.0f} J/kgK mu={w['mu']*1e6:.0f}uPa.s")
    print(f"  Bit dia {C.BIT_DIAMETER*1000:.0f} mm (A={L['A_bit']*1e4:.0f} cm^2), ROP {ROP*3600:.2f} m/hr, MSE {MSE/1e6:.0f} MPa")
    print("-" * 72)
    print(f"  HEAT LOADS at the face:")
    print(f"    Q_cond  (wall conduction, steady) : {L['Q_cond']:9.0f} W   {100*L['Q_cond']/L['Q_total']:4.1f}%")
    print(f"    Q_cut   (cutting work -> heat)    : {L['Q_cut']:9.0f} W   {100*L['Q_cut']/L['Q_total']:4.1f}%")
    print(f"    Q_spall (sensible heat of cuttings): {L['Q_spall']:9.0f} W   {100*L['Q_spall']/L['Q_total']:4.1f}%")
    print(f"    -------------------------------------------------------")
    print(f"    Q_total                            : {L['Q_total']:9.0f} W")
    print(f"    [transient conduction on fresh face, t_exp={t_exp:.0f}s: {Q_trans:.0f} W]")
    print("-" * 72)
    print(f"  Required coolant flow for dT_allow={dT_allow:.0f} K:")
    print(f"    m_dot = {m_dot*1000:.1f} g/s   (= {m_dot:.4f} kg/s, {m_dot/w['rho']*1000:.2f} L/s)")
    # film / convective check: is the metal near the fluid temp?
    # crude: needed h so film drop < few K over bit face wetted area
    h_needed = L["Q_total"] / (L["A_bit"] * 5.0)   # for 5 K film drop
    print(f"    convective coeff for <5K film drop over face: h ~ {h_needed:.0f} W/m^2K  (jet impingement easily 10^4-10^5)")
    print("=" * 72)
    return dict(depth=depth, P=P, loads=L, m_dot=m_dot, water=w)


if __name__ == "__main__":
    # Baseline: superhot target, modest mechanical ROP
    report()
    print()
    # Sensitivity: crank ROP to the spallation-demo rate to see the face load explode
    print("### Sensitivity: ROP raised to spallation-demo rate (30 mm/min) ###")
    report(ROP=30e-3/60.0, MSE=50e6)  # spallation: low effective MSE, high ROP
