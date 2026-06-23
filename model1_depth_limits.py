"""
model1_depth_limits.py
=====================
How deep can we go, and what does depth buy -- using the models we already have?

Depth is a far better lever than diameter:
  * contact length grows LINEARLY with depth  -> heat ~ linear in L
  * rock gets hotter (T_rock = T_surf + G z)   -> bigger dT AND higher-grade heat

But three ceilings bind, and our models locate them:
  (1) BRITTLE-DUCTILE TRANSITION (Model 2's T_BDT ~400 C): hotter rock flows
      rather than fractures -> the drilling mechanism fails AND an open hole
      creeps shut. This caps the achievable ROCK TEMPERATURE, hence depth:
          z_BDT = (T_BDT - T_surf)/G
  (2) TOOL SURVIVAL (Model 1): can active cooling keep the bit < ceiling over a
      longer, hotter descent?
  (3) PARASITIC PUMP POWER (Model 1): friction grows with length; pump power
      must stay a small fraction of thermal output.

Plus a non-modelled mechanical limit: hydrostatic pressure ~10 MPa/km
(=> 150-200 MPa casing class at 15-20 km), reported for context.

Reality anchor: the Kola superdeep borehole reached 12.26 km / ~180 C (1989) and
stalled partly because hotter-than-expected rock behaved plastically -- i.e. the
PRACTICAL open-hole ceiling can arrive BELOW the 400 C textbook BDT. We treat
T_BDT as optimistic and flag it.
"""
import numpy as np
import geo_constants as C
import model1_coupled as m1

T_BDT = 400.0       # degC, optimistic open-hole/brittle ceiling (see Model 2)


def depth_to_temp(T, G, T_surf=C.SURFACE_TEMP):
    return (T - T_surf) / G


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 78)
    print("(A) The brittle-ductile depth ceiling vs geothermal gradient")
    print("    z_BDT = depth where rock hits ~400 C (hole won't drill/stay open below)")
    print("=" * 78)
    print(f"{'G[K/km]':>8}{'z_BDT[km]':>11}{'P_hydro[MPa]':>14}{'z for 374C[km]':>16}")
    for Gk in (20, 25, 35, 50, 70):
        G = Gk / 1000
        zb = depth_to_temp(T_BDT, G)
        z374 = depth_to_temp(374, G)
        print(f"{Gk:8.0f}{zb/1000:11.1f}{C.HYDROSTATIC_GRAD*zb/1e6:14.0f}{z374/1000:16.1f}")
    print("  -> LOW gradient lets you drill DEEPER before going ductile (but it's a")
    print("     longer, colder, costlier drill); HIGH gradient hits the ceiling shallow.")

    print("\n" + "=" * 78)
    print("(B) Model 1 at the BDT-limited depth for each gradient (400 C bottom)")
    print("    Same peak grade (400 C rock) -> compare what depth/length buys")
    print("=" * 78)
    print(f"{'G[K/km]':>8}{'depth[km]':>10}{'m_dot':>7}{'T_bit':>7}{'T_surf_out':>11}"
          f"{'MW_th':>7}{'pump_kW':>8}{'pump/MW%':>9}{'survive':>8}")
    for Gk in (20, 25, 35, 50, 70):
        G = Gk / 1000
        r = m1.solve(m_dot=10.0, k_ins=0.02, G=G, target_rock_T=T_BDT,
                     Q_face=25000.0, t_years=1.0, verbose=False)
        pp, _ = m1.pump_power(10.0, r["L"])
        ratio = 100 * pp / max(r["Q_product"], 1)
        print(f"{Gk:8.0f}{r['L']/1000:10.1f}{10.0:7.1f}{r['T_bottom_delivered']:7.0f}"
              f"{r['T_return_surface']:11.0f}{r['Q_product']/1e6:7.2f}{pp/1e3:8.1f}"
              f"{ratio:9.2f}{str(r['T_bottom_delivered']<200):>8}")
    print("  -> deeper (low-G) hole = more length at the SAME grade = more MW. The")
    print("     conduction limit still applies per-metre; depth adds metres + heat-grade.")

    print("\n" + "=" * 78)
    print("(C) Push a single gradient DEEPER (hypothetical: if rock stayed brittle)")
    print("    G=25 K/km -- isolates the energy benefit of depth from the BDT wall")
    print("=" * 78)
    print(f"{'bottomT':>8}{'depth[km]':>10}{'P[MPa]':>8}{'T_bit':>7}{'T_out':>7}"
          f"{'MW_th':>7}{'pump_kW':>8}{'>BDT?':>7}")
    G = 0.025
    for Tb in (300, 374, 400, 450, 500):
        r = m1.solve(m_dot=10.0, k_ins=0.02, G=G, target_rock_T=Tb,
                     Q_face=25000.0, verbose=False)
        pp, _ = m1.pump_power(10.0, r["L"])
        flag = "DUCTILE" if Tb > T_BDT else "ok"
        print(f"{Tb:8.0f}{r['L']/1000:10.1f}{C.HYDROSTATIC_GRAD*r['L']/1e6:8.0f}"
              f"{r['T_bottom_delivered']:7.0f}{r['T_return_surface']:7.0f}"
              f"{r['Q_product']/1e6:7.2f}{pp/1e3:8.1f}{flag:>7}")
    print("  -> the 450-500 C rows are physically off-limits (ductile): you cannot")
    print("     keep an open hole there. The ACTIVE-COOLING angle is the one lever that")
    print("     could push this wall deeper -- see note below.")
