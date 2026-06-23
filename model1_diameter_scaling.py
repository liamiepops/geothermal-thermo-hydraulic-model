"""
model1_diameter_scaling.py
==========================
How does closed-loop thermal power scale with borehole DIAMETER?

Two things fight each other as the hole gets bigger:
  (+) wall contact area per length  ~  2*pi*r_w        (LINEAR in radius)
  (-) the rock can only SUPPLY heat by radial conduction, whose per-length
      conductance is   q' = 2*pi*k_rock*dT / ln(r_inf/r_w)   with the thermal
      penetration radius r_inf ~ 2*sqrt(alpha*t) set by DIFFUSION, not by the
      hole. So the supply scales only as 1/ln(r_inf/r_w)  (LOGARITHMIC).

The conduction supply is the binding constraint at production timescales, so the
net scaling of extractable power with radius is ~logarithmic = badly sub-linear.
We show this analytically and confirm it by re-running the full Model 1 BVP at
several borehole diameters.
"""
import numpy as np
import geo_constants as C
import model1_coupled as m1

ALPHA = C.K_ROCK / (C.RHO_ROCK * C.CP_ROCK)


def analytic_qprime(r_w, dT=200.0, t_years=1.0):
    """Quasi-steady conduction-limited heat per unit length [W/m]."""
    r_inf = r_w + 2.0 * np.sqrt(ALPHA * t_years * 3.1536e7)
    return 2 * np.pi * C.K_ROCK * dT / np.log(r_inf / r_w)


def set_geometry(scale):
    """Scale the whole coaxial geometry by `scale` (bigger borehole)."""
    m1.r_ii = C.R_INNER_PIPE_IN * scale
    m1.r_io = C.R_INNER_PIPE_OUT * scale
    m1.r_w = C.R_WELL * scale
    m1.A_pipe = np.pi * m1.r_ii ** 2
    m1.A_ann = np.pi * (m1.r_w ** 2 - m1.r_io ** 2)
    m1.Dh_pipe = 2 * m1.r_ii
    m1.Dh_ann = 2 * (m1.r_w - m1.r_io)


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    print("=" * 74)
    print("ANALYTIC conduction-limited heat per length vs borehole radius")
    print("  (r_inf set by 1-yr thermal diffusion ~ %.1f m)" %
          (2 * np.sqrt(ALPHA * 3.1536e7)))
    print("=" * 74)
    print(f"{'dia[m]':>8}{'r_w[m]':>8}{'wall area/len':>14}{'q-prime[W/m]':>14}{'vs 8.5in':>9}")
    q0 = analytic_qprime(C.R_WELL)
    for dia in (0.216, 0.3, 0.4, 0.6, 1.0, 2.0):
        r = dia / 2
        qp = analytic_qprime(r)
        print(f"{dia:8.2f}{r:8.3f}{2*np.pi*r:14.3f}{qp:14.1f}{qp/q0:9.2f}")
    print("  -> 10x diameter (0.2->2.0 m) raises conduction supply only ~%.1fx" %
          (analytic_qprime(1.0) / analytic_qprime(0.108)))

    print("\n" + "=" * 74)
    print("FULL MODEL 1 BVP at several diameters (G=35 K/km, 450C@12.4km)")
    print("=" * 74)
    base_dia = 2 * C.R_WELL
    print(f"{'dia[m]':>7}{'scale':>6}{'m_dot':>7}{'T_bit[C]':>9}{'T_surf':>8}"
          f"{'MW_th':>7}{'MW/dia':>8}{'pump_kW':>8}")
    print("--- fixed m_dot = 10 kg/s ---")
    for scale in (1.0, 1.5, 2.0, 3.0):
        set_geometry(scale)
        r = m1.solve(m_dot=10.0, k_ins=0.02, t_years=1.0, verbose=False)
        dia = base_dia * scale
        pp, _ = m1.pump_power(10.0, r["L"])
        print(f"{dia:7.2f}{scale:6.1f}{10.0:7.1f}{r['T_bottom_delivered']:9.0f}"
              f"{r['T_return_surface']:8.0f}{r['Q_product']/1e6:7.2f}"
              f"{r['Q_product']/1e6/dia:8.2f}{pp/1e3:8.1f}")
    print("--- m_dot scaled with annulus area (constant velocity) ---")
    for scale in (1.0, 1.5, 2.0, 3.0):
        set_geometry(scale)
        md = 10.0 * scale ** 2
        r = m1.solve(m_dot=md, k_ins=0.02, t_years=1.0, verbose=False)
        dia = base_dia * scale
        pp, _ = m1.pump_power(md, r["L"])
        print(f"{dia:7.2f}{scale:6.1f}{md:7.1f}{r['T_bottom_delivered']:9.0f}"
              f"{r['T_return_surface']:8.0f}{r['Q_product']/1e6:7.2f}"
              f"{r['Q_product']/1e6/dia:8.2f}{pp/1e3:8.1f}")

    print("\n" + "=" * 74)
    print("THE REAL LEVERS (linear, unlike diameter) -- for comparison")
    print("=" * 74)
    set_geometry(1.0)
    base = m1.solve(m_dot=10.0, k_ins=0.02, t_years=1.0, verbose=False)
    print(f"  baseline single 12.4 km hole, 10 kg/s : {base['Q_product']/1e6:.2f} MW_th")
    print(f"  -> 5 such wells (linear in count)      : {5*base['Q_product']/1e6:.2f} MW_th")
    print(f"  -> +contact length via laterals/EGS    : ~linear in contact area")
    print(f"     (Eavor uses long laterals; EGS opens fracture area -- both beat")
    print(f"      diameter, which is only logarithmic)")
