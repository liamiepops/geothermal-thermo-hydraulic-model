"""
model1_coupled.py
=================
MODEL 1, STAGE 2 -- coupled 1-D counterflow borehole heat exchanger.

Coaxial well (Eavor-like), z measured DOWNWARD from surface (z=0) to bit (z=L):
  * COLD fluid flows DOWN an insulated central pipe:   T_d(z)
  * HOT  fluid returns UP the annulus, against the rock: T_u(z)
  * rock wall temperature: T_rock(z) = T_surf + G*z

This routing is the crux of the "one loop does both" thesis: insulate the
downcomer so fluid reaches the bit COLD (tool survival + quench), while the
annulus harvests the rock's heat on the way up (the energy product).

Energy balances (per unit length), m_dot = mass flow in each leg:

  m_dot*cp_d * dT_d/dz =  UAi(z) * (T_u - T_d)
  m_dot*cp_u * dT_u/dz = -UAo(z) * (T_rock(z) - T_u) + UAi(z) * (T_u - T_d)

  UAi(z) [W/m/K] = centre<->annulus conductance per length (set by INSULATION)
  UAo(z) [W/m/K] = rock<->annulus conductance per length (rock conduction + film)

BCs:
  T_d(0) = T_inj
  T_u(L) = T_d(L) + Q_face/(m_dot*cp_L)          (face heat load at turnaround)

Water properties: IAPWS-95 via a precomputed table (water_table.py),
fully vectorized. P(z) hydrostatic (+1 atm surface).
"""
import numpy as np
from scipy.integrate import solve_bvp
import geo_constants as C
import water_table as wt

P_SURF = 1.0e5  # Pa, atmospheric at wellhead
P_OP = 5.0e6    # Pa, loop operating pressure for the surface heat-product calc
                # (keeps the hot return liquid; a pumped loop is pressurized)

# Geometry
r_ii = C.R_INNER_PIPE_IN
r_io = C.R_INNER_PIPE_OUT
r_w = C.R_WELL
A_pipe = np.pi * r_ii ** 2
A_ann = np.pi * (r_w ** 2 - r_io ** 2)
Dh_pipe = 2 * r_ii
Dh_ann = 2 * (r_w - r_io)


def P_of_z(z):
    return P_SURF + C.HYDROSTATIC_GRAD * z


def h_dittus(m_dot, area, Dh, T_C, P_Pa, n=0.4):
    """Vectorized convective coefficient [W/m^2/K]. Laminar floor Nu=4.36."""
    rho, cp, mu, k = wt.props(T_C, P_Pa)
    Gflux = m_dot / area
    Re = Gflux * Dh / mu
    Pr = cp * mu / k
    Nu = np.where(Re < 2300.0, 4.36, 0.023 * np.power(Re, 0.8) * np.power(Pr, n))
    return Nu * k / Dh


def conductances(z, T_d, T_u, m_dot, t_years, k_ins, k_rock):
    """Vectorized per-length UAi (centre<->annulus), UAo (rock<->annulus)."""
    P = P_of_z(z)
    h_d = h_dittus(m_dot, A_pipe, Dh_pipe, T_d, P)
    h_a = h_dittus(m_dot, A_ann, Dh_ann, T_u, P)

    R_in = 1.0 / (h_d * 2 * np.pi * r_ii)
    R_wall = np.log(r_io / r_ii) / (2 * np.pi * k_ins)
    R_ao = 1.0 / (h_a * 2 * np.pi * r_io)
    UAi = 1.0 / (R_in + R_wall + R_ao)

    alpha = k_rock / (C.RHO_ROCK * C.CP_ROCK)
    t = max(t_years, 1e-3) * 3.1536e7
    r_inf = r_w + 2.0 * np.sqrt(alpha * t)
    R_aw = 1.0 / (h_a * 2 * np.pi * r_w)
    R_rock = np.log(r_inf / r_w) / (2 * np.pi * k_rock)
    UAo = 1.0 / (R_aw + R_rock)
    return UAi, UAo, P


def solve(m_dot=2.0, T_inj=40.0, G=C.GEOTHERM_GRADIENT, T_surf=C.SURFACE_TEMP,
          target_rock_T=C.TARGET_ROCK_TEMP, Q_face=28000.0, t_years=1.0,
          k_ins=C.K_PIPE_INSULATION, k_rock=C.K_ROCK, n_nodes=160, verbose=True,
          geotherm=None, target_depth=None):
    # geotherm: optional callable z[m]->T_rock[degC] (e.g. a layered site profile).
    # If given, target_depth sets well length; else linear T_surf + G z to target_rock_T.
    if geotherm is not None:
        L = target_depth
        Trock = geotherm
    else:
        L = (target_rock_T - T_surf) / G
        Trock = lambda zz: T_surf + G * zz
    z = np.linspace(0, L, n_nodes)

    def odes(zz, y):
        Td, Tu = y
        UAi, UAo, P = conductances(zz, Td, Tu, m_dot, t_years, k_ins, k_rock)
        cpd = wt.CP(np.column_stack([np.clip(Td, 1, 480), np.clip(P, 1e5, 220e6)]))
        cpu = wt.CP(np.column_stack([np.clip(Tu, 1, 480), np.clip(P, 1e5, 220e6)]))
        dTd = UAi * (Tu - Td) / (m_dot * cpd)
        dTu = (-UAo * (Trock(zz) - Tu) + UAi * (Tu - Td)) / (m_dot * cpu)
        return np.vstack([dTd, dTu])

    def bc(ya, yb):
        cp_L = float(wt.cp(np.array([yb[0]]), np.array([P_of_z(L)]))[0])
        return np.array([ya[0] - T_inj,
                         yb[1] - yb[0] - Q_face / (m_dot * cp_L)])

    y0 = np.vstack([np.linspace(T_inj, T_inj + 30, n_nodes),
                    np.linspace(T_inj + 30, target_rock_T * 0.7, n_nodes)])
    sol = solve_bvp(odes, bc, z, y0, max_nodes=20000, tol=1e-4)

    zz = np.linspace(0, L, 400)
    Td, Tu = sol.sol(zz)
    # Heat delivered to the surface plant. The loop is PRESSURIZED, so the hot
    # return stays liquid (evaluating cp at 1 atm would wrongly treat >100 C
    # return as steam, cp~2000, halving the result). Use the loop operating
    # pressure P_OP and the mean temperature of the product stream.
    cp_prod = float(wt.cp(np.array([0.5 * (Tu[0] + Td[0])]), np.array([P_OP]))[0])
    Q_product = m_dot * cp_prod * (Tu[0] - Td[0])

    Trock_arr = np.array([Trock(zi) for zi in zz]) if geotherm is not None else T_surf + G * zz
    res = dict(sol=sol, z=zz, Td=Td, Tu=Tu, Trock=Trock_arr, L=L,
               T_bottom_delivered=Td[-1], T_return_surface=Tu[0],
               Q_product=Q_product, m_dot=m_dot, success=sol.success,
               P_bottom=P_of_z(L))
    if verbose:
        _print(res, T_inj, Q_face, t_years, k_ins)
    return res


def pump_power(m_dot, L, T_avg=80.0):
    P = P_of_z(L / 2)
    rho, cp, mu, k = (float(x[0]) for x in
                      wt.props(np.array([T_avg]), np.array([P])))
    dP_total = 0.0
    for area, Dh, length in [(A_pipe, Dh_pipe, L), (A_ann, Dh_ann, L)]:
        v = m_dot / (rho * area)
        Re = rho * v * Dh / mu
        f = 64 / Re if Re < 2300 else 0.316 * Re ** -0.25
        dP_total += f * (length / Dh) * 0.5 * rho * v ** 2
    return dP_total * (m_dot / rho), dP_total


def _print(r, T_inj, Q_face, t_years, k_ins):
    print("=" * 72)
    print("COUPLED 1-D COUNTERFLOW BOREHOLE  (Model 1, stage 2)")
    print("=" * 72)
    print(f"  Depth L                 : {r['L']/1000:.1f} km   (bottom rock {r['Trock'][-1]:.0f} C)")
    print(f"  Mass flow (per leg)     : {r['m_dot']:.2f} kg/s")
    print(f"  Injection temp          : {T_inj:.0f} C ; insulation k={k_ins} W/mK ; t={t_years} yr")
    print(f"  Face heat load Q_face   : {Q_face/1000:.1f} kW")
    print(f"  BVP converged           : {r['success']}  (P_bottom {r['P_bottom']/1e6:.0f} MPa)")
    print("-" * 72)
    pp, dP = pump_power(r['m_dot'], r['L'])
    ok = r['T_bottom_delivered'] < C.BHA_SURVIVAL_TEMP
    print(f"  >> Bit-delivered temp T_d(L) : {r['T_bottom_delivered']:6.1f} C   "
          f"(ceiling {C.BHA_SURVIVAL_TEMP:.0f} C)  {'OK' if ok else 'FAIL'}")
    print(f"  >> Return temp at surface    : {r['T_return_surface']:6.1f} C")
    print(f"  >> Thermal power exported    : {r['Q_product']/1e6:6.2f} MW_th")
    print(f"  >> Pump (friction) power     : {pp/1e3:6.1f} kW   (dP~{dP/1e6:.1f} MPa)")
    if pp > 0:
        print(f"  >> Thermal/pump ratio        : {r['Q_product']/max(pp,1):6.0f}x")
    print("=" * 72)


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    solve(m_dot=2.0, T_inj=40.0, Q_face=28000.0, t_years=1.0)
    print("\n### mass-flow sweep (insulated downcomer, t=1yr) ###")
    print(f"{'m_dot':>7}{'T_bit[C]':>10}{'T_surf[C]':>11}{'Q_MWth':>9}{'pump_kW':>9}{'conv':>6}")
    for md in (0.5, 1.0, 2.0, 5.0, 10.0, 20.0):
        r = solve(m_dot=md, t_years=1.0, verbose=False)
        pp, _ = pump_power(md, r['L'])
        print(f"{md:7.1f}{r['T_bottom_delivered']:10.1f}{r['T_return_surface']:11.1f}"
              f"{r['Q_product']/1e6:9.2f}{pp/1e3:9.1f}{str(r['success']):>6}")
    print("\n### downcomer insulation sensitivity (m_dot=2 kg/s) ###")
    print(f"{'k_ins':>8}{'T_bit[C]':>10}{'Q_MWth':>9}")
    for ki in (0.02, 0.10, 1.0, 45.0):
        r = solve(m_dot=2.0, k_ins=ki, verbose=False)
        print(f"{ki:8.2f}{r['T_bottom_delivered']:10.1f}{r['Q_product']/1e6:9.2f}")
    print("\n### operating-time sensitivity (rock cooldown, m_dot=2 kg/s) ###")
    print(f"{'years':>8}{'T_bit[C]':>10}{'Q_MWth':>9}")
    for ty in (0.1, 1.0, 5.0, 30.0):
        r = solve(m_dot=2.0, t_years=ty, verbose=False)
        print(f"{ty:8.1f}{r['T_bottom_delivered']:10.1f}{r['Q_product']/1e6:9.2f}")
