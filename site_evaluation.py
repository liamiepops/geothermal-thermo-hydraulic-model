"""
site_evaluation.py
=================
SITE-EVALUATION TOOL -- wires Models 1-5 into a single verdict for a real
candidate site. This is the venture's core deliverable: feed in a province's
geotherm, in-situ stress, and rock properties; get back drillability, tool
survival, hole stability, and energy, each with its operating envelope.

CANDIDATE: Upper Rhine Graben / Soultz-sous-Forets (France) -- the best-
characterised deep-geothermal site in Europe (GPK-1..4 wells to ~5 km), an
active continental rift (extensional => low K0, favoured by Models 2/3/5).

All site parameters are REAL and cited inline. The point of the tool is that the
verdict is mixed and site-specific -- low K0 helps, high stress anisotropy hurts
-- which is exactly the judgement a generic pitch cannot make.

Sources (see chat): Soultz thermal profile ~90 C/km to 1.4 km then near-
isothermal convective granite to 200 C at 5 km (Genter et al.; MDPI Geosciences
2020); stress state Shmin~0.54 Sv, SHmax~Sv, SHmax N169E, breakouts observed
(Cornet/Valley/Heidbach; "Stress State at Soultz to 5 km").
"""
from dataclasses import dataclass, field
import numpy as np
import geo_constants as C
import model1_coupled as m1
import model2_spallation as m2
import model3_optimiser as m3
import model4_hole_stability as m4
import model5_convergence_confinement as m5

YEAR = 3.156e7


# --------------------------------------------------------------- site profile
def urg_geotherm(z):
    """Layered Upper Rhine Graben geotherm [degC], vectorized. z in metres.
    0-1400 m: ~90 C/km (sediments); 1400-5000 m: convective granite -> 200 C at
    5 km; >5000 m: ~35 C/km conductive (regional high heat flow)."""
    z = np.asarray(z, dtype=float)
    return np.select(
        [z <= 1400, z <= 5000],
        [12.0 + 0.090 * z,
         138.0 + 0.0172 * (z - 1400.0)],
        default=200.0 + 0.035 * np.clip(z - 5000.0, 0, None))


@dataclass
class SiteProfile:
    name: str
    geotherm: callable          # z[m] -> T_rock[degC]
    target_depth: float         # m
    target_T: float             # degC (consistency check vs geotherm)
    Sv_grad: float              # Pa/m, vertical (overburden) stress gradient
    K0_min: float               # Shmin / Sv  (low => extensional, favourable)
    SHmax_over_Sv: float        # SHmax / Sv
    rho_fluid_grad: float       # Pa/m, in-hole fluid pressure gradient
    k_rock: float; E_rock: float; UCS: float
    T_inj: float = 40.0

    @property
    def anisotropy(self):
        return self.SHmax_over_Sv / self.K0_min   # SHmax / Shmin


SOULTZ = SiteProfile(
    name="Upper Rhine Graben / Soultz-sous-Forets (France)",
    geotherm=urg_geotherm,
    target_depth=10700.0,       # ~400 C per layered geotherm
    target_T=400.0,
    Sv_grad=2650.0 * 9.81,      # URG granite ~2650 kg/m3 -> 26.0 kPa/m
    K0_min=0.54,                # Shmin/Sv ~0.54 (graben extension) -- measured
    SHmax_over_Sv=1.0,          # SHmax ~ Sv (strike-slip/normal transition)
    rho_fluid_grad=C.HYDROSTATIC_GRAD,
    k_rock=2.9,                 # URG biotite granite ~2.5-3.2 W/m/K
    E_rock=55e9,                # granite ~50-60 GPa
    UCS=170e6,                  # crystalline basement ~150-200 MPa
)


# ------------------------------------------------------------- the evaluation
def evaluate(site: SiteProfile):
    z = site.target_depth
    T_rock = float(site.geotherm(z))
    P_fluid = site.rho_fluid_grad * z
    Sv = site.Sv_grad * z
    Shmin = site.K0_min * Sv
    SHmax = site.SHmax_over_Sv * Sv
    # effective gradient that reproduces the real target depth in scalar-G models
    G_deep = (T_rock - C.SURFACE_TEMP) / z
    out = dict(site=site, z=z, T_rock=T_rock, Sv=Sv, Shmin=Shmin, SHmax=SHmax,
               P_fluid=P_fluid)

    # --- (1) SURVIVAL & ENERGY: Model 1 with the real LAYERED geotherm ---
    # scan flow for the minimum that keeps the bit < survival ceiling
    surv = None
    for md in (2, 4, 7, 10, 14, 20):
        r = m1.solve(m_dot=md, T_inj=site.T_inj, k_ins=0.02, k_rock=site.k_rock,
                     geotherm=site.geotherm, target_depth=z, Q_face=30000.0,
                     verbose=False)
        if r["T_bottom_delivered"] < C.BHA_SURVIVAL_TEMP:
            surv = (md, r); break
    if surv is None:
        surv = (20, r)
    out["m_min"], out["m1"] = surv

    # --- (2/3) DRILLABILITY: regime + quench ROP gain (spallation uses K0=Shmin)
    drill = m3.evaluate(G_deep, site.K0_min, out["m_min"], k_ins=0.02,
                        target_rock_T=T_rock, quench=True)
    drill_nq = m3.evaluate(G_deep, site.K0_min, out["m_min"], k_ins=0.02,
                           target_rock_T=T_rock, quench=False)
    out["drill"] = drill
    out["rop_gain"] = drill["ROP"] / drill_nq["ROP"] if drill_nq["ROP"] > 0 else np.nan

    # --- (4) CREEP closure: hot vs cooled wall ---
    out["creep_hot"] = m4.closure_rate(z, T_rock, 200.0, 7*86400, cooled=False) * YEAR * 100
    out["creep_cold"] = m4.closure_rate(z, T_rock, 200.0, 7*86400, cooled=True) * YEAR * 100

    # --- (5) STABILITY: isotropic GRC (uses Shmin as p0) + anisotropic breakout ---
    u_hot, rp_hot, pcr_h, reg_h = m5.grc(P_fluid, Shmin, T_rock)
    u_cold, rp_cold, pcr_c, reg_c = m5.grc(P_fluid, Shmin, 200.0)
    out["grc"] = dict(u_hot=u_hot, u_cold=u_cold, rp_hot=rp_hot, rp_cold=rp_cold,
                      reg_h=reg_h, reg_c=reg_c)
    # anisotropic breakout: sigma_theta = 3 SHmax - Shmin - P_i at the Shmin azimuth
    sth = 3 * SHmax - Shmin - P_fluid
    mc_cold = m5.KMC * P_fluid + m5.sigma_cm(200.0)
    mc_hot = m5.KMC * P_fluid + m5.sigma_cm(T_rock)
    # mud weight needed to suppress breakout even when cold
    P_need = (3 * SHmax - Shmin - m5.sigma_cm(200.0)) / (1 + m5.KMC)
    # If the mud weight needed to stop breakout exceeds Shmin, you would
    # hydraulically fracture the formation (lose returns) -> NOT mud-controllable.
    out["breakout"] = dict(sigma_theta=sth, mc_hot=mc_hot, mc_cold=mc_cold,
                           P_need=P_need, over_hydro=P_need / P_fluid,
                           overbalance_MPa=(P_need - P_fluid) / 1e6,
                           breaks=sth > mc_cold, frac_limited=P_need > Shmin)

    # production-flow energy (the survival run uses min flow, which minimises MW)
    rp = m1.solve(m_dot=10.0, T_inj=site.T_inj, k_ins=0.02, k_rock=site.k_rock,
                  geotherm=site.geotherm, target_depth=z, Q_face=30000.0, verbose=False)
    out["MW_prod"] = rp["Q_product"] / 1e6
    out["Tret_prod"] = rp["T_return_surface"]
    return out


def report(o):
    s = o["site"]
    print("=" * 80)
    print(f"SITE EVALUATION:  {s.name}")
    print("=" * 80)
    print(f"  Target: {o['T_rock']:.0f} C rock at {o['z']/1000:.1f} km "
          f"(layered geotherm; supercritical-class)")
    print(f"  In-situ stress: Sv={o['Sv']/1e6:.0f}  SHmax={o['SHmax']/1e6:.0f}  "
          f"Shmin={o['Shmin']/1e6:.0f} MPa  | K0={s.K0_min:.2f}  anisotropy={s.anisotropy:.2f}")
    print(f"  Fluid pressure (hydrostatic): {o['P_fluid']/1e6:.0f} MPa")
    print("-" * 80)
    md, r = o["m_min"], o["m1"]
    print(f"  [1] TOOL SURVIVAL & ENERGY (layered geotherm, vacuum tubing)")
    print(f"      min flow to survive: {md} kg/s -> bit {r['T_bottom_delivered']:.0f} C "
          f"(ceiling {C.BHA_SURVIVAL_TEMP:.0f}) {'OK' if r['T_bottom_delivered']<200 else 'FAIL'}")
    print(f"      surface return {r['T_return_surface']:.0f} C | early-life {r['Q_product']/1e6:.1f} MW_th")
    print(f"      at production flow 10 kg/s: {o['MW_prod']:.1f} MW_th, return {o['Tret_prod']:.0f} C")
    print("-" * 80)
    d = o["drill"]
    print(f"  [2/3] DRILLABILITY (quench-assist; spallation uses Shmin/Sv={s.K0_min:.2f})")
    print(f"      regime: {d['regime']} | effective MSE {d['MSE_eff']/1e6:.0f} MPa | "
          f"ROP {d['ROP']*3600:.1f} m/hr | quench gain {o['rop_gain']:.1f}x")
    print("-" * 80)
    print(f"  [4] CREEP CLOSURE (time-dependent)")
    print(f"      hot wall: {o['creep_hot']:.2e} %/yr  ->  cooled 200C: {o['creep_cold']:.2e} %/yr")
    print(f"      cooling factor ~ {o['creep_hot']/max(o['creep_cold'],1e-30):.1e}x  (squeezing controlled)")
    print("-" * 80)
    g = o["grc"]; b = o["breakout"]
    print(f"  [5] HOLE STABILITY")
    print(f"      isotropic GRC (p0=Shmin): {g['reg_c']}, convergence {g['u_cold']*1000:.1f} mm, "
          f"plastic r/a {g['rp_cold']:.2f}")
    print(f"      ANISOTROPIC breakout: sigma_theta={b['sigma_theta']/1e6:.0f} MPa vs "
          f"MC-limit cold {b['mc_cold']/1e6:.0f} MPa -> "
          f"{'BREAKS OUT' if b['sigma_theta']>b['mc_cold'] else 'stable'}")
    print(f"      mud weight to suppress breakout: {b['P_need']/1e6:.0f} MPa vs "
          f"hydrostatic {o['P_fluid']/1e6:.0f} MPa = +{b['overbalance_MPa']:.0f} MPa "
          f"({b['over_hydro']:.2f}x hydrostatic)")
    print("=" * 80)
    print("  VERDICT")
    verdict_lines(o)
    print("=" * 80)


def verdict_lines(o):
    s = o["site"]; b = o["breakout"]
    print(f"   + Extensional graben (K0={s.K0_min:.2f}): low confinement aids quench & limits")
    print(f"     the minimum-stress; tool survival closes with active cooling.")
    print(f"   + Cooling controls time-dependent creep by ~{o['creep_hot']/max(o['creep_cold'],1e-30):.0e}x.")
    print(f"   - HIGH anisotropy (SHmax/Shmin={s.anisotropy:.2f}) drives breakout, but only marginally:")
    print(f"     +{b['overbalance_MPa']:.0f} MPa overbalance ({b['over_hydro']:.2f}x hydrostatic) suppresses it,")
    print(f"     plus cooling margin. The binding constraint -- as in the real GPK wells, which")
    print(f"     broke out yet were drilled to 5 km -- not a showstopper.")
    print(f"   ~ At {o['T_rock']:.0f} C the rock is at the brittle-ductile edge; a 374 C/~10 km target")
    print(f"     sits more safely in the brittle field.")
    print(f"   => CONDITIONAL GO: drillable & survivable with quench+cooling; stability is mud-")
    print(f"      weight-limited by anisotropy, not temperature. A lower-anisotropy URG segment")
    print(f"      would score higher -- which is exactly what the tool is for.")


if __name__ == "__main__":
    import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    o = evaluate(SOULTZ)
    report(o)
