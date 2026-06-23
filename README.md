# Deep / superhot-rock geothermal: a first-principles model stack + site-evaluation tool

A self-contained set of coupled physical models, in Python, for the engineering
questions that decide whether deep and superhot-rock geothermal wells are
drillable, survivable, stable, and worth the energy — and a tool that scores a
real candidate **site** against all of them.

Built from first principles with real material data (IAPWS-95 water properties,
cited rock and in-situ-stress parameters). Every constant is sourced and exposed
as a knob; every model states its own assumptions and limits. The goal was not a
polished product but a **defensible, pressure-tested analysis** — including
results that constrain or contradict the starting hypothesis.

> **Status / honesty note.** These are order-of-magnitude / directional models:
> 1-D or axisymmetric, mostly quasi-steady, with several coupling parameters that
> are uncertain and flagged as such. They are meant to *locate* the binding
> constraints and quantify trade-offs, not to replace a full reservoir/geomech
> simulator. The limitations are listed per model and summarised at the bottom.

---

## The questions, and what the models found

| # | Model | Question | Headline result |
|---|-------|----------|-----------------|
| 1 | `model1_coupled.py` | Can a circulating coolant loop keep the bottom-hole assembly survivable while exporting heat at 10–20 km? | Yes, in a flow/insulation window. Single closed loop is **conduction-limited (~3–4 MW)**; water stays dense/supercritical (no flashing) above ~2 km. |
| 2 | `model2_spallation.py` | Can a cold-coolant **quench** fracture hot rock at the cutting face? | Tensile spallation is **confinement-limited**: feasible only in low-K0 (extensional) crust below the brittle-ductile transition. Elsewhere it acts as an **MSE-reducer**, not a stand-alone driller. |
| 3 | `model3_optimiser.py` | Couple drilling and thermal balance: where does the integrated loop win? | Quench buys **1.5–2.3× ROP** (assisted regime), up to **7.5×** in a narrow spallation island. Produces a siting/operating map. |
| 4 | `model4_hole_stability.py` | Does cooling the hot ductile rock keep the hole open? | Cooling suppresses time-dependent **creep closure by 10⁵–10⁶×**. (Convective heat-resupply to the cooled zone is negligible — Péclet ≪ 1.) |
| 5 | `model5_convergence_confinement.py` | Is "breakout" at depth a collapse or a manageable yielded zone? | **Breakout ≠ collapse**: convergence is mm-scale and fluid-pressure-controlled. Cooling adds **+3–6 km** of stable depth. **Stress anisotropy** (SHmax/Shmin), not temperature, is the real limiter. |

### Integrating tool
- `site_evaluation.py` — runs Models 1–5 against a real `SiteProfile` (geotherm,
  in-situ stress, rock properties) and returns a consolidated verdict + envelope.
- `comparative_sites.py` — scores four European provinces side by side.

### Comparative result (European superhot candidates, target 400 °C)

| Site | Depth to 400 °C | SHmax/Shmin | Verdict | Why |
|------|-----------------|-------------|---------|-----|
| Pannonian Basin (HU) | 9.7 km | 1.50 | **GO** | high gradient → shallower; moderate stress |
| Larderello (IT) | 2.9 km | 1.36 | **GO** | extreme gradient → very shallow; lowest power |
| Upper Rhine Graben (FR) | 10.7 km | 1.85 | **CONDITIONAL** | marginal breakout (+3 MPa mud) |
| United Downs (UK) | 12.5 km | 2.55 | **NO-GO** | strike-slip: breakout exceeds frac gradient |

The governing axes drop out of the physics: **geothermal gradient sets
depth-to-superhot; stress anisotropy sets stability.** Figure:
`comparative_sites.png`.

---

## Key engineering conclusions

- **Closed-loop power is conduction-limited and diameter-insensitive.** 10×
  borehole diameter → only ~1.9× heat (logarithmic), while cutting cost scales as
  r². The levers are contact *length* (laterals), well count, or fracture area —
  not diameter. (`model1_diameter_scaling.py`)
- **Depth is a grade lever, not a high-power lever**, and is capped by the
  brittle-ductile transition (~400 °C → 5.5–19 km depending on gradient).
  Active cooling removes *tool survival* as the limit; the rock does the limiting.
  (`model1_depth_limits.py`)
- **The novelty is the integrated computational capability, not the mechanism.**
  Cold-quench rock fracture is prior art (Holtzman et al., *GRC Transactions* 47,
  2023; cryogenic/thermal-shock EGS stimulation). The contribution here is a
  coupled, real-property, site-scoring model stack that reproduces known
  trade-offs and flags site-specific constraints.

---

## Running it

```bash
pip install -r requirements.txt
python water_table.py        # one-time: build the cached IAPWS-95 property table (~2.5 min)
python site_evaluation.py    # full verdict for the Upper Rhine Graben / Soultz site
python comparative_sites.py  # four-site comparison table
```

Each model file is runnable standalone and prints its own analysis
(`python model1_coupled.py`, etc.). Figure scripts (`*_figures.py`) regenerate
the PNGs.

### Layout
```
geo_constants.py          sourced constants (geotherm, rock, drilling, geometry)
water_props.py            IAPWS-95 wrapper (point queries)
water_table.py            cached, vectorised IAPWS-95 property table (T,P grid)
model1_coupled.py         coupled 1-D counterflow borehole heat exchanger (BVP)
model1_steadystate.py     single-depth bottom-hole heat balance
model1_diameter_scaling.py / model1_depth_limits.py   scaling analyses
model2_spallation.py      thermoelastic quench fracture vs confining stress
model3_optimiser.py       coupled drilling<->thermal optimiser + siting map
model4_hole_stability.py  creep closure & convection (Péclet) check
model5_convergence_confinement.py   ground-reaction-curve hole stability
site_evaluation.py        SiteProfile + integrated verdict
comparative_sites.py      four European provinces
*_figures.py              figure generators
```

---

## Limitations (read before trusting any number)

- **Dimensionality:** Model 1 is 1-D; Models 2/5 are axisymmetric/plane-strain.
  Real stress anisotropy and 3-D geometry are partially captured (breakout uses
  SHmax/Shmin) but not fully.
- **Time dependence:** the rock-coupling and stability models are mostly
  quasi-steady. Multi-year reservoir thermal decline (production economics) is
  *not* modelled — deliberately, since the drilling face sees fresh rock.
- **Uncertain coupling parameters:** the quench→effective-MSE coupling
  (`χ_macro`, `χ_micro`) and the high-temperature rock flow law are
  literature-typical and **require experimental calibration**; results are given
  as sensitivities, not predictions.
- **Stress data:** deep in-situ stress (especially anisotropy below measured
  depths) is the dominant input uncertainty in the site verdicts.
- **Intact-rock assumption:** faulted/altered rock will behave worse than the
  intact-strength values used.

---

## Author

Applied mathematician (Research M.Sc.) and software engineer (~15 yrs,
TypeScript/Node, test automation → full-stack), working at the intersection of
applied maths, computational physics, and dependable software. This repository is
a self-directed first-principles study; contact details in the accompanying note.
