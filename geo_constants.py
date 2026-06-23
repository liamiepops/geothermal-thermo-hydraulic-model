"""
geo_constants.py
================
Real, sourced constants for the deep/superhot-rock geothermal drilling models.

EVERY number here has a comment with (a) the value, (b) the range/uncertainty,
(c) a source you can cite and defend. Where a value varies strongly with
temperature/pressure it is flagged TEMP-DEP and handled by a function elsewhere
(e.g. water properties come from IAPWS-95, not from a frozen constant).

Design philosophy: these are *defensible baseline* numbers, not best cases.
When a number is uncertain, the model exposes it as a knob so you can run the
sensitivity yourself and own the conclusion.
"""

# ----------------------------------------------------------------------------
# SITE / GEOTHERM
# ----------------------------------------------------------------------------
# Geothermal gradient G [K/m]. Continental crust:
#   - Stable cratons:            ~15-25 K/km   (avoid; deep & cold)
#   - Average continental:       ~25-30 K/km   (Fridleifsson/IEA baseline)
#   - Elevated / extensional
#     (Basin & Range, where the
#      DOE Utah FORGE site sits): ~35-50 K/km sustained to mid-crust
#   - High-enthalpy volcanic
#     near-surface (Iceland):     60-100 K/km but DECLINES with depth
# We take an elevated-but-generic continental value as the defensible "drill
# almost anywhere with a decent gradient" design point, and sweep it.
# Utah FORGE measured ~200 C at 2.5 km (~75 K/km near surface, basin-enhanced);
# that near-surface number does NOT extrapolate linearly to 12 km, so we do not
# use it as G. Source: Utah FORGE resource assessment; Cascade Institute
# "Drilling for Superhot Geothermal Energy" (2024).
GEOTHERM_GRADIENT = 0.035          # K/m  (35 K/km) baseline design point
GEOTHERM_GRADIENT_RANGE = (0.025, 0.050)  # K/m sweep range
SURFACE_TEMP = 15.0                # degC  mean annual surface temperature

# Target bottom-hole rock temperature [degC]. "Superhot" is conventionally
# >374 C (water critical T). We design to 450 C as a representative superhot
# target (Holtzman et al. 2023 work the 400-600 C band; Quaise targets ~500 C).
TARGET_ROCK_TEMP = 450.0           # degC

# Hydrostatic pressure gradient [Pa/m]. A water column averages ~9.5-10 MPa/km
# depending on density profile (cold dense near surface, less dense at depth).
# We use 10 MPa/km as a round, slightly conservative (high) value; the model
# can also integrate rho(z)*g for a self-consistent column.
HYDROSTATIC_GRAD = 10.0e6 / 1000.0  # Pa/m  (10 MPa/km)
# Lithostatic (overburden) gradient, for the confining-stress reality check in
# Model 2: rho_rock * g ~ 2700*9.81 = 26.5 kPa/m ~ 26.5 MPa/km.
LITHOSTATIC_GRAD = 2700.0 * 9.81    # Pa/m  (~26.5 MPa/km)

g = 9.81                            # m/s^2

# ----------------------------------------------------------------------------
# ROCK (granite / crystalline basement baseline)
# ----------------------------------------------------------------------------
# Thermal conductivity k_rock [W/m/K]. Granite ~2.5-3.5 at room T; DECREASES
# with temperature, ~2.0-2.5 at 400-500 C. Use 2.5 (matches your envelope calc).
K_ROCK = 2.5                       # W/m/K  TEMP-DEP (declines at high T)
RHO_ROCK = 2700.0                  # kg/m^3   granite
CP_ROCK = 1000.0                   # J/kg/K   granite specific heat (~790-1100; rises with T)
# Thermal diffusivity alpha_rock = k/(rho*cp) ~ 0.93e-6 m^2/s
ALPHA_ROCK = K_ROCK / (RHO_ROCK * CP_ROCK)

# Thermo-mechanical (for Model 2 thermal-stress / spallation):
E_ROCK = 50.0e9                    # Pa   Young's modulus, granite 40-70 GPa; falls at high T & with cracking
NU_ROCK = 0.25                     # -    Poisson ratio 0.2-0.3
ALPHA_THERMAL = 8.0e-6             # 1/K  linear thermal expansion; granite 7-10e-6, RISES with T
TENSILE_STRENGTH = 10.0e6          # Pa   granite tensile strength 7-15 MPa; FALLS with T and pre-cracking
UCS = 200.0e6                      # Pa   unconfined compressive strength 150-250 MPa (for context only)

# ----------------------------------------------------------------------------
# DRILLING / BOTTOM-HOLE ASSEMBLY (BHA)
# ----------------------------------------------------------------------------
BIT_DIAMETER = 0.216               # m    8.5" bit (common deep-hole size)
CAVITY_RADIUS = 0.108              # m    ~ bit radius, hemispherical bottom-hole cavity
# Mechanical Specific Energy MSE [J/m^3 == Pa]: energy to remove unit volume.
# Hard crystalline rock: 100-400 MPa typical, up to ~1 GPa with dull bits.
# (Teale 1965 definition; Pessier & Fear; common drilling-optimization range.)
MSE = 300.0e6                      # Pa  (300 MJ/m^3) baseline hard-rock
MSE_RANGE = (100.0e6, 1000.0e6)
# Rate of penetration ROP [m/s]. Hard granite conventional: 1-5 m/hr.
# Quench-spallation demo (Holtzman 2023): ~30 mm/min = 1.8 m/hr at 480 C,
#   ambient pressure (optimistic, unconfined).
ROP = 2.0 / 3600.0                 # m/s  (2 m/hr) baseline
ROP_RANGE = (1.0/3600.0, 30.0e-3/60.0)  # 1 m/hr .. 30 mm/min

# Survival temperature of the bottom-hole assembly [degC].
# - Conventional MWD/LWD electronics: rated ~150-175 C.
# - High-temperature commercial tools: ~200 C.
# - Special/experimental electronics: ~300 C (very limited).
# - Above ~300 C effectively NO commercial downhole electronics survive.
# This is THE hard bottleneck the "always-cooling" thesis must beat.
BHA_SURVIVAL_TEMP = 200.0          # degC  design ceiling (high-T tool class)
BHA_SURVIVAL_TEMP_HARD = 300.0     # degC  absolute ceiling for any electronics

# ----------------------------------------------------------------------------
# WELLBORE GEOMETRY (coaxial / concentric, Eavor-like)
# ----------------------------------------------------------------------------
# Cold fluid DOWN an insulated central pipe -> reaches the bit cold;
# hot fluid UP the annulus -> picks up heat from the rock wall (the product).
# This routing is what lets one loop both cool the bit and harvest heat.
R_INNER_PIPE_IN = 0.050            # m  inner radius of central (downcomer) pipe
R_INNER_PIPE_OUT = 0.060           # m  outer radius of central pipe (wall+insulation)
R_WELL = 0.108                     # m  borehole / casing inner radius (annulus outer)
# Insulation of the central pipe is decisive: a poorly insulated centre pipe
# lets the cold downflow and hot upflow short-circuit thermally.
K_PIPE_INSULATION = 0.10           # W/m/K  vacuum-insulated tubing ~0.02-0.1; bare steel ~45 (bad)

# ----------------------------------------------------------------------------
# WATER / COOLANT  -> see water_props.py (IAPWS-95, real EOS). No frozen cp here
# on purpose: at 120 MPa, cp runs 3.9-4.9 kJ/kg/K over 15-450 C.
# ----------------------------------------------------------------------------
CP_WATER_NOMINAL = 4200.0          # J/kg/K  ONLY for quick hand-checks; models use IAPWS

if __name__ == "__main__":
    print(f"alpha_rock = {ALPHA_ROCK:.3e} m^2/s")
    d = (TARGET_ROCK_TEMP - SURFACE_TEMP) / (GEOTHERM_GRADIENT * 1000)
    print(f"Depth to reach {TARGET_ROCK_TEMP} C at {GEOTHERM_GRADIENT*1000:.0f} K/km: {d:.1f} km")
    for G in (0.025, 0.035, 0.050):
        d = (TARGET_ROCK_TEMP - SURFACE_TEMP) / (G * 1000)
        print(f"   G={G*1000:.0f} K/km -> {d:5.1f} km ; P_hydro~{HYDROSTATIC_GRAD*d*1000/1e6:.0f} MPa")
