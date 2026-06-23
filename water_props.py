"""
water_props.py
==============
Thin wrapper over IAPWS-95 (the full Helmholtz-energy EOS, valid to high
pressure) giving the four transport/thermo properties the borehole models need
as functions of (T [degC], P [Pa]):

    rho  [kg/m^3]   density
    cp   [J/kg/K]   isobaric specific heat
    mu   [Pa.s]     dynamic viscosity
    k    [W/m/K]    thermal conductivity

Why IAPWS-95 and not IF97: deep-well pressures reach 100-200 MPa. IAPWS-IF97 is
only validated to 100 MPa; IAPWS-95 covers the full range and the supercritical
region continuously, so there is no spurious "phase boundary" in the model.

A simple in-memory cache keeps the BVP solver from re-solving the EOS millions
of times.
"""
from functools import lru_cache
from iapws import IAPWS95

KELVIN = 273.15


@lru_cache(maxsize=200_000)
def _props_cached(T_C_round, P_Pa_round):
    # Clamp into the IAPWS-95 validity / sane borehole domain.
    T_C_round = min(max(T_C_round, 1.0), 1000.0)
    P_Pa_round = min(max(P_Pa_round, 1.0e5), 1000.0e6)  # >=1 atm, <=1 GPa
    T = T_C_round + KELVIN
    P_MPa = P_Pa_round / 1.0e6
    w = IAPWS95(T=T, P=P_MPa)
    # IAPWS95 occasionally returns None on a property near a boundary; nudge P.
    if w.cp is None or w.mu is None or w.k is None or w.rho is None:
        w = IAPWS95(T=T, P=P_MPa * 1.001 + 1e-3)
    # IAPWS95 returns cp in kJ/kg/K, mu in Pa.s, k in W/m/K, rho in kg/m^3
    return (w.rho, w.cp * 1000.0, w.mu, w.k, str(w.phase))


def water(T_C, P_Pa):
    """Return dict of properties. T in degC, P in Pa. Rounded for caching."""
    # round to 0.05 C and 0.05 MPa: smooth enough for these models
    Tr = round(min(max(T_C, 1.0), 1000.0) * 20) / 20.0
    Pr = round(min(max(P_Pa, 1.0e5), 1000.0e6) / 5.0e4) * 5.0e4
    rho, cp, mu, k, phase = _props_cached(Tr, Pr)
    return dict(rho=rho, cp=cp, mu=mu, k=k, phase=phase)


def cp(T_C, P_Pa):
    return water(T_C, P_Pa)["cp"]


def rho(T_C, P_Pa):
    return water(T_C, P_Pa)["rho"]


if __name__ == "__main__":
    print(f"{'T[C]':>6}{'P[MPa]':>8}{'rho':>9}{'cp[J/kgK]':>11}{'mu[uPa.s]':>11}{'k':>8}  phase")
    for T in (15, 100, 200, 300, 374, 450):
        for P in (5e6, 22.1e6, 60e6, 120e6):
            w = water(T, P)
            print(f"{T:6.0f}{P/1e6:8.1f}{w['rho']:9.1f}{w['cp']:11.0f}"
                  f"{w['mu']*1e6:11.1f}{w['k']:8.3f}  {w['phase']}")
