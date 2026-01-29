"""
Drift Velocity Calculator for Charged Particles in Earth's Magnetosphere

This code calculates the gradient and curvature drift velocities and the
azimuthal drift period of protons and electrons in the Earth's dipolar
magnetic field.

Inputs:
- Particle energy (keV or MeV)
- McIlwain L parameter
- Particle type (proton or electron)

Outputs:
- Drift velocity [m/s]
- Drift period [minutes]

For detailed physical background, assumptions, and equations, please
refer to the README.md file in this repository.
"""


#%%
import numpy as np

# -----------------------------
# Physical constants
# -----------------------------
q_e = 1.602176634e-19      # C, elementary charge
R_E = 6.371e6              # m, Earth radius
B0 = 3.12e-5               # T, equatorial magnetic field at Earth's surface
C_EARTH = 63.78            # constant derived from Earth parameters


def drift_velocity(E, L, particle="proton", energy_unit="keV"):
    """
    Computes the gradient and curvature drift velocity for a charged particle
    with 90° pitch angle at the magnetic equator.

    Parameters
    ----------
    E : float
        Particle energy.
    L : float
        McIlwain L parameter.
    particle : str
        'proton' or 'electron'.
    energy_unit : str
        'keV' or 'MeV'.

    Returns
    -------
    V_CG : float
        Drift velocity [m/s]
    """

    # Charge sign
    if particle.lower() == "proton":
        q = q_e
    elif particle.lower() == "electron":
        q = -q_e
    else:
        raise ValueError("particle must be 'proton' or 'electron'")

    # Energy conversion to Joules
    if energy_unit.lower() == "kev":
        T = E * 1.602176634e-16  # J
    elif energy_unit.lower() == "mev":
        T = E * 1.602176634e-13  # J
    else:
        raise ValueError("energy_unit must be 'keV' or 'MeV'")

    # Drift velocity (Eq. 8)
    V_CG = (T * L**2) / (q * C_EARTH)

    return V_CG


def drift_period(E, L, particle="proton", energy_unit="keV"):
    """
    Computes the azimuthal drift period around the Earth.

    Returns
    -------
    T_D : float
        Drift period [s]
    """

    V_CG = drift_velocity(E, L, particle, energy_unit)

    # Drift path length (circumference)
    r = L * R_E
    d = 2 * np.pi * r

    # Drift period
    T_D = d / np.abs(V_CG)

    return T_D


# -----------------------------
# Example
# -----------------------------
if __name__ == "__main__":

    E = 100       # keV
    L = 5.0
    particle = "electron"

    V = drift_velocity(E, L, particle)
    T_D = drift_period(E, L, particle)

    print(f"Particle: {particle}")
    print(f"Energy: {E} keV")
    print(f"L-shell: {L}")
    print(f"Drift velocity: {V:.2e} m/s")
    print(f"Drift period: {T_D/3600:.2f} hours")

# %%
