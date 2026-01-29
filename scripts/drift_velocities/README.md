# Drift Velocities of Charged Particles in the Earth's Magnetosphere

## Overview
This code calculates the drift velocities of charged particles in the Earth's magnetic field.
The calculations are based on classical guiding center theory and are intended for applications
in radiation belt and particle injection studies.

## Important Considerations
The following assumptions and approximations are adopted:

- The Earth's magnetic field is approximated as a dipole.
- Particles are treated under the non-relativistic approximation.
- Only first-order gradient and curvature drift velocities are considered.

## Particle Populations
The particle populations of interest are protons and electrons with low to medium energies.
This energy range is characteristic of injected particle populations associated with
magnetospheric dynamics.

## References
- Roederer, J. G. (1970). *Dynamics of Geomagnetically Trapped Radiation*.  
  Physics and Chemistry in Space, Vol. 4. Springer, Berlin, Heidelberg.  
  https://doi.org/10.1007/978-3-642-65422-6_1

- Koskinen, H. E. J., & Kilpua, E. K. J. (2021). *Physics of Earth's Radiation Belts*.


## 1. Gradient and Curvature Drift Velocity

Under the non-relativistic approximation, the combined gradient and curvature
drift velocity \( \mathbf{V}_{CG} \) of a charged particle in a curved magnetic
field is given by:

$$
V_{CG} = \frac{T}{q\,B\,R_c}\left(1 + \cos^2 \alpha \right)
$$

where:
- \( T \) is the particle kinetic energy (J),
- \( q \) is the particle electric charge (C),
- \( B \) is the magnetic field strength (T),
- \( R_c \) is the radius of curvature of the magnetic field line (m),
- \( \alpha \) is the particle pitch angle.

### Special Case: \( \alpha = 90^\circ \)

For particles with a pitch angle of \( 90^\circ \), the expression simplifies to:

$$
V_{CG}^{90^\circ} = \frac{T}{q\,B\,R_c}
$$

This configuration is commonly used as a reference case in radiation belt studies
and corresponds to particles with purely perpendicular motion.


## 2. Magnetic Field Magnitude in a Dipolar Approximation

To compute the drift velocity, the magnetic field magnitude \( B \) and the
radius of curvature \( R_c \) must be expressed in terms of the parameters
of a dipolar magnetic field.

In a dipole approximation, the magnetic field strength at a distance \( r \)
from the center of the Earth and magnetic latitude \( \lambda \) is given by:

$$
B(r, \lambda) = \frac{k_0}{r^3}\sqrt{1 + 3\sin^2 \lambda}
$$

where:
- \( k_0 \) is the magnetic field strength at the equator on the Earth's surface,
- \( r \) is the distance from the center of the Earth (m),
- \( \lambda \) is the magnetic latitude (rad).

### Equatorial Magnetic Field

At the magnetic equator (\( \lambda = 0 \)), the field magnitude simplifies to:

$$
B_{\text{eq}} = \frac{k_0}{r^3}
$$

### Relation with the Magnetic Dipole Moment

The parameter \( k_0 \) can be expressed in terms of the Earth's magnetic dipole
moment \( M \) as:

$$
k_0 = \frac{\mu_0 M}{4\pi}
$$

where:
- \( \mu_0 \) is the permeability of free space  
  (\( \mu_0 = 4\pi \times 10^{-7} \, \text{T·m/A} \)),
- \( M \) is the Earth's magnetic dipole moment  
  (\( M \approx 7.96 \times 10^{22} \, \text{A·m}^2 \)).

Using these values, one obtains:

$$
k_0 \approx 8 \times 10^{15} \, \text{T·m}^3
$$

### Expression in Terms of the McIlwain \( L \)-Parameter

The radial distance \( r \) can be written as:

$$
r = L\,R_E
$$

where \( R_E = 6.371 \times 10^6 \, \text{m} \) is the Earth's radius.

Substituting into the equatorial magnetic field expression yields:

$$
B_{\text{eq}}(L) = \frac{8 \times 10^{15}}{L^3 R_E^3}
$$




## 3. Radius of Curvature of Dipolar Magnetic Field Lines

In a dipolar magnetic field, the radius of curvature \( R_c \) of a magnetic
field line at a given magnetic latitude \( \lambda \) and radial distance \( r \)
is given by:

$$
R_c(r,\lambda) =
\frac{r}{3}
\frac{\cos \lambda \left( 1 + 3\sin^2 \lambda \right)^{3/2}}
{2 - \cos^2 \lambda}
$$

where:
- \( r \) is the radial distance from the center of the Earth (m),
- \( \lambda \) is the magnetic latitude (rad).

### Equatorial Approximation

At the magnetic equator (\( \lambda = 0 \)), the expression simplifies to:

$$
R_{c,\text{eq}} = \frac{r}{3}
$$

### Expression in Terms of the McIlwain \( L \)-Parameter

Writing the radial distance as:

$$
r = L\,R_E
$$

the equatorial radius of curvature becomes:

$$
R_{c,\text{eq}}(L) = \frac{L\,R_E}{3}
$$


## 4. Final Expression for the Drift Velocity

For particles with a pitch angle of \( 90^\circ \), the combined gradient and
curvature drift velocity is given by:

$$
V_{CG}^{90^\circ} = \frac{T}{q\,B\,R_c}
$$

where \( T \) is the kinetic energy and \( q \) is the particle charge.

### Substitution of Dipole Field Quantities

At the magnetic equator, the magnetic field strength can be written as:

$$
B_{\text{eq}}(L) = \frac{B_0}{L^3}
$$

where \( B_0 = 3.12 \times 10^{-5} \, \text{T} \) is the equatorial magnetic field
strength at the Earth's surface.

The equatorial radius of curvature is:

$$
R_{c,\text{eq}}(L) = \frac{L\,R_E}{3}
$$

where \( R_E = 6.371 \times 10^6 \, \text{m} \) is the Earth's radius.

### Final Drift Velocity Expression

Substituting the expressions for \( B_{\text{eq}}(L) \) and \( R_{c,\text{eq}}(L) \)
into the drift velocity equation yields:

$$
V_{CG}(L) =
\frac{T\,L^2}{q}
\left(
\frac{3}{B_0 R_E}
\right)
$$

Evaluating the constant term using Earth parameters results in:

$$
V_{CG}(L) \approx \frac{T\,L^2}{q \times 63.78}
$$

This expression highlights the quadratic dependence of the drift velocity on the
McIlwain \( L \)-parameter and is valid under the assumptions described in the
previous sections.
