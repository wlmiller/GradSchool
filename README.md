Grad School Code
===

This repository contains example code from my time working toward my Ph.D. from Columbia University in Chemical Physics.

## [Janus_MD.c](Janus_MD.c)

This code runs a simulation of a generalization of so-called "Janus" particles.  Janus particles are particles which are one-half attractive and one-half repulsive; our generalization allows that ratio to vary from one-half.  For more information, see [1].

### Usage
```
./Janus_MD thetaMax e0 range
```
Where:

* `thetaMax` = The maximum angle of the attractive region (Janus particles are recovered for `thetaMax = pi`.
* `e0` = The dept of the attractive well (in unitless energy).
* `range` = The range of the attractive part of the interaction potential as compared to the repulsive part; `range = 1` is a standard Lennard-Jones potential.

## [aspherical_coexistence.c]
Performs a Monte Carlo of a box containing one-half fcc crystalline and one-half liquid-phase hard aspherical particles. Particles are constructed as described in [2], by placing some number of balls in an overlapping configuration.  The number of balls and degree of overlap are adjustable, allowing for control of the degree of asphericity.  See [3] for more information.

## References

1. W. L. Miller and A. Cacciuto. Hierarchical self-assembly of asymmetric amphiphatic spherical colloidal particles. Phys. Rev. E 80, 021404 (2009).
2. W. L. Miller, B. Bozorgui, and A. Cacciuto. Crystallization of hard aspherical particles. J. Chem. Phys. 132, 134901 (2010).
3. W. L. Miller and A. Cacciuto. On the phase behavior of hard aspherical particles. J. Chem. Phys. 133, 234903 (2010).
