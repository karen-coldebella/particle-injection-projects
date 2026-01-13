# Particle Injection Projects

This repository contains analysis tools for identifying and studying particle injection events in the Earth's magnetosphere using THEMIS and Van Allen Probes data.

## Objectives
- Identify particle injections using magnetic field dipolarization signatures.
- Combine THEMIS FGM and SST observations.
- Combine Van Allen Probes HOPE and Mageis observations
- Support case studies and statistical analyses of plasma sheet injections.

## Data
The analysis uses Level-2 THEMIS data:
- FGM (Fluxgate Magnetometer)
- SST (Solid State Telescope)

And Van Allen Probes data:
- HOPE
- Mageis

Data files are not versioned in this repository.

## Project structure

```text
.
├── injection_themis.py
├── Data/
│   ├── FGM/
│   └── SST/
├── Figures/
└── README.md
