## Overview

This repository contains a work-in-progress 1D electrochemical model of a Li-O_2 battery cathode based primarily on:

K. Jiang et al., *Parameter sensitivity analysis and cathode structure optimization of a non-aqueous Li-O_2 battery model*, Journal of Power Sources 451 (2020) 227821.
https://doi.org/10.1016/j.jpowsour.2020.227821

Note: This is an independent reimplementation, and I am not affiliated with the authors.

## Status

This model is under active development.

- Core physics implemented: Butler-Volmer kinetics, O_2 transport, Li_2O_2 deposition, porosity evolution
- Model runs under no initial Li_2O_2, but fails otherwise

## Known issues:

- Differential-algebraic system is overdetermined
- Initialization fails when Li_2O_2 precipitation is seeded

## Goals:

- Achieve a well-posed cathode + separator formulation
- Reproduce Jiang et al. discharge curves
- Use the model as a learning platform for Li-air battery electrochemistry

## Environment:

Developed and tested with:
- Python 3.10.16 (packaged by Anaconda)
- PyBaMM 25.4.0
