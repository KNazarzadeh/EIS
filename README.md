# V2G
This repository is dedicated to Vehicle-to-Grid (V2G) technology, encompassing a technical documentation and corresponding simulation files.
V2G stands for Vehicle-to-Grid, a technology that allows electric vehicles (EVs) to not only receive power from the grid but also send electricity back to it.

## Schematic diagram of the V2G system
<img width="1322" height="897" alt="ChatGPT Image Nov 6, 2025, 08_48_25 AM" src="https://github.com/user-attachments/assets/9cb8252e-473a-4ea9-a001-aa8555216c0e" />

## File Structure
```bash
├───V2G
│  ├───datapreprocessing
│  ├───models
│  │   └───lithium_ion
│  │       ├───aging_models
│  │       ├───battery_models
│  │           ├───ESPM
│  │       └───submodels
│  │           ├───CCCV
│  │           └───FVM
│  ├───parameters
│  │   ├───ConcentrationParameters
│  │   ├───ConstantParameters
│  │   ├───DefaultParameters
│  │   ├───ElectricalParameters
│  │   ├───GeometricParameters
│  │   ├───KineticParameters
│  │   ├───ThermodynamicParameters
│  │   └───TransportParameters
│  ├───plots
│  ├───simulation
│  ├───Test
│  └───utils
```
# License
This project is licensed under the FZJ-Battery Modelling License
