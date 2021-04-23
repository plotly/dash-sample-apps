# Dash: Vehicle Geometry with OpenFOAM

by [Nathan Rooy](https://github.com/nathanrooy)

## Background

This project initially began with vehicle geometry that was originally created for graphic design and rendering purposes. Therefore, an initial cleaning process was required before any CFD could take place. This cleaning phase included sealing body panels, improving triangle quality, and creating geometry that wasn’t included in the initial model (intake plenum, cooling ductwork, contact patches). Once the model was watertight and suitable for simulation, snappyHexMesh was employed for creating the mesh. Given that the final output needs to be run in a browser, the mesh resolution was purposely kept coarse in an attempt to keep the memory/compute requirements low. 

Next, the boundary conditions were set and the equations of fluid dynamics were then discretized across the mesh volume and solved using OpenFOAM. The aerodynamics of open-wheeled vehicles such as this, tend to produce  a substantial amount of flow separation which requires resolving turbulence down to both a very fine spatial and temporal domain. However, given the aforementioned browser constraints, a simpler approximation was needed. Because of this, the RANS equations using the SIMPLE approach with the SST k-ω turbulence model was used.

Results were then post-processed in ParaView (paraFoam) and exported. From here, a simple Python app was built using Plotly Dash while leveraging VTK. The end result is a web app in which surface velocity, surface pressure, and a Cp isosurface can be visualized. Vehicle geometry visibility can be toggled as well.

## Installation and Usage
Install the dependencies found in `requirements.txt`:
```
pip install -r requirements.txt
```

Next, run `app.py` which will launch a local Dash server:
```
python app.py
```

## References
- Original repository: https://github.com/nathanrooy/dash-vtk-cfd
- SIMPLE algorithm: https://en.wikipedia.org/wiki/SIMPLE_algorithm
- RANS: https://en.wikipedia.org/wiki/Reynolds-averaged_Navier%E2%80%93Stokes_equations
- SST k-ω turbulence model: https://www.cfd-online.com/Wiki/SST_k-omega_model
- OpenFOAM: https://openfoam.org/
