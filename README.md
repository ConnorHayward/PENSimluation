# LOSim

Working branch for GEANT4 simulaiton of 6 PMT set up for measuring the light yield of PEN.

Includes predefined geometry for existing 6 PMT set up at MPI, including all 6 PMTs, 3d printed holding strucutre, table and aluminium breadboard.
Also includes model of the 3d printed "sandwhich" used to encase the Cs137 source, which is also modeled.

The following data files are produced:

Light Output - number of photons detected at the chosen detector type.
Energy deposited in target, in MeV.
Light Yield - number of photons produced in an event.
Ratio of Light Output to Light Yield.
Number of photons leaving the target. Geometry independent version of light output.

To use: Navigate to build folder and enter:

cmake -DGeant4_DIR=/opt/geant4/lib64/Geant4-10.3.0 ..

This creates the makefile as needed. Then, run make to create the simulation. ./PEN runs the program.

Input values for PEN in the simulation:

| Property | Value | Reference |
|:--------:|:-----:|:---------:|
|Light Yield|10500 photon/MeV|Nakamura|
