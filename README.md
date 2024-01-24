# REPROBUS

REPROBUS (Reactive Processes Ruling the Ozone Budget in the Stratosphere) is a 3D chemistry-transport model calculating the temporal evolution of 55 chemical species via 147 chemical reactions (Lefèvre et al, 1998). Heterogeneous reactions are also taken into account. Liquid PSCs (Polar Stratospheric Clouds) are described using the model of Carslaw et al (1995). Wind and temperature fields from ECMWF analyses data are imposed every 3 hours. The time step of the chemical calculations is 15 minutes.

<p align="center">
    <img width="300" src="reprobus_logo.png" alt="REPROBUS logo">
</p>

## Requirements
The REPROBUS tool is containerized into a Singularity container so one must have Singularity installed on the host system intended for simulations.

## Installation
```
git clone https://github.com/aeris-data/reprobus.git
sudo singularity build ./reprobus.sif ./reprobus-ifort-container.def
```
The `singularity build` command will build the container `reprobus.sif` from its definition file, using the source files got from the git repo; so for the build it is important to call the command from the git repo directory that one has made. ⚠️ ***The build requires sudo rights.*** Afterwards, the sif image can be placed anywhere (even on another system) independently of the source files. To run the image no sudo rights are required.

## Usage
The main script is `reprobus-simulation.sh` which needs the input configuration file reprobus-user.conf (which can be renamed, the name is not important). This bash script handles user's input parameters, launch simulations and post-process simulation results. The main usage is 
```
./reprobus-simulation.sh --config reprobus-user.conf
```
The script must be launched inside the Singularity container. The outputs of the simulation are : simulation results for the 40 observation stations, the daily REPROBUS recovery files (also needed for the ulterior simulations) and the PNG files of the results visualisation . More details about input/output and folder structure are in the manual `xxx.pdf`.

## Input meteorological data extraction
The input data for the simulations is meteorological data : wind, temperature, humidity and logarithm of surface pressure, coming from the ECMWF database. To extract and prepare the data in the correct format, the script `reprobus-extract-ecmwf.sh` should be used:
```
./reprobus-extract-ecmwf.sh --config ./ecmwf.conf
```
This script handles the data extraction and its formatting for the REPROBUS tool. The user can configure the start and end date of the data, as well as the spatial resolution. The script must be launched on the ECMWF MARS server (ecs, hpc or other). The data extraction was tested with a member-state user account. Other more public accounts might customize the script based on the MARS services or APIs available for their type of user. The data is extracted and stored in the directory requested in the input configuration; afterwards, the data can be used for the simulation.