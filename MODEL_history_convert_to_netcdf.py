from scipy.io import FortranFile
import numpy as np
import struct
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import datetime as dt
import os
import netCDF4 as nc
import glob
import logging
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import math
import sys
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import locale

SPECIES_1 = ['N2O','CH4','H2O','NOy','HNO3','N2O5','Cly','Ox','CO','OClO','Passive Ox','H2SO4','HCl','ClONO2','HOCl','Cl2','H2O2','ClNO2','HBr','BrONO2','NOx','HNO4','ClOx','BrOx','Cl2O2','HOBr','BrCl','CH2O','CH3O2','CH3O2H','CFC-11','CFC-12','CFC-113','CCl4','CH3CCl3*','CH3Cl','HCFC-22*','CH3Br','H-1211*','H-1301','Bry','CH2Br2*','HNO3 GAS']
SPECIES_2 = ['O(1D)','OH','Cl','O(3P)','O3','HO2','NO2','NO','Br','N','ClO','BrO','NO3','H','CH3']

PLOT_LEVELS = {
"N2O" : [0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,220.,240.,260.,280.,300.,320.,340.,360.],
"HCl" : [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0],
"ClONO2" : [0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2],
"NOx" : [0,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.,4.25,4.5,4.75],
"ClOx" : [0.,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0],
"BrOx" : [0,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.],
"HNO3" : [0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0],
"Surface_Area" : [0,1.,2.,3.,4.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.],
"O3loss" : [-100.,-90.,-80.,-70.,-60.,-50.,-40,-35.,-30.,-25.,-20.,-15.,-10.,-5.,0.,5.],
"O3" : [0.,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.,6.4],
"NO2" : [0,0.25,0.5,0.75,1.,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.,4.25,4.5,4.75]}

PLOT_COEFFS = {
"N2O" : 1e-9,
"HCl" : 1e-9,
"ClONO2" : 1e-9,
"NOx" : 1e-9,
"ClOx" : 1e-9,
"BrOx" : 1e-12,
"HNO3" : 1e-9,
"Surface_Area" : 1e-8,
"O3loss" : 1e-2,
"O3" : 1e-6,
"NO2" : 1e-9}

PLOT_UNITS = {
"N2O" : "(ppbv)",
"HCl" : "(ppbv)",
"ClONO2" : "(ppbv)",
"NOx" : "(ppbv)",
"ClOx" : "(ppbv)",
"BrOx" : "(pptv)",
"HNO3" : "(ppbv)",
"Surface_Area" : "(1e-8)",
"O3loss" : "(%)",
"O3" : "(ppmv)",
"NO2" : "(ppbv)"}

PLOT_COLORS = {
"N2O" : len(PLOT_LEVELS["N2O"])-1,
"HCl" : len(PLOT_LEVELS["HCl"])-1,
"ClONO2" : len(PLOT_LEVELS["ClONO2"])-1,
"NOx" : len(PLOT_LEVELS["NOx"])-1,
"ClOx" : len(PLOT_LEVELS["ClOx"])-1,
"BrOx" : len(PLOT_LEVELS["BrOx"])-1,
"HNO3" : len(PLOT_LEVELS["HNO3"])-1,
"Surface_Area" : len(PLOT_LEVELS["Surface_Area"])-1,
"O3loss" : len(PLOT_LEVELS["O3loss"])-1,
"O3" : len(PLOT_LEVELS["O3"])-1,
"NO2" : len(PLOT_LEVELS["NO2"])-1
}

COLOR_PALETTES = {
15 : np.array([[0,221,154,6,5,5,5,5,179,251,252,251,252,251,252,251,255],
[0,141,39,3,124,204,253,253,254,253,189,126,61,3,35,75,255],
[0,252,247,249,249,251,206,8,1,7,8,26,54,115,202,251,255]]),
16 : np.array([[0,221,181,154,6,5,5,5,5,179,251,252,251,252,251,252,251,255],
[0,141,61,39,3,124,204,253,253,254,253,189,126,61,3,35,75,255],
[0,252,251,247,249,249,251,206,8,1,7,8,26,54,115,202,251,255]]),
17 : np.array([[0,221,181,154,6,5,5,5,5,179,251,252,251,252,251,251,252,252,255],
[0,141,61,39,3,124,204,253,253,254,253,189,126,61,3,3,35,35,255],
[0,252,251,247,249,249,251,206,8,1,7,8,26,54,97,115,202,202,255]]),
18 : np.array([[0,221,181,154,6,5,5,5,4,5,179,251,252,251,252,251,251,252,251,255],
[0,141,61,39,3,124,204,253,253,253,254,253,189,126,61,3,3,35,75,255],
[0,252,251,247,249,249,251,206,119,8,1,7,8,26,54,97,115,202,251,255]]),
19:np.array([[0,221,181,154,6,5,5,5,4,5,179,251,252,251,252,251,251,252,254,251,255],
[0,141,61,39,3,124,204,253,253,253,254,253,189,126,61,3,3,35,42,75,255],
[0,252,251,247,249,249,251,206,119,8,1,7,8,26,54,97,115,202,214,251,255]]),
20:np.array([[0,221,181,154,122,6,5,5,5,4,5,179,251,255,251,252,251,251,252,254,251,255],
[0,141,61,39,11,3,124,204,253,253,253,254,253,239,126,61,3,3,35,42,75,255],
[0,252,251,247,249,249,249,251,206,119,8,1,7,0,26,54,97,115,202,214,251,255]]),
21:np.array([[0,221,181,154,122,6,5,5,5,4,5,179,251,255,252,251,252,251,251,252,254,251,255],
[0,141,61,39,11,3,124,204,253,253,253,254,253,239,189,126,61,3,3,35,42,75,255],
[0,252,251,247,249,249,249,251,206,119,8,1,7,0,8,26,54,97,115,202,214,251,255]]),
}

def MODEL_post_processing(restart_file: str) -> None:
    """
    Post processing of the MODEL_history binary file into a netcdf file

    Args:
        date (str): date of the MODEL_history file
        restart_dirpath (str): directory where to find this file
    """
    fortran_file_filepath = restart_file
    restart_dirpath = os.path.dirname(restart_file)
    data_time_string = os.path.basename(fortran_file_filepath).split("_")[2]
    data_time = dt.datetime.strptime(data_time_string, "%Y%m%d%H") - dt.datetime.strptime("1970-01-01 00:00", "%Y-%m-%d %H:%M")

    with open(fortran_file_filepath, 'rb') as file:
        # Read the 6-character string
        data_string = file.read(10)
        data_integers = struct.unpack('>iiiii', file.read(20))
        pj1 = np.fromfile(file, dtype='>f8', count=180*91)
        pj1 = pj1.reshape(91, 180)
        matrices_4d = {"uj1":[], "vj1":[], "alt":[], "tj1":[]}
        for key in ["uj1","vj1","alt","tj1"]:
            matrix_4d = np.fromfile(file, dtype='>f8', count=180*91*137)
            matrix_4d = matrix_4d.reshape(137, 91, 180)
            matrices_4d[key] = matrix_4d
        qj1 = np.fromfile(file, dtype='>f8', count=180*91*137*43)
        qj1 = qj1.reshape(43, 137, 91, 180)
        surfarea = np.fromfile(file, dtype='>f8', count=180*91*137)
        surfarea = surfarea.reshape(137, 91, 180)
        hc = np.fromfile(file, dtype='>f8', count=180*91*137*15)
        hc = hc.reshape(15, 137, 91, 180)

    pression_levels_filepath = "/usr/local/REPROBUS/src/ecmwf_levels_table.csv"
    df = pd.read_csv(pression_levels_filepath)
    pression_levels = [float(elem) for elem in df["pf [hPa]"].values[1:]]

    netcdf_output_filepath = f"{restart_dirpath}/{os.path.basename(fortran_file_filepath)}.nc"

    with nc.Dataset(netcdf_output_filepath, "w", format="NETCDF4") as ncfile:
        # ------------------------------------------------------------------------------------
        ncfile.conventions              = "CF-1.0"
        ncfile.netcdf_version_id        = nc.__netcdf4libversion__
        ncfile.standard_name_vocabulary = "NetCDF Standard"
        # ------------------------------------------------------------------------------------
        ncfile.title                    = "REPROBUS (Reactive Procecesses Ruling the Ozone Budget in the Stratosphere)"
        ncfile.summary                  = "This file contains RESTART output results of the REPROBUS simulation destined for the long runs."
        ncfile.institution              = "OMP / ESPRI / Magellium"
        # ------------------------------------------------------------------------------------
        ncfile.source_filename          = os.path.basename(fortran_file_filepath)
        # ------------------------------------------------------------------------------------
        ncfile.westernmost_longitude    = 0
        ncfile.easternmost_longitude    = 360
        ncfile.southernmost_latitude    = -90
        ncfile.northernmost_latitude    = 90
        # ------------------------------------------------------------------------------------
        ncfile.createGroup("isobaric")
        ncfile.createGroup("isentropic")
        # ------------------------------------------------------------------------------------
        ncfile.createDimension("time", 1)
        ncfile.createDimension("pressure", 137)
        ncfile.createDimension("latitude", 91)
        ncfile.createDimension("longitude", 180)
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("time", "f4", ("time"))
        var.units = "days since 1970-01-01 00:00:00"
        var.standard_name = 'time'
        var.long_name = 'Data date and time'
        var.calendar = "standard"
        var[:] = data_time.days + data_time.seconds/(24*3600)
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("pressure", "f8", ("pressure"))
        var.units = "hPa"
        var.long_name = "Pressure on model ECMWF levels"
        var.standard_name = "levels_pressure"
        var[:] = pression_levels
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("latitude", "f4", ("latitude"))
        var.units = "degrees_north"
        var.long_name = "Latitude of the simulation points"
        var.standard_name = "latitude"
        var[:] = np.arange(-90,91,2)
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("longitude", "f4", ("longitude"))
        var.units = "degrees_east"
        var.long_name = "Longitude of the simulation points"
        var.standard_name = "longitude"
        var[:] = np.arange(0,360,2)
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("surface_pressure", "f8", ("time","latitude","longitude"))
        var.units = "hPa"
        var.long_name = "Surface pressure in simulation points"
        var.standard_name = "surface_pressure"
        var[:] = pj1
        # ------------------------------------------------------------------------------------
        var = ncfile.createVariable("temperature", "f8", ("time","pressure","latitude","longitude"), zlib=True, complevel=4, shuffle=True)
        var.units = "K"
        var.long_name = "Temperature on the ECMWF model levels in the simulation points"
        var.standard_name = "temperature"
        var[:] = matrices_4d["tj1"][np.newaxis, :]
        # ------------------------------------------------------------------------------------
        for ii,var_name in enumerate(SPECIES_1):
            var = ncfile.createVariable("/isobaric/"+var_name, "f8", ("time","pressure","latitude","longitude"), zlib=True, complevel=4, shuffle=True)
            var.units = "ppv"
            var.long_name = var_name
            var.standard_name = var_name
            var[:] = qj1[ii,np.newaxis,:]
        # ------------------------------------------------------------------------------------
        for ii,var_name in enumerate(SPECIES_2):
            var = ncfile.createVariable("/isobaric/"+var_name, "f8", ("time","pressure","latitude","longitude"), zlib=True, complevel=4, shuffle=True)
            var.units = "ppv"
            var.long_name = var_name
            var.standard_name = var_name
            var[:] = hc[ii,np.newaxis,:]

MODEL_post_processing("/sedoo/resos/reprobus/RESTART/MODEL_history_2023092012_001442")
MODEL_post_processing("/sedoo/resos/reprobus/RESTART/MODEL_history_2023092012_001442_my_result")
