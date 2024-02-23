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
"O3loss" : -1e-7,
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


def start_log() -> logging.Logger:
    log_handlers = [logging.StreamHandler()]
    logging.basicConfig(format="%(asctime)s   [%(levelname)s]   %(message)s",
                        datefmt="%d/%m/%Y %H:%M:%S",
                        handlers=log_handlers)
    logger = logging.getLogger('my_log')
    logger.setLevel(logging.INFO)
    return logger

def MODEL_post_processing(date: str, restart_dirpath: str) -> None:
    """
    Post processing of the MODEL_history binary file into a netcdf file

    Args:
        date (str): date of the MODEL_history file
        restart_dirpath (str): directory where to find this file
    """
    fortran_file_filepath = glob.glob(f"{restart_dirpath}/MODEL_history_{date}12_??????")[0]
    data_time_string = os.path.basename(fortran_file_filepath).split("_")[2]
    data_time = dt.datetime.strptime(data_time_string, "%Y%m%d%H") - dt.datetime.strptime("1970-01-01 00:00", "%Y-%m-%d %H:%M")

    LOGGER.info(f"Creating netCDF vesrion of {os.path.basename(fortran_file_filepath)}")
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

    pression_levels_filepath = "/home/resos/git/reprobus/table.csv"
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
            var = ncfile.createVariable("isobaric/"+var_name, "f8", ("time","pressure","latitude","longitude"), zlib=True, complevel=4, shuffle=True)
            var.units = "ppv"
            var.long_name = var_name
            var.standard_name = var_name
            var[:] = qj1[ii,np.newaxis,:]
        # ------------------------------------------------------------------------------------
        for ii,var_name in enumerate(SPECIES_2):
            var = ncfile.createVariable("isobaric/"+var_name, "f8", ("time","pressure","latitude","longitude"), zlib=True, complevel=4, shuffle=True)
            var.units = "ppv"
            var.long_name = var_name
            var.standard_name = var_name
            var[:] = hc[ii,np.newaxis,:]

def stations_post_processing(date: str, result_dirpath: str) -> None:
    """

    This function transforms stations ASCII files into netcdf files
    ('stations_' and 'reprobus_[station_name]' files)

    Args:
        date (str): date of the files to use
        result_dirpath (str): directory where to find these files
    """
    # *********************************************************************************************
    files = glob.glob(f"{result_dirpath}/reprobus_*_{date}12_??????")
    for input_filename in files:
        LOGGER.info(f"Creating netCDF version of {os.path.basename(input_filename)}")
        output_filename = f"{result_dirpath}/{os.path.basename(input_filename)}.nc"
        with open(input_filename, "r") as f:
            lines = f.readlines()
        parts = [i for i, x in enumerate(lines) if x == lines[0]]
        start_index = parts
        end_index = parts[1:] + [len(lines)]
        with nc.Dataset(output_filename, "w", format="NETCDF4") as ds:
            ds.conventions              = "CF-1.0"
            ds.netcdf_version_id        = nc.__netcdf4libversion__
            ds.standard_name_vocabulary = "NetCDF Standard"
            ds.title = f"REPROBUS results from the {os.path.basename(input_filename)} file"
            ds.run = lines[0].split(" ")[-1].strip()
            # ds.julian_day = lines[2].split("=")[-1].strip()
            # ds.datetime = dt.datetime.strftime(dt.datetime.strptime(lines[3].strip(), "%d %B %Y %H 0UT"), "%d/%m/%Y %H:%M")
            ds.station_name = lines[5].strip()
            ds.station_latitude = lines[7].split("LATITUDE")[-1].split("=")[-1].strip()
            ds.station_longitude = lines[7].split("LATITUDE")[0].split("=")[-1].strip()
            ds.latitude = lines[8].split("ILAT")[-1].split("=")[-1].strip()
            ds.longitude = lines[8].split("ILAT")[0].split("=")[-1].strip()
            ds.solar_zenith_angle = lines[10].split("=")[-1].strip()
            ds.createDimension("time", len(parts))
            ds.createDimension("altitude", 137)
            ds.createVariable("time", "f4", ("time"))
            ds.createVariable("altitude", "f4", ("altitude"))
            ds.createVariable("pressure", "f4", ("altitude"))
            ds.createVariable("temperature", "f4", ("altitude"))
            ds.createVariable("theta", "f4", ("altitude"))
            ds.createVariable("density", "f4", ("altitude"))
            for varname in ["N2O","CH4","H2O","NOy","HNO3","N2O5"]:
                ds.createVariable(varname, "f4", ("time","altitude"))
            for ii in range(len(parts)):
                part = lines[start_index[ii]:end_index[ii]]
                alt = [float(line.strip().split(" ")[0]) for line in part[14:(14+137)]]
                pressure = [float(line.strip().split(" ")[2]) for line in part[14:(14+137)]]
                temperature = [float(line.strip().split(" ")[4]) for line in part[14:(14+137)]]
                theta = [float(line.strip().split(" ")[6]) for line in part[14:(14+137)]]
                density = [float(line.strip().split(" ")[8]) for line in part[14:(14+137)]]
                N2O = [float(line.strip().split(" ")[10]) for line in part[14:(14+137)]]
                ds["time"][ii] = float(part[2].split("=")[-1].strip())
                if ii == 0:
                    ds["altitude"][:] = alt
                    ds["pressure"][:] = pressure
                    ds["temperature"][:] = temperature
                    ds["theta"][:] = theta
                    ds["density"][:] = density
                ds["N2O"][ii,:] = N2O
                for jj,var in enumerate(["CH4","H2O","NOy","HNO3","N2O5"]):
                    arr = [float(line.strip().split(" ")[(jj+1)*2]) for line in part[154:(154+137)]]
                    ds[var][ii,:] = arr
    # *********************************************************************************************
    files = glob.glob(f"{result_dirpath}/stations_{date}??_??????")
    for input_file in files:
        LOGGER.info(f"Creating netCDF version of {os.path.basename(input_file)}")
        df = pd.DataFrame({"Station name": [],
                        "Julian day": [],
                        "Date": [],
                        "O3_1": [],
                        "O3_2": [],
                        "NO2": [],
                        "OClO": []})
        df_res = pd.DataFrame()
        with open(input_file, "r") as f:
            lines = f.readlines()
        df["Station name"] = [line[:13].strip() for line in lines]
        df["Julian day"] = [float(line[13:19].strip()) for line in lines]
        df["Date"] = pd.to_datetime(df["Julian day"], unit='D', origin=f"{os.path.basename(input_file).split('_')[1][:4]}-01-01")
        df["O3_1"] = [float(line[20:26].strip()) for line in lines]
        df["O3_2"] = [float(line[30:35].strip()) for line in lines]
        df["NO2"] = [float(line[41:51].strip()) for line in lines]
        df["OClO"] = [float(line[62:72].strip()) for line in lines]
        df_res = pd.concat((df_res, df))
        start_date = dt.datetime.strftime(df_res["Date"].min().date(), "%Y%m%d")
        end_date   = dt.datetime.strftime(df_res["Date"].max().date(), "%Y%m%d")
        df_res = df_res.sort_values(by=["Station name", "Julian day"])
        stations = df_res["Station name"].unique()
        times    = df_res["Julian day"].unique()
        output_filename = f"{result_dirpath}/${os.path.basename(input_file)}.nc"
        with nc.Dataset(output_filename, "w", format="NETCDF4") as ds:
            ds.conventions              = "CF-1.0"
            ds.netcdf_version_id        = nc.__netcdf4libversion__
            ds.standard_name_vocabulary = "NetCDF Standard"
            ds.title = f"REPROBUS results from the {os.path.basename(input_file)} file"
            ds.createDimension("time", len(times))
            ds.createDimension("species", 4)
            ds.createVariable("time", "f4", ("time"))
            ds["time"][:] = times
            gr = ds.createGroup("O3_1")
            gr.description = "Integrated O3 column"
            gr.units = "DU"
            gr = ds.createGroup("O3_2")
            gr.description = "Integrated passive O3 column"
            gr.units = "DU"
            gr = ds.createGroup("NO2")
            gr.description = "Integrated NO2"
            gr.units = "mol.cm-2"
            gr = ds.createGroup("OClO")
            gr.description = "Integrated OClO"
            gr.units = "mol.cm-2"
            for varname in stations:
                ds["O3_1"].createVariable(varname, "f4", ("time"))
                ds["O3_2"].createVariable(varname, "f4", ("time"))
                ds["NO2"].createVariable(varname, "f4", ("time"))
                ds["OClO"].createVariable(varname, "f4", ("time"))
                ds["O3_1/"+varname][:] = df_res.loc[df_res["Station name"]==varname]["O3_1"]
                ds["O3_2/"+varname][:] = df_res.loc[df_res["Station name"]==varname]["O3_2"]
                ds["NO2/"+varname][:] = df_res.loc[df_res["Station name"]==varname]["NO2"]
                ds["OClO/"+varname][:] = df_res.loc[df_res["Station name"]==varname]["OClO"]

def compute_on_theta_levels(date: str, restart_dirpath: str) -> xr.Dataset:
    # =============================================================================
    # COEFFICIENTS A ET B DES NIVEAUX HYBRIDES
    # =============================================================================
    LOGGER.info("Computing a and b coefficients")
    df = pd.read_csv("/usr/local/REPROBUS/src/ecmwf_137_levels.txt",
                     sep="\t",
                     skiprows=2,
                     names = ["n","a","b","ph[hPa]","pf[hPa]"])
    aaa = (df["a"].values[:-1] + df["a"].values[1:])*0.5*0.01
    bbb = (df["b"].values[:-1] + df["b"].values[1:])*0.5
    # =============================================================================
    # ISOTHETA ; ISENTROPE CHOISIE POUR LE TRACE
    # =============================================================================
    LOGGER.info("Reading data")
    nbcon = 44
    ncm = 15
    fortran_file_filepath = glob.glob(f"{restart_dirpath}/MODEL_history_{date}12_??????")[0]
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
        qj1 = np.fromfile(file, dtype='>f8', count=180*91*137*nbcon)
        qj1 = qj1.reshape(nbcon, 137, 91, 180)
        hc = np.fromfile(file, dtype='>f8', count=180*91*137*ncm)
        hc = hc.reshape(ncm, 137, 91, 180)
    ds = xr.Dataset(data_vars={"pj1":(["lat", "lon"], pj1),
                               "uj1":(["niv", "lat", "lon"], matrices_4d["uj1"]),
                               "vj1":(["niv", "lat", "lon"], matrices_4d["vj1"]),
                               "alt":(["niv", "lat", "lon"], matrices_4d["alt"]),
                               "tj1":(["niv", "lat", "lon"], matrices_4d["tj1"]),
                               "qj1":(["dim1", "niv", "lat", "lon"], qj1),
                               "hc":(["dim2", "niv", "lat", "lon"], hc),
                               "aaa":(["niv"], aaa),
                               "bbb":(["niv"], bbb)})
    # =============================================================================
    #   BOUCLER SUR THETA
    # =============================================================================
    nlat = 91
    nlon = 180
    nlonp1 = nlon+1
    lon = np.arange(nlonp1)*362/nlonp1
    lat = 90 - np.arange(nlat)*180/(nlat-1)
    
    rter  = 6.371229e+06
    gg    = 9.80665
    omega = 7.292e-05
    rascp = 2./7.
    p0    = 1000.

    theta_arr = [435., 475., 550.]

    nlong = 9
    nshort = 2
    il_arr = np.array([1,13,14,21,23,24,43,44,11]) - 1
    il_names = ["N2O","HCl","ClONO2","NOx","ClOx","BrOx","HNO3","Surface_Area","O3loss"]
    is_arr = np.array([5, 7])-1
    is_names = ["O3","NO2"]

    output_ds = xr.Dataset(coords={"theta":(["theta"], theta_arr),
                                    "lat":(["lat"], lat),
                                    "lon":(["lon"], lon)})
    species_da = {}
    for i_theta, theta_value in enumerate(theta_arr):
        LOGGER.info(f"Computing variables for theta {theta_value} K")
        # =============================================================================
        # Calcul de theta à chaque point de grille
        # =============================================================================
        LOGGER.info(" -- Calcul de theta")
        pmb = ds["aaa"] + ds["pj1"]*ds["bbb"]
        theta = ds["tj1"]*((p0/pmb)**rascp)
        # =============================================================================
        # Detection des niveaux encadrant isotheta, calcul des poids
        # =============================================================================
        LOGGER.info(" -- Calcul des niveaux")
        index = (theta<theta_value).argmax(dim="niv")
        cinf  = (theta_value - theta.isel(niv=index)) / (theta.isel(niv=index-1) - theta.isel(niv=index))
        csup  = 1 - cinf
        # =============================================================================
        # Calcul de dthetadp aux niveaux qui encadrent isotheta,
        # puis interpolation verticale
        # =============================================================================
        dthetadpsup = (theta.isel(niv=index-2) - theta.isel(niv=index)) / ((pmb.isel(niv=index-2) - pmb.isel(niv=index))*100)
        dthetadpinf = (theta.isel(niv=index-1) - theta.isel(niv=index+1)) / ((pmb.isel(niv=index-1) - pmb.isel(niv=index+1))*100)
        dthetadp    = dthetadpsup*cinf + dthetadpinf*csup
        # =============================================================================
        # Interpolation verticale des vents sur la surface isotheta
        # =============================================================================
        LOGGER.info(" -- Interpolation verticale")
        utheta = cinf*ds["uj1"].sel(niv=index-1) + csup*ds["uj1"].sel(niv=index)
        vtheta = cinf*ds["vj1"].sel(niv=index-1) + csup*ds["vj1"].sel(niv=index)
        # =============================================================================
        # Calcul du tourbillon sur la surface isotheta
        # =============================================================================
        LOGGER.info(" -- Calcul du tourbillon")
        xpi   = 2.*math.asin(1.)
        xpih  = xpi/2.
        xpi2  = xpi*2.
        ddlat = xpi/(nlat - 1)
        ddlon = xpi2/(nlon - 1)
        
        dlat     = xpih - (np.arange(1,nlat+1)-1)*ddlat
        cosphi   = np.cos(dlat)
        coriolis = 2 * omega * np.sin(dlat)
        
        dvdlam, dudphi, eta = xr.DataArray(np.full((nlat, nlon), np.nan), dims=["lat","lon"]), \
                            xr.DataArray(np.full((nlat, nlon), np.nan), dims=["lat","lon"]), \
                            xr.DataArray(np.full((nlat, nlon), np.nan), dims=["lat","lon"])
        dvdlam[1:-1,0]      = (vtheta[1:-1,1] - vtheta[1:-1,-1]) / (2*ddlon)
        dvdlam[1:-1,-1]     = (vtheta[1:-1,0] - vtheta[1:-1,-2]) / (2*ddlon)
        dvdlam[1:-1,1:-1]   = (vtheta[1:-1,2:] - vtheta[1:-1,0:-2]) / (2*ddlon)
        dudphi[1:-1,:]      = (utheta[0:-2,:]*cosphi[0:-2, np.newaxis] - utheta[2:,:]*cosphi[2:, np.newaxis]) / (2*ddlat)
        eta[1:-1,:]         = (dvdlam[1:-1,:] - dudphi[1:-1,:]) / (rter * cosphi[1:-1, np.newaxis])
        # =============================================================================
        # Calcul du tourbillon potentiel proprement dit
        # =============================================================================
        pv = -gg * (eta + coriolis[:,np.newaxis]) * dthetadp
        # =============================================================================
        # Traitement des poles
        # =============================================================================
        LOGGER.info(" -- Traitement des pôles")
        pvpolen  = np.sum(pv[1,:]/nlon)
        pvpoles  = np.sum(pv[-2,:]/nlon)
        pv[0,:]  = pvpolen
        pv[-1,:] = pvpoles
        # =============================================================================
        # Temperature
        # =============================================================================
        LOGGER.info(" -- Temperature")
        num_id = 1
        gridfi_t = cinf*ds["tj1"].sel(niv=index-1) + csup*ds["tj1"].sel(niv=index)
        gridfi_t = xr.concat([gridfi_t, gridfi_t[:,0]], dim="lon").expand_dims("theta")
        if i_theta==0:
            gridfi_t_da = gridfi_t
        else:
            gridfi_t_da = xr.concat([gridfi_t_da, gridfi_t], dim="theta")
        # =============================================================================
        # PV
        # =============================================================================
        LOGGER.info(" -- PV")
        num_id = 2
        gridfi_pv = xr.concat([pv, pv[:,0]], dim="lon").expand_dims("theta")
        if i_theta==0:
            gridfi_pv_da = gridfi_pv
        else:
            gridfi_pv_da = xr.concat([gridfi_pv_da, gridfi_pv], dim="theta")
        # =============================================================================
        # Especes transportees
        # =============================================================================
        for ii in range(nlong):
            LOGGER.info(f" -- Espèce transporté {ii+1}/{nlong}")
            nc = il_arr[ii]
            num_id = nc+3
            if nc==11:
                gridfi_espece = ((cinf*ds["qj1"].sel(niv=index-1, dim1=8) + csup*ds["qj1"].sel(niv=index, dim1=8)) - \
                                (cinf*ds["qj1"].sel(niv=index-1, dim1=11) + csup*ds["qj1"].sel(niv=index, dim1=11))) / \
                                (cinf*ds["qj1"].sel(niv=index-1, dim1=11) + csup*ds["qj1"].sel(niv=index, dim1=11))
            elif nc==23:
                gridfi_espece = cinf*ds["qj1"].sel(niv=index-1, dim1=nc) + \
                                csup*ds["qj1"].sel(niv=index, dim1=nc) + \
                                2*(cinf*ds["qj1"].sel(niv=index-1, dim1=25)) + \
                                csup*ds["qj1"].sel(niv=index, dim1=25)                              
            else:
                gridfi_espece = cinf*ds["qj1"].sel(niv=index-1, dim1=nc) + csup*ds["qj1"].sel(niv=index, dim1=nc)
            gridfi_espece = xr.concat([gridfi_espece, gridfi_espece[:,0]], dim="lon").expand_dims("theta")
            if i_theta==0:
                species_da[il_names[ii]] = gridfi_espece
            else:
                species_da[il_names[ii]] = xr.concat([species_da[il_names[ii]], gridfi_espece], dim="theta")
        # =============================================================================
        # Especes a l'equilibre
        # =============================================================================
        for ii in range(nshort):
            LOGGER.info(f" -- Espèce à l'équilibre {ii+1}/{nshort}")
            nc = is_arr[ii]
            num_id = nbcon + nc + 3
            gridfi_espece = cinf*ds["hc"].sel(niv=index-1, dim2=nc) + csup*ds["hc"].sel(niv=index, dim2=nc)
            gridfi_espece = xr.concat([gridfi_espece, gridfi_espece[:,0]], dim="lon").expand_dims("theta")
            if i_theta==0:
                species_da[is_names[ii]] = gridfi_espece
            else:
                species_da[is_names[ii]] = xr.concat([species_da[is_names[ii]], gridfi_espece], dim="theta")
    output_ds = output_ds.assign({"t":gridfi_t_da, "pv":gridfi_pv_da}).assign(species_da)
    return output_ds

def create_theta_plots(dataset: xr.Dataset, im_dir: str) -> None:
    # vars_to_plot = ["N2O","HCl","ClONO2","NOx","ClOx","BrOx","HNO3","Surface_Area","O3loss", "O3","NO2"]
    vars_to_plot = ["O3loss"]
    theta_arr = dataset["theta"].values
    for var in vars_to_plot:
        cmap = COLOR_PALETTES[PLOT_COLORS[var]].T/255
        cmap /= cmap.max()
        cmap = [tuple(line) for line in cmap[1:-1,:]]
        custom_cmap = ListedColormap(cmap)
        for ii,theta_val in enumerate(theta_arr):
            LOGGER.info(f"Creating figure {var} on {theta_val} K")
            # LOGGER.info((dataset[var][ii,:,:]/PLOT_COEFFS[var]).values)
            fig = plt.figure()
            p = (dataset[var][ii,:,:]/PLOT_COEFFS[var]).plot.contourf(
                                                                    transform=ccrs.PlateCarree(),
                                                                    subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                                                    cmap=custom_cmap,
                                                                    levels=np.array(PLOT_LEVELS[var]))
            p1 = (dataset[var][ii,:,:]/PLOT_COEFFS[var]).plot.contour(
                                                                    transform=ccrs.PlateCarree(),
                                                                    subplot_kws={"projection": ccrs.SouthPolarStereo()},
                                                                    colors="k", linestyles="-", linewidths=0.5,
                                                                    levels=np.array(PLOT_LEVELS[var]))
            p.axes.coastlines(color="w", linewidth=1.5)
            obj = p.axes.gridlines(linestyle="--", linewidth=0.5)
            p.axes.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
            title = f"Reprobus : {var} {PLOT_UNITS[var]}\n{theta_val} K"
            p.axes.set_title(title, fontsize=15)
            p.colorbar.ax.yaxis.label.set_fontsize(15)
            p.colorbar.set_ticks(np.array(PLOT_LEVELS[var]))
            p.colorbar.ax.tick_params(labelsize=10)
            fig.savefig(f"{im_dir}/{var}_{int(theta_val)}.png", dpi=200, bbox_inches='tight')
            plt.close(fig)

def create_figures(date: str, restart_dirpath: str, im_dir: str) -> None:
    # ds = compute_on_theta_levels(date, restart_dirpath)
    # ds.to_netcdf(f"{im_dir}/image_data.nc")
    ds = xr.open_dataset(f"{im_dir}/image_data.nc")
    create_theta_plots(ds, im_dir)


if __name__=="__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Post-processing of the REPROBUS Fortran output files",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--date", type=str, help="End date of the simulation in YYMMDD format")
    parser.add_argument("--restart-dir", type=str, help="Path to the directory with the MODEL_history Fortran binary files")
    parser.add_argument("--res-dir",  type=str, help="Path to the directory where the station result files are stored")
    parser.add_argument("--image-dir", type=str, help="Path to the directory where to save the figures")
    args = parser.parse_args()

    global LOGGER
    LOGGER = start_log()

    LOGGER.info("Starting post-processing of the REPROBUS output")
    # MODEL_post_processing(args.date, args.restart_dir)
    # stations_post_processing(args.date, args.res_dir)
    create_figures(args.date, args.restart_dir, args.image_dir)

