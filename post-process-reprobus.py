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

SPECIES_1 = ['N2O','CH4','H2O','NOy','HNO3','N2O5','Cly','Ox','CO','OClO','Passive Ox','H2SO4','HCl','ClONO2','HOCl','Cl2','H2O2','ClNO2','HBr','BrONO2','NOx','HNO4','ClOx','BrOx','Cl2O2','HOBr','BrCl','CH2O','CH3O2','CH3O2H','CFC-11','CFC-12','CFC-113','CCl4','CH3CCl3*','CH3Cl','HCFC-22*','CH3Br','H-1211*','H-1301','Bry','CH2Br2*','HNO3 GAS']
SPECIES_2 = ['O(1D)','OH','Cl','O(3P)','O3','HO2','NO2','NO','Br','N','ClO','BrO','NO3','H','CH3']

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


if __name__=="__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Post-processing of the REPROBUS Fortran output files",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--date", type=str, help="End date of the simulation in YYMMDD format")
    parser.add_argument("--restart-dir", type=str, help="Path to the directory with the MODEL_history Fortran binary files")
    parser.add_argument("--res-dir",  type=str, help="Path to the directory where the station result files are stored")
    args = parser.parse_args()

    global LOGGER
    LOGGER = start_log()

    LOGGER.info("Starting post-processing of the REPROBUS output")
    MODEL_post_processing(args.date, args.restart_dir)
    stations_post_processing(args.date, args.res_dir)
