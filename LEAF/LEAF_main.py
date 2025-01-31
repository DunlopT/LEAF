# LEAF main code
'''
Copyright notice

This script comprises the main code for the LEAF model version 1.0 (dated 31.01.2025) as outlined in the following publication:
Dunlop, T., Felder, S., Glamore, W. 2025. A mangrove Lifecycle Ecosystem Analysis and Forecasting (LEAF) model. Environmental Modelling & Software.

The LEAF model and associated input files were developed as part of the PhD research by Thomas Dunlop conducted at UNSW Sydney.
Correspondence regarding the model should be directed to:

    Thomas Dunlop
    t.dunlop@unsw.edu.au
    Water Research Laboratory
    School of Civil and Environmental Engineering, UNSW Sydney
    110 King St, Manly Vale NSW 2093 Australia

This code has been developed such that key logical expressions and input values can be updated based on the latest available data from the field and laboratory.
This code can be modified to tailor the code to the user's project site.
The LEAF model has been developed to enable practitioners and researchers without detailed Python knowledge to use this tool to track mangrove ecosystem trajectories.

The required inputs for the main code are listed in the sections "FOLDER SETUP" and "INPUTS" below.
Where logical expressions can be updated in conjunction with the separate input file, the following comment is made next to the code: "NOTE: OPTIONAL USER UPDATE".

For simplicity in LEAF version 1.0, all lifecycle calculations are conducted in this code.
'''
########################################################################
# MODEL SETUP
########################################################################
# ======================================================================
# PACKAGE IMPORT
# ======================================================================
# Standard library imports
import numpy as np
import os
import json
import math
import ctypes
import matplotlib.pyplot as plt
import datetime
import time
import pandas
import netCDF4 as nc
import warnings
import faulthandler
faulthandler.enable()
from matplotlib.colors import Normalize
from datetime import datetime, timedelta, date
from mdu import read_deltares_ini
from dateutil import parser
from math import floor
# Third party imports
from bmi.wrapper import BMIWrapper # Used to connect to and run D3D

# Matplotlib error prevention line
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"

# Suppress plot warnings
warnings.filterwarnings("ignore", category=UserWarning, message=".*Setting the 'color' property will override the edgecolor or facecolor properties.*")

# For labelling output files based on current date and time
now = time.localtime()
today = time.strftime("%Y-%m-%d_%H-%M", now) # Format of YYYY-MM-DD_HH-MM

# ======================================================================
# FOLDER SETUP
# ======================================================================
# Hydro-morphodynamic folder setup
# LEAF has been coupled to Delft3D Flexible Mesh (DFM)
# These are the paths to the main Delft3D Flexible Mesh folders and dll files
D3D_HOME = os.path.join('C:\\Program Files (x86)','Deltares','Delft3D Flexible Mesh Suite HMWQ (2021.03)','plugins','DeltaShell.Dimr','kernels','x64')                      
dflowfm_path = os.path.join(D3D_HOME,'dflowfm','bin','dflowfm.dll')
dwaves_path = os.path.join(D3D_HOME,'dwaves','bin','wave.dll')
dimr_path = os.path.join(D3D_HOME,'dimr','bin','dimr_dll.dll')

# These are the paths to the working directory and associated DFM files, i.e., where the DFM project is saved
workdir = os.path.join('C:\\Users','Name','Documents','DFM MODEL RUNS','PROJECTNAME') # Update filepath for the folder in which the DFM project is saved
config_file = os.path.join(workdir,'dimr_config.xml') # Update filepath for the DFM configuration file
mdu_file = os.path.join(workdir,'dflowfm','FlowFM.mdu') # Update filepath for MDU file from DFM
mdw_file = os.path.join(workdir,'wave','Waves.mdw') # Update filepath for MDW file from DFM (for coupled D-Waves model)
grid_file = os.path.join(workdir,'dflowfm','FlowFM_net.nc') # Update filepath for grid file
mor_file = os.path.join(workdir,'dflowfm','FlowFM.mor') # Update filepath for morphology file
wl_file = os.path.join(workdir,'dflowfm','WaterLevel.bc') # Update filepath for water level file
map_file = os.path.join(workdir,'PROJECTNAME.dsproj_data','FlowFM','output','FlowFM_map.nc') # Update filepath for map file
# This is the filepath to the inputs file for the LEAF model
input_file = os.path.join(workdir,'dflowfm','LEAF_inputs.json') # User defined inputs file

# Paths to sub-folders for saving csv files and figures
csvpath = f'{today}_csv' # Name of folder to store csv files
plotpath = f'{today}_plots' # Name of folder to store plots / figures
csvfolder_path = os.path.join(workdir,'dflowfm', csvpath) # Path to csv sub-folder
plotfolder_path = os.path.join(workdir,'dflowfm',plotpath) # Path to plot sub-folder
# Create the subfolder if it doesn't exist
os.makedirs(csvfolder_path, exist_ok=True)
os.makedirs(plotfolder_path, exist_ok=True)

# ======================================================================
# ENVIRONMENT AND MODEL INITIALISATION
# ======================================================================
# Update the environment variable PATH in the BMI wrapper python file
os.environ['PATH'] = os.path.join(D3D_HOME,'share','bin')\
+ ";" + os.path.join(D3D_HOME,'dflowfm','bin')\
+ ";" + os.path.join(D3D_HOME,'dimr','bin') \
+ ";" + os.path.join(D3D_HOME,'dwaves','bin') \
+ ";" + os.path.join(D3D_HOME,'esmf','scripts') \
+ ";" + os.path.join(D3D_HOME,'swan','scripts')

# Define DFM dll and DIMR paths using ctypes.WinDLL function and setting winmode to 0. Bypass for a bug in the ctypes package init module.
_func1 = ctypes.WinDLL('C:\\Program Files (x86)\\Deltares\\Delft3D Flexible Mesh Suite HMWQ (2021.03)\\plugins\\DeltaShell.Dimr\\kernels\\x64\\dflowfm\\bin\\dflowfm.dll', winmode = 0)
_func2 = ctypes.WinDLL('C:\\Program Files (x86)\\Deltares\\Delft3D Flexible Mesh Suite HMWQ (2021.03)\\plugins\\DeltaShell.Dimr\\kernels\\x64\\dimr\\bin\\dimr_dll.dll', winmode = 0)

# Define DFM wrapper
model_dfm = BMIWrapper(engine=dflowfm_path,configfile=mdu_file)

# Define and initalise DIMR wrapper
model_dimr = BMIWrapper(engine=dimr_path,configfile=config_file)
model_dimr.initialize()
print('model initialised')

# ======================================================================
# BMI PARAMETER EXCHANGE
# ======================================================================
# Get pointers to important model variables from DFM (see parameter list here:
# https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/dflowfm/packages/dflowfm_lib/include/bmi_get_var.inc
bedlevel = model_dfm.get_var('bl') # Bed level retrieved here to inform bed level for initial vegetation (if included)
bedlevel = np.array(bedlevel)

# Get pointer for surface area
ba = model_dfm.get_var('ba') # surface area of the boxes (bottom area) {"location": "face", "shape": ["ndx"]}
ba = np.array(ba)

# Retrieve pointers for vegetation (mangrove) parameters to update the model
stemheight = model_dfm.get_var('stemheight') # Retrieve stem height
pandas.DataFrame(stemheight).to_csv(os.path.join(csvfolder_path,'initial_stemheight.csv'), index=False) # output stemheight to csv
stemheight = np.array(stemheight) # Create numpy array for stemheight
rnveg = model_dfm.get_var('rnveg') # Retrieve vegetation density
rnveg = np.array(rnveg) # Create numpy array for vegetation density
print('Length of vegetation parameter arrays: ', len(rnveg))
diaveg = model_dfm.get_var('diaveg') # Retrieve stem diameter
diaveg = np.array(diaveg) # Create numpy array for stem diameter

# Retrieve number of grid cells 
ndx = model_dfm.get_var('ndx') # no. of cells including boundary cells. Note that this is trimmed when plotting
ndxi = model_dfm.get_var('ndxi') # no. of cells excluding boundary cells

cdvegsp = model_dfm.get_var('Cdvegsp')
cdvegsp = np.array(cdvegsp) # Create numpy array for vegetation drag coefficient
cdvegsp_updated = cdvegsp[:ndx]

# Retrieve morphological factor
getmor = read_deltares_ini(mor_file)
morfac = float(getmor[(getmor['section'] == 'Morphology') & (getmor['key'] == 'MorFac')]['value'].values[0])
print('morfac:',morfac)

## Obtain start and end dates from the model
getmdu = read_deltares_ini(mdu_file) # This refers to the mdu file as specified above
refdate = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'RefDate')]['value']
tstart = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'TStart')]['value']
tend = getmdu[(getmdu['section'] == 'time') & (getmdu['key'] == 'TStop')]['value']

# ======================================================================
# INPUTS
# ======================================================================
# Outside of the JSON input file, the following inputs are required
# Plotting parameters
xsecycoords = [112.5, 137.5] # Update y-coordinates of the cross-sections that are plotted
xplotminbl = -1.0 # Update the minimum bed level for cross-section plot
includestatus = 'Y' # Y/N. Option to colour code mangrove canopies in the cross-section plot based on their health status. Default value is set to Y
includeages = 'N' # Y/N. Option to include stem age in the canopies in the cross-section plot. Default value is set to N

# The following section loads the inputs from the JSON input file. No action here is required unless additional inputs are added or inputs are removed
# Read parameters from the JSON input file
with open(input_file,'r') as f:
    parameters = json.load(f)

# Define parameters from the JSON file
# NOTE: OPTIONAL USER UPDATE - Update this list with any additional parameters added to the JSON input file
# 1. Time
mtpervt = parameters['1. Time']['mtpervt'] # 1A. Timestep duration (s)
nt = parameters['1. Time']['nt'] # 1B. Number of timesteps (-)
seedwinstart = datetime.strptime(parameters['1. Time']['seedwinstart'], '%d/%m/%Y') # 1C. Fruiting window start date (DD/MM/YYYY)
seedwinend = datetime.strptime(parameters['1. Time']['seedwinend'], '%d/%m/%Y') # 1D. Fruiting window end date (DD/MM/YYYY)
seedwinrep = parameters['1. Time']['seedwinrep'] # 1E. Fruiting window annual recurrence (Y/N)
# 2. Plotting
plotint = parameters['2. Plotting']['plotint'] # 2A. Plotting interval (Number of timesteps)
# 3. Initial mangroves
iniveg = parameters['3. Initial mangroves']['iniveg'] # 3A. Initial mangrove presence (Y/N)
iniden = parameters['3. Initial mangroves']['iniden'] # 3B. Stem density (stems/m2)
inistemdia = parameters['3. Initial mangroves']['inistemdia'] # 3C. Stem diameter (m)
inirtlen = parameters['3. Initial mangroves']['inirtlen'] # 3D. Cable root length (m)
inirtdia = parameters['3. Initial mangroves']['inirtdia'] # 3E. Cable root diameter (m)
inipnelen = parameters['3. Initial mangroves']['inipnelen'] # 3F. Pneumatophore length (m)
inipnedia = parameters['3. Initial mangroves']['inipnedia'] # 3G. Pneumatophore diameter (m)
iniage = parameters['3. Initial mangroves']['iniage'] # 3H. Age (years)
# 4. Mangrove extrema
maxstemsmsq = parameters['4. Mangrove extrema']['maxstemsmsq'] # 4A. Maximum stem density (stems/m2)
maxstemht = parameters['4. Mangrove extrema']['maxstemht'] # 4B. Maximum stem height (m)
maxstemdia = parameters['4. Mangrove extrema']['maxstemdia'] # 4C. Maximum stem diameter (m)
maxrtlen = parameters['4. Mangrove extrema']['maxrtlen'] # 4D. Maximum cable root length (m)
maxrtdia = parameters['4. Mangrove extrema']['maxrtdia'] # 4E. Maximum cable root diameter (m)
maxrts = parameters['4. Mangrove extrema']['maxrts'] # 4F. Maximum number of cable roots (-)
maxpnedia = parameters['4. Mangrove extrema']['maxpnedia'] # 4G. Maximum pneumatophore diameter (m)
# 5. Other mangrove parameters
sapht = parameters['5. Other mangrove parameters']['sapht'] # 5A. Sapling height threshold (m)
fecundity = parameters['5. Other mangrove parameters']['fecundity'] # 5B. Fecundity age (years)
cablertdepth = parameters['5. Other mangrove parameters']['cablertdepth'] # 5C. Cable root depth (m)
pneumint = parameters['5. Other mangrove parameters']['pneumint'] # 5D. Pneumatophore spacing (m)
autoch = parameters['5. Other mangrove parameters']['autoch'] # 5E. Autochthonous accretion rate (m/year)
startrt = parameters['5. Other mangrove parameters']['startrt'] # 5F. Starting cable/tap root length on establishment (m)
startdia = parameters['5. Other mangrove parameters']['startdia'] # 5G. Starting stem diameter on establishment (m)
startpnedia = parameters['5. Other mangrove parameters']['startpnedia'] # 5H. Starting pneumatophore diameter on growth (m)
# 6. Establishment
estactive = parameters['6. Establishment']['estactive'] # 6A. Establishment stage activation (Y/N)
seedtrvldist = parameters['6. Establishment']['seedtrvldist'] # 6B. Seed travel distance (m)
inundfreepd = parameters['6. Establishment']['inundfreepd'] # 6C. Inundation free period (hrs)
estchance = parameters['6. Establishment']['estchance'] # 6D. Chance of establishment (-)
taucrit = parameters['6. Establishment']['taucrit'] # 6E. Critical bed shear stress (N/m2 or formula)
vegtrkrt_thresh = parameters['6. Establishment']['vegtrkrt_thresh'] # 6F. Cable root length requirement for stem growth (m or formula)
erocrit = parameters['6. Establishment']['erocrit'] # 6G. Critical erosion limit (m or formula)
# 7. Growth
growthactive = parameters['7. Growth']['growthactive'] # 7A. Growth stage activation (Y/N)
growthlogic = parameters['7. Growth']['growthlogic'] # 7B. Growth logic (main code to be adjusted below)
vegtrkdiagrowth = parameters['7. Growth']['vegtrkdiagrowth'] # 7C. Seedling stem diameter growth rate (m/day or formula)
seedht2dia = parameters['7. Growth']['seedht2dia'] # 7D. Seedling stem height to diameter relationship (number or formula)
vegtrkrt_seedling_growth = parameters['7. Growth']['vegtrkrt_seedling_growth'] # 7E. Seedling cable root growth rate (m/day or formula)
seedrtdia2len = parameters['7. Growth']['seedrtdia2len'] # 7F. Seedling root diameter to length relationship (number or formula)
vegtrkdia_adult_growth = parameters['7. Growth']['vegtrkdia_adult_growth'] # 7G. Sapling/adult stem diameter growth rate (m/day or formula)
dbh2ht = parameters['7. Growth']['dbh2ht'] # 7H. Sapling/adult stem height to diameter relationship (number or formula)
vegtrkrt_adult_growth = parameters['7. Growth']['vegtrkrt_adult_growth'] # 7I. Sapling/adult cable root growth rate
vegtrkrtdia_adult_growth = parameters['7. Growth']['vegtrkrtdia_adult_growth'] # 7J. Sapling/adult cable root diameter to length relationship (number or formula)
pnelenreq = parameters['7. Growth']['pnelenreq'] # 7K. Pneumatophore growth condition (proportion of pneumatophore height)
vegtrkpnedia_adult_growth = parameters['7. Growth']['vegtrkpnedia_adult_growth'] # 7L. Pneumatophore base diameter growth rate (m/day or formula)
pnelen2pnedia = parameters['7. Growth']['pnelen2pnedia'] # 7M. Ratio of pneumatophore height to base diameter (number or formula)
# 8. Recovery and mortality
# 8A. General
mortactive = parameters['8. Recovery and mortality']['mortactive'] # 8A1. Recovery and mortality stage activation (Y/N)
# 8B. Inundation
inundactive = parameters['8. Recovery and mortality']['inundactive'] # 8B1. Inundation mortality activation (Y/N)
inundhtthresh = parameters['8. Recovery and mortality']['inundhtthresh'] # 8B2. Minimum depth of inundation to cause mortality (m)
inundstemht = parameters['8. Recovery and mortality']['inundstemht'] # 8B3. Proportion of seedling stem height inundated to cause stress/mortality (-)
inundpneumlen = parameters['8. Recovery and mortality']['inundpneumlen'] # 8B4. Proportion of sapling/adult pneumatophore height inundated to cause stress/mortality (-)
inundtimeseed1 = parameters['8. Recovery and mortality']['inundtimeseed1'] # 8B5(i). Duration of inundation to cause stress for seedlings (number of timesteps)
inundtimesap1 = parameters['8. Recovery and mortality']['inundtimesap1'] # 8B5(ii). Duration of inundation to cause stress for saplings/adults (number of timesteps)
inundtimeseed2 = parameters['8. Recovery and mortality']['inundtimeseed2'] # 8B6(i). Duration of inundation to cause mortality for seedlings (number of timesteps)
inundtimesap2 = parameters['8. Recovery and mortality']['inundtimesap2'] # 8B6(ii). Duration of inundation to cause mortality for saplings/adults (number of timesteps)
recinundtimeseed = parameters['8. Recovery and mortality']['recinundtimeseed'] # 8B7(i). Duration of recovery from inundation stress for seedlings (number of timesteps)
recinundtimesap = parameters['8. Recovery and mortality']['recinundtimesap'] # 8B7(ii). Duration of recovery from inundation stress for saplings/adults (number of timesteps)
addinundtimeseed = parameters['8. Recovery and mortality']['addinundtimeseed'] # 8B8(i). Additional recovery time due to second event occurring prior to recovery for seedlings (number of timesteps)
addinundtimesap = parameters['8. Recovery and mortality']['addinundtimesap'] # 8B8(ii). Additional recovery time due to second event occurring prior to recovery for saplings/adults (number of timesteps)
maxinundevents = parameters['8. Recovery and mortality']['maxinundevents'] #8B9. Maximum number of inundation events limiting recovery timeframe (-)
# 8C. Desiccation
desiccactive = parameters['8. Recovery and mortality']['desiccactive'] # 8C1. Desiccation mortality activation (Y/N)
wdthresh = parameters['8. Recovery and mortality']['wdthresh'] # 8C2. Water depth threshold for desiccation (m)
desicctimeseed1 = parameters['8. Recovery and mortality']['desicctimeseed1'] # 8C3(i). Duration of desiccation to initiate stress for seedlings (number of timesteps)
desicctimesap1 = parameters['8. Recovery and mortality']['desicctimesap1'] # 8C3(ii). Duration of desiccation to initiate stress for saplings/adults (number of timesteps)
desicctimeseed2 = parameters['8. Recovery and mortality']['desicctimeseed2'] # 8C4(i). Duration of desiccation to initiate mortality for seedlings (number of timesteps)
desicctimesap2 = parameters['8. Recovery and mortality']['desicctimesap2'] # 8C4(ii). Duration of desiccation to initiate mortality for saplings/adults (number of timesteps)
recdesicctimeseed = parameters['8. Recovery and mortality']['recdesicctimeseed'] # 8C5(i). Duration of recovery from desiccation stress for seedlings (number of timesteps)
recdesicctimesap = parameters['8. Recovery and mortality']['recdesicctimesap'] # 8C5(ii). Duration of recovery from desiccation stress for saplings/adults (number of timesteps)
adddesicctimeseed = parameters['8. Recovery and mortality']['adddesicctimeseed'] # 8C6(i). Additional recovery time due to second event occurring prior to recovery for seedlings (number of timesteps)
adddesicctimesap = parameters['8. Recovery and mortality']['adddesicctimesap'] # 8C6(ii). Additional recovery time due to second event occurring prior to recovery for saplings/adults (number of timesteps)
maxdesiccevents = parameters['8. Recovery and mortality']['maxdesiccevents'] # 8C7. Maximum number of desiccation events limiting recovery timeframe (-)
# 8D. Burial
burialactive = parameters['8. Recovery and mortality']['burialactive'] # 8D1. Burial mortality activation (Y/N)
burialeventtime = parameters['8. Recovery and mortality']['burialeventtime'] # 8D2. Maximum number of timesteps over which sedimentation is assumed to be sudden (number of timesteps)
burialseedht1 = parameters['8. Recovery and mortality']['burialseedht1'] # 8D3. Proportion of seedling height buried to induce stress (-)
burialseedht2 = parameters['8. Recovery and mortality']['burialseedht2'] # 8D4. Proportion of seedling height buried to cause mortality (-)
burialpneht1 = parameters['8. Recovery and mortality']['burialpneht1'] # 8D5. Proportion of pneumatophore height buried to induce stress (-)
burialpneht2 = parameters['8. Recovery and mortality']['burialpneht2'] # 8D6. Proportion of pneumatophore height buried to cause mortality (-)
burialmorttime = parameters['8. Recovery and mortality']['burialmorttime'] # 8D7. Duration of burial required to cause mortality (number of timesteps)
# 8E. Senescence
senesactive = parameters['8. Recovery and mortality']['senesactive'] # 8E1. Senescence mortality activation (Y/N)
maxage = parameters['8. Recovery and mortality']['maxage'] # 8E2. Maximum age (years)
# 8F. System removal
sysremactive = parameters['8. Recovery and mortality']['sysremactive'] # 8F1. System removal activation (Y/N)
sysremthresh = parameters['8. Recovery and mortality']['sysremthresh'] # 8F2. Stem height threshold for removal of dead mangroves from the system (m)
# 9. Functionality
funcactive = parameters['9. Functionality']['funcactive'] # 9A. Functionality stage activation (Y/N)
projareastem = parameters['9. Functionality']['projareastem'] # 9B. Projected area (stems) (m2 or formula)
projareapneum = parameters['9. Functionality']['projareapneum'] # 9C. Projected area (pneumatophores) (m2 or formula)
projvolstem = parameters['9. Functionality']['projvolstem'] # 9D. Projected volume (stems) (m3 or formula)
projvolpneum = parameters['9. Functionality']['projvolpneum'] # 9E. Projected volume (pneumatophores) (m3 or formula)
dragcoeff = parameters['9. Functionality']['dragcoeff'] # 9F. Drag coefficient (number or formula)
percentileveg = parameters['9. Functionality']['percentileveg'] # 9G. Percentile of stem heights and diameters in each grid cell that are returned to DFM (-)
biomassactive = parameters['9. Functionality']['biomassactive'] # 9H. Biomass activation (Y/N)
agb = parameters['9. Functionality']['agb'] # 9I. Above-ground biomass (kg or formula)
bgb = parameters['9. Functionality']['bgb'] # 9J. Below-ground biomass (kg or formula)
autochactive = parameters['9. Functionality']['autochactive'] # 9K. Autochthonous accretion activation (Y/N)

# ======================================================================
# READ WATER LEVEL DATA
# ======================================================================
# Function to read water level data for manual manipulation to retrieve statistics that are unavailable in DFM
def read_water_level_data(file_path):
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            if line.strip() and not line.startswith('[') and not '=' in line:
                time, water_level = map(float, line.split())
                data.append((time, water_level))
        return np.array(data)

# Read water level time series data
print('Reading water level time series data...')
wldata = read_water_level_data(wl_file) # This array comprises the number of seconds since the start of the model run, and the water level

# Function for calculating the minimum water level for each timestep interval
def calculate_min_water_levels(data, mtpervt, nt):  
    min_levels = []
    for i in range(nt):
        start_time = i * mtpervt
        end_time = start_time + mtpervt
        interval_data = data[(data[:, 0] >= start_time) & (data[:, 0] < end_time)]
        if interval_data.size > 0:
            min_level = np.min(interval_data[:, 1])
            min_levels.append((i, min_level))
        else:
            min_levels.append((i, np.nan))
    return np.array(min_levels)

# Calculate minimum water level data
minwl = calculate_min_water_levels(wldata, mtpervt, nt)

# ======================================================================
# DFM COORDINATES AND GRID CELLS
# ======================================================================
# Retrieving x and y coordinates from DFM
dataset = nc.Dataset(map_file)
x_coord = dataset.variables['mesh2d_face_x'][:] # The length of this dataset corresponds to the number of grid cells
y_coord = dataset.variables['mesh2d_face_y'][:] # The length of this dataset corresponds to the number of grid cells
x_coord2 = dataset.variables['mesh2d_node_x'][:] # The length of this dataset corresponds to the number of grid cells + boundary cells
y_coord2 = dataset.variables['mesh2d_node_y'][:] # The length of this dataset corresponds to the number of grid cells + boundary cells
x_coord3 = dataset.variables['mesh2d_edge_x'][:] # The length of this dataset corresponds to the number of edges
y_coord3 = dataset.variables['mesh2d_edge_y'][:] # The length of this dataset corresponds to the number of edges
# Output coordinates to csv files
pandas.DataFrame(x_coord).to_csv(os.path.join(csvfolder_path,'x_coord.csv'),index=False)
pandas.DataFrame(x_coord2).to_csv(os.path.join(csvfolder_path,'x_coord2.csv'),index=False)
pandas.DataFrame(x_coord3).to_csv(os.path.join(csvfolder_path,'x_coord3.csv'),index=False)
pandas.DataFrame(y_coord).to_csv(os.path.join(csvfolder_path,'y_coord.csv'),index=False)
pandas.DataFrame(y_coord2).to_csv(os.path.join(csvfolder_path,'y_coord2.csv'),index=False)
pandas.DataFrame(y_coord3).to_csv(os.path.join(csvfolder_path,'y_coord3.csv'),index=False)
# Define number of grid cells
numcells = len(y_coord) # This is equal to the number of entries in the y-coordinate variable
print("Number of grid cells: ", numcells)
# Define max. stem carrying capacity of each grid cell
gridcellsize = ba[0] # This assumes that every grid cell in the model is the same size. Future revisions of the model will incorporate varying grid cell sizes
print('Grid cell size:', gridcellsize)
maxstemsgrid = gridcellsize * maxstemsmsq
maxcapacity = int(maxstemsgrid)
print("Max. capacity of grid cells:", maxstemsgrid, "stems")

# ======================================================================
# PRINT MODEL TIMES
# ======================================================================
# Obtain wet and dry cell threshold from the model
epshu = float(getmdu[(getmdu['section'] == 'numerics') & (getmdu['key'] == 'Epshu')]['value'].values[0])
# Parse start and end dates from DFM's time, string format, into the datetime format in Python
refdatet = parser.parse(refdate.iloc[0])
tstartt = timedelta(seconds=float(tstart.iloc[0]))
tendd = timedelta(seconds=float(tend.iloc[0]))
# Reference times
reftime = refdatet + tstartt
timeend = refdatet + tendd
timendwmorf = refdatet + (tendd * morfac)
print('Start simulation time:', reftime)
print('End of simulation in DFM:', timeend)
print('End of simulation with MorFac:', timendwmorf)
# Timesteps in days, months, and years
model_days = int(tend) / (nt*86400) # This is to define the number of days for the time between mangrove model updates (i.e., per timestep)
model_months = model_days / 30 # This is to define the number of months per timestep
model_years = model_days / 365 # This is to define the number of years per timestep
print('Model spans ', model_days*nt, ' days, or ', model_months*nt, ' months')
print('Number of mangrove model timesteps: ', nt)
print('Timestep length: ', mtpervt, ' seconds, or ', model_days, ' days, or ', model_months, ' months, or ', model_years, ' years')
# Initial restoration / seeding window
print("Initial seeding window start: ",seedwinstart.strftime("%x"))
print("Initial seeding window end: ",seedwinend.strftime("%x"))

# ======================================================================
# MANGROVE ATTRIBUTE MATRIX CREATION
# ======================================================================
# The following matrices are developed to monitor the physical characteristics, age, and number of stems in each grid cell. These are preallocated with zero
# Note that the DFM variables refer to "faces" and provide data at the centre of the grid cell (with additional entries for boundaries)

# Mangrove attributes
# General
vegtrk = np.zeros((maxcapacity,ndx)) # Tracker matrix for individual mangroves (including anchored seeds)
# A 1 is included in the matrix vegtrk if a mature stem exists, and a 0 if it doesn't.
stemtrk = np.zeros((maxcapacity,ndx)) # Tracker purely for stems to inform stem density calcs (excluding anchored seeds)
vegcheck = np.zeros(ndx) # Create matrix for the number of grid cells to check if mangroves exist
vegtrkadult = np.zeros(len(x_coord)) # Tracker for whether each grid cell containts 'adult' mangroves (used to determine if neighbouring cells contain fecund mangroves)
# Note that this array does not include boundary cells to maintain consistency with x and y co-ordinates.
vegtrkage = np.zeros((maxcapacity,ndx)) # As above but for age
vegtrktime = np.zeros((maxcapacity,ndx)) # To track the timestep that vegetation established
# Stems
vegtrkdia = np.zeros((maxcapacity,ndx)) # Tracker for stem diameter
vegtrkht = np.zeros((maxcapacity,ndx)) # Tracker for stem height
# Cable roots
vegtrkrt = np.zeros((maxcapacity,ndx)) # Tracker for cable root length
vegtrkrtdia = np.zeros((maxcapacity,ndx)) # Tracker for cable root diameter
# Pneumatophores
maxpneums = int(maxrtlen / pneumint) # Maximum number of pneumatophores per cable root for each stem (based on max. cable root length and spacing between pneumatophores)
vegtrkpnechk = np.zeros((maxcapacity,ndx,maxpneums)) # Tracker for pneumatophores presence
# A 1 is included in the matrix vegtrkpnechk if a pneumatophore exists, and a 0 if it doesn't.
vegtrkpnelen = np.zeros((maxcapacity,ndx,maxpneums)) # As above but for pneumatophore length
vegtrkpnedia = np.zeros((maxcapacity,ndx,maxpneums)) # Tracker for pneumatophore base diameter
agpnelen = np.zeros((maxcapacity,ndx,maxpneums)) # For aboveground pneumatophore length
# Bed level
vegtrkaut = np.zeros((maxcapacity,ndx)) # As above but for autochthonous accretion
vegtrkbed = np.zeros((maxcapacity,ndx)) # As above but for tracking the bed level that the vegetation was established on
autochgrd = np.zeros(ndx) # Array to represent the total autochthonous accretion in each grid cell. This is just 1 x number of grid cells
autochgrdavg = np.zeros(ndx) # Array to represent the average autochthonous accretion in each grid cell (to be passed back to bed level)
blh = np.zeros(ndx) # Setting up the matrix for saving historical bed levels (used when comparing current bed level to previous) for each grid cell
# Create lists to set all values to 0 when vegetation dies
vegtrklist = [vegtrk,stemtrk,vegtrkrt,vegtrkrtdia,vegtrkdia,vegtrkage,vegtrkaut,vegtrkbed,vegtrktime] # List of all arrays with 2 dimensions that should be reset when mangroves die

# Mortality and Recovery
# Damage regimes for individual mangroves - Set with 0s for healthy state. These values will change as per the relevant damage regime
dmgregin = np.zeros((maxcapacity,ndx)) # Damage regime tracker for inundation
dmgregde = np.zeros((maxcapacity,ndx)) # Damage regime tracker for desiccation
dmgregbu = np.zeros((maxcapacity,ndx)) # Damage regime tracker for burial
dmgregse = np.zeros((maxcapacity,ndx)) # Damage regime tracker for senescence
dmgregall = np.zeros((maxcapacity,ndx)) # Damage regime tracker across all events. This takes the value of the highest number from the above regime trackers
dmgreglist = {} # Initialise damage regime dictionary to store lists for all stem, grid combinations
# Damage regimes for each grid cell - The value for each grid cell is set to the maximum value (worst case) of all mangroves in the cell for plotting purposes
dmgregingc = np.zeros(ndx) # Damage regime tracker for inundation in each grid cell
dmgregdegc = np.zeros(ndx) # Damage regime tracker for desiccation in each grid cell
dmgregbugc = np.zeros(ndx) # Damage regime tracker for burial in each grid cell
dmgregsegc = np.zeros(ndx) # Damage regime tracker for senescence in each grid cell
dmgregallgc = np.zeros(ndx) # Damage regime tracker for all events in each grid cell. This is the average (or percentile) health of all stems in the grid cell based on the overall tracker above
# Inundation
# Seedling
recinseed = np.zeros((maxcapacity,ndx)) # Tracker for recovery time for individual mangroves
recinseed = (recinseed + recinundtimeseed).astype(int) # Assign input value for recovery
# Sapling/Adult
recinsap = np.zeros((maxcapacity,ndx)) # Tracker for recovery time for individual mangroves
recinsap = (recinsap + recinundtimesap).astype(int) # Assign input value for recovery
# Desiccation
# Seedling
recdeseed = np.zeros((maxcapacity,ndx)) # Tracker for recovery time for individual mangroves
recdeseed = (recdeseed + recdesicctimeseed).astype(int) # Assign input value for recovery
# Sapling/Adult
recdesap = np.zeros((maxcapacity,ndx)) # Tracker for recovery time for individual mangroves
recdesap = (recdesap + recdesicctimesap).astype(int) # Assign input value for recovery

# Historical variable trackers per grid cell
bldiff = np.zeros(ndx) # Tracker to store the difference in bed levels from one timestep to the next
bldiffcum = np.zeros(ndx) # Tracker to store the cumulative bed level difference during burial
bltab = np.zeros(ndx) # Tracker to store the timestep at which accretion starts for the grid cell. This is reset to 0 if accretion stops
blgc = np.zeros(ndx) # Tracker to store the bed level per grid cell for plotting purposes
maxwd = np.zeros(ndx) # Tracker to store maximum water depth
minwd = np.zeros(ndx) # Tracker to store minimum water depth
maxbss = np.zeros(ndx) # Tracker to store maximum bed shear stress

# Functionality parameters per grid cell
vmpneum = np.zeros(ndx) # Tracker for storing pneumatophore volume per grid cell
ampneum = np.zeros(ndx) # Tracker for storing pneumatophore area per grid cell
vmstem = np.zeros(ndx) # Tracker for storing stem volume per grid cell
amstem = np.zeros(ndx) # Tracker for storing stem area per grid cell
vm = np.zeros(ndx) # Tracker for storing total submerged volume per grid cell
am = np.zeros(ndx) # Tracker for storing total submerged area per grid cell
vcontrol = np.zeros(ndx) # Tracker for storing control volume per grid cell
L = np.zeros(ndx) # Tracker for storing characteristic length per grid cell

# Attribute trackers for plotting
n_pne = np.zeros(ndx) # Tracker for average pneumatophore density
hv_stem = np.zeros(ndx) # Tracker for average stem height per grid cell
hv_pneum = np.zeros(ndx) # Tracker for average pneum height per grid cell
stemdia = np.zeros(ndx) # Tracker for average stem diameter per grid cell
pneumdia = np.zeros(ndx) # Tracker for average pneum diameter per grid cell
vegbed = np.zeros(ndx) # Tracker for average initial mangrove bed level per grid cell

# Carbon Storage and Biomass
agb_stem = np.zeros((maxcapacity,ndx)) # AGB for stems
agb_pne = np.zeros((maxcapacity,ndx,maxpneums)) # AGB for pneumatophores
agb_pne_sum = np.zeros((maxcapacity,ndx)) # AGB sum for pneumatophores for each stem
agb_total = np.zeros((ndx)) # AGB total per grid cell
bgb_stem = np.zeros((maxcapacity,ndx)) # BGB for each stem
bgb_total = np.zeros(ndx) # BGB total per grid cell

# ======================================================================
# INITIAL MANGROVE SETUP
# ======================================================================
# This is the setup for mangroves that already exist at the start of the model run
# Set trackers to match initial planting / existing mangrove attributes
for j in range(ndx): # For every grid cell
    if stemheight[j] > 0: # If a stem height exists
        if iniage * 365 * 24 * 3600 / mtpervt >= fecundity * 365 * 24 * 3600 / mtpervt: # If the initial age is >= the age for fruiting, set the mature mangrove tracker to 1 for the grid cell
            vegtrkadult[j] = 1
        inistemsgrid = int(iniden * ba[j]) # Set the initial number of stems per grid cell based on the initial density
        for b in range(inistemsgrid): # For every initial stem
                vegtrk[b,j] = 1 # Set the mangrove tracker to be 1 at this position (i.e., a stem / established seed exists here)
                stemtrk[b,j] = 1 # Set the stem tracker to be 1 at this position (i.e., there is a stem here)
                vegtrkbed[b,j] = bedlevel[j] # Set bed level of each initial mangrove to match the grid cell bed level
                vegtrkht[b,j] = stemheight[j] # Set mangrove height
                vegtrkrt[b,j] = inirtlen # Set the cable root length
                vegtrkrtdia[b,j] = inirtdia # Set the cable root diameter
                vegtrkdia[b,j] = inistemdia # Set the stem diameter
                vegtrkage[b,j] = iniage * 365 * 24 * 3600 / mtpervt # Set the age in terms of no. of time steps
                if floor(vegtrkrt[b,j]/pneumint) >= 1: # If the cable root reaches next interval and tracker matrix has no pneumatophore
                    for r in range(floor(vegtrkrt[b,j]/pneumint)): # For all intervals up to the length of the cable root that pneumatophores could exist
                        vegtrkpnechk[b,j,r] = 1 # Assign a value of 1 to the tracker matrix to highlight the presence of a pneumatophore (at each interval)
                        vegtrkpnelen[b,j,r] = inipnelen # Give an initial length to the pneumatophore at timestep 0
                        vegtrkpnedia[b,j,r] = inipnedia # Give an initial base diameter to the pneumatophore

# Function to check whether a non-zero value exists in the mangrove presence tracker (i.e., mangroves are present)
def check_column(matrix, column_index):
    column = matrix[:, column_index]  # Extract the specified column
    nonzero_elements = np.nonzero(column)  # Find the indices of non-zero elements
    return bool(nonzero_elements[0].size)  # Return True if non-zero elements exist, False otherwise
# Populate tracking array with 1 where initial mangroves exist
for p in range(0,ndx-1):
    vegcheck[p] = check_column(vegtrk,p)

# ======================================================================
# ESTABLISHMENT SETUP
# ======================================================================
# Function to check if the current date in the model falls within the fruiting window each year
def date_within_range(current_date, start_date, end_date):
    if seedwinrep == 'Y':
        # Adjust the start and end dates to have the same year (in case of the event where the start date falls at the end of one year, and the end date at the start of the next)
        adjusted_start_date = datetime(current_date.year,start_date.month,start_date.day)
        adjusted_end_date = datetime(current_date.year,end_date.month,end_date.day)
        # Output if the current date falls within the range from the initial seeding window
        if adjusted_start_date <= adjusted_end_date: # When the initial seeding window dates are in the same year
            return adjusted_start_date <= current_date <= adjusted_end_date
        else: # When the initial seeding window dates span two calendar years
            return current_date >= adjusted_start_date or current_date <= adjusted_end_date
    else:
        return start_date <= current_date <= end_date

# Function for evaluating neighbouring cells (assuming mangroves only establish in cells adjacent to, or in the same cell as, those with existing fecund mangroves)
# Note that this assumes square cells
vegexist = np.column_stack((x_coord,y_coord,vegtrkadult))
def neighbourveg(j, threshold):
    cell_x, cell_y, cell_veg = vegexist[j]
    for neighbour_data in vegexist:
        neighbour_x, neighbour_y, neighbour_veg = neighbour_data
        distance = np.sqrt((cell_x - neighbour_x)**2 + (cell_y - neighbour_y)**2)
        if neighbour_veg and distance < threshold:
            return True
    return False

# WoO 1 - Inundation
# Set inundation free period
inundfreestep = int(inundfreepd * 3600 / mtpervt) # Required inundation free period in terms of the number of timesteps 
print("Inundation free period: ",inundfreepd,"hrs (",inundfreestep," timesteps)")

# ======================================================================
# PLOTTING SETUP
# ======================================================================
# Function for calculating a user-defined percentile of non-zero values for key parameters to include in plotting
# This simplifies the output per grid cell
def calculate_percentile(original_variable, percentile_value, output_variable):
    # Find non-zero indices
    nz_indices = np.nonzero(original_variable)
    if nz_indices[0].size > 0:
        # Calculate percentile for each unique value along the second dimension
        unique_indices = np.unique(nz_indices[1]) # Note that the second dimension refers to each grid cell. This should be changed if the plotting is to consider other dimensions
        for idx in unique_indices:
            # Find non-zero values corresponding to current unique index
            idx_mask = nz_indices[1] == idx
            non_zero_values = original_variable[nz_indices[0][idx_mask], nz_indices[1][idx_mask]]
            non_zero_values = non_zero_values[non_zero_values != 0] # Filter out zero values
            # Calculate percentile for non-zero values
            percentile_value = np.percentile(non_zero_values, percentile_value)
            # Update corresponding entry in output variable
            output_variable[idx] = percentile_value
    return output_variable

# Colour parameters for cross-section plotting
# Mortality status colour coding
colour_healthy = '#7b9a6e'
colour_recovery = '#e9f59a'
colour_stress = '#fad182'
colour_mortality = '#fc9488'
# Function to map mortality values to colours
def map_mortality_to_colour(value):
    if value == 0:
        return colour_healthy
    elif 0 < value <= 1:
        return colour_recovery
    elif 1 < value < 2:
        return colour_stress
    elif value == 2:
        return colour_mortality
    return colour_healthy # Default colour for values that are outside limits

# ======================================================================
# HISTORICAL VARIABLE TRACKERS
# ======================================================================
# These matrices are created to track the duration for which logical expressions are true
# These counters increase by 1 for every timestep that the associated conditions are met
# When the conditions are no longer met, the counter resets to 0
# WoO1
c_woo1 = np.zeros((maxcapacity, ndx)) # Counter for inundation free steps
# Inundation
c_instr = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is stressed (multiple of 1). This reverts to 0 if healthy.
c_inrec = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is in recovery (multiple of 1). This reverts to 0 if healthy.
# Desiccation
c_destr = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is stressed (multiple of 1). This reverts to 0 if healthy.
c_derec = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is in recovery (multiple of 1). This reverts to 0 if healthy.
# Burial
c_bustr = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is stressed (multiple of 1). This reverts to 0 if healthy.
c_burec = np.zeros((maxcapacity,ndx)) # Counter for number of timesteps the stem is in recovery (multiple of 1). This reverts to 0 if healthy.

########################################################################
# START MODEL RUN
########################################################################
# ======================================================================
# TIMESTEP LOOP
# ======================================================================
# Update and run model for each mangrove timestep
for i in range(0,nt): # Where i represents the time step
    
    # ======================================================================
    # INITIAL SETUP
    # ======================================================================
    # Timesteps
    model_dimr.update(mtpervt) # Update the model for the timeframe indicated by mtpervt (i.e., run the hydrodynamics in DFM for this length of time)
    print("Timestep: ", i) # This print message is to track the current timestep as the code runs
    currenttime = reftime + timedelta(seconds=mtpervt * i)
    print('Current time:',currenttime)
    
    # BMI parameter exchange during model run
    bedlevel = model_dfm.get_var('bl') # Bed level is retrieved here again because this should change every timestep
    bedlevel = np.array(bedlevel)
    is_maxvalsnd = model_dfm.get_var('is_maxvalsnd') # Retrieve maximum variables from DFM (this creates a 3d array with 1) max bed shear stress, 2) max velocity, 3) max water depth)
    maxwd = is_maxvalsnd[range(ndx),2] # Create an array with max water depth values (for the timestep) for each non-boundary grid cell   
    maxbss = is_maxvalsnd[range(ndx),0] # Create an array with max bed shear stress (for the timestep) for each non-boundary grid cell

    # Calculate minimum water depth using calculated minimum water levels and the bed level from DFM
    # Min water depth for every grid cell
    minwd = minwl[i,1] - bedlevel # The minimum water depth (with rows for every grid cell) is calculated as the minimum water level for the timestep minus the bed level for the timestep
    minwd = np.maximum(minwd,0)
    
    # Store historical variables
    bl = model_dfm.get_var('bl') # Bed level variable retrieved again but not set as a numpy array
    if i == 0:
        blini = bl # Set initial bed level at the start of the model run
    blh = bl # Store the bed level for all grid cells at this timestep in a numpy array blh
    if i > 0:
        bldiff = blh - blprev # Single column with difference in bed levels from the previous time step
    # Set previous bed level for the next timestep
    blprev = bl

    # ======================================================================
    # GRID CELL LOOP
    # ======================================================================
    # This is a for-loop for every grid cell in the model, which then runs through various conditional statements over the mangrove lifecycle     
    for j in range(ndx): # Where j represents the grid cell
              
        # Update age of existing vegetation before the rest of the calculations take place
        if i > 0:
            stem_exists = vegtrk[:,j] != 0
            vegtrkage[stem_exists,j] += 1

        # Calculate grid cell area
        cellarea = ba[j]

        # Set the vegetation tracker for each grid cell to be equal to 1 if there is a stem within the grid cell that is of an appropriate age for fruiting
        if j < len(x_coord): # Limiting the calculation to non-boundary cells
            if any(vegtrkage[b,j] >= fecundity * 365 * 24 * 3600 / mtpervt for b in range(maxcapacity)): # If the mangrove age in timesteps is greater than fecundity age in timesteps
                vegtrkadult[j] = 1 # Set fecund mangrove (adult) tracker to 1 for the grid cell

        # ======================================================================
        # ESTABLISHMENT
        # ======================================================================
        if estactive == 'Y': # Check if user has activated the establishment stage
            if j < len(x_coord) and iniveg == 'Y' and neighbourveg(j, seedtrvldist) or iniveg == 'N': # If 'Y', then neighbouring vegetation rule applies. If 'N', then seeds are assumed to be available
                for b in range(maxcapacity): # For every stem in the grid cell
                    if vegtrkht[b,j] < sapht: # This ensures the stem is still in the seedling stage
                        
                        ## Windows of Opportunity ##
                        # WoO 1 - Inundation #
                        if vegtrk[b,j] == 0 and date_within_range(currenttime,seedwinstart,seedwinend): # If veg doesn't exist and the time lies within the fruiting window
                            if maxwd[j] < epshu: # If the max water depth for the grid cell is lower than the dry cell threshold
                                c_woo1[b,j] += 1 # Add 1 to the WoO1 counter (i.e., the grid cell is inundation free for another timestep)
                            if c_woo1[b,j] >= inundfreestep: # If the counter reaches a value greater than or equal to the required inundation free period for establishment
                                # Inundation free period is sufficient for propagules to anchor
                                c_woo1[b,j] = 0 # Reset counter to enable future re-growth should the newly anchored seed fail in the future
                                # Random chance of establishment                                                                                                                                                                                                                                                                       
                                estrandom = np.random.uniform(0,1,1) # Generates random value for comparison at each stem iteration (instead of creating a matrix of random values)
                                if estrandom < estchance: # Include restriction on random chance of establishment (by including this within the for-loop, more than one stem can establish per timestep)
                                    # Storing mangrove variables                   
                                    vegtrk[b,j] = 1 # Set mangrove tracker to include a 1 when the propagule has anchored
                                    vegtrktime[b,j] = i # Set the mangrove time tracker to be the timestep that the mangrove established. Note that this will revert to 0 if washed away
                                    vegtrkbed[b,j] = blh[j] # Assign the bed level at the current time step for the grid cell, to the tracker matrix to store the bed level of establishment
                            if maxwd[j] > epshu: # If the max water depth for the grid cell is above the dry cell threshold
                                c_woo1[b,j] = 0 # Reset counter
                        
                        # WoO 2 - Hydrodynamics #
                        # Check that seeds are not washed away between beginning root and stem growth
                        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for the critical bed shear stress calculations
                        if taucrit == "balke2015": # A. alba propagules in mangrove mud (Balke et al., 2015)
                            if vegtrkrt[b,j] >= 0.058/(0.4135*100): # Check root length is positive and within validity of the equation
                                taudisl = 0.4135 * vegtrkrt[b,j]*100 - 0.058 # Bed shear stress at dislodgement. Root length in m.
                        elif isinstance(taucrit, (int, float)): # If the user has defined a single numerical value instead of an equation
                            taudisl = taucrit
                        if vegtrk[b,j] == 1: # If the propagule has rooted
                            if 'taudisl' in locals():
                                if maxbss[j] >= taudisl: # If the bed shear stress is greater than the dislodgement stress
                                    for item in vegtrklist: # Set all mangrove attributes to zero (including the tracker)
                                        item[b,j] = 0
                                    vegtrkht[b,j] = 0
                                    vegtrkpnechk[b,j,:] = 0
                                    vegtrkpnelen[b,j,:] = 0
                                    vegtrkpnedia[b,j,:] = 0

                        ### WoO 3 - Erosion ###
                        # If there is erosion and the absolute difference in bed level is greater than the critical erosion, set veg tracker to 0.
                        if stemtrk[b,j] == 1: # If the seedling has a stem
                            # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for the critical erosion depth calculations
                            if erocrit == "balke2015": # Avicennia marina, Balke et al., 2015
                                if vegtrkht[b,j] >= math.exp(5.0584/3.5058) / 100: # Check stem height is positive and within validity of the critial erosion equation
                                    ecrit = 3.5058 * math.log(vegtrkht[b,j]*100) - 5.0584 # Critical erosion limit in centimetres
                                    ecrit = ecrit*100 # Critical erosion limit in metres (positive values indicate positive depth of erosion)
                            elif isinstance(erocrit, (int, float)): # If the user has defined a single numerical value instead of an equation
                                ecrit = erocrit # metres
                            if 'ecrit' in locals():
                                if ecrit > 0 and bldiff[j] < 0 and np.abs(bldiff[j]) > ecrit: # If there is erosion and the absolute difference in bed level is greater than the critical erosion, set veg tracker to 0.
                                    # Note that the burial here is only calculated using the difference in bed level between two timesteps, rather than over a long period
                                    # print('ero check')
                                    for item in vegtrklist:
                                        item[b,j] = 0 # Set all veg properties to zero (including the veg tracker if the bed shear stress exceeds critical erosion or burial thresholds)
                                    vegtrkht[b,j] = 0
                                    vegtrkpnechk[b,j,:] = 0
                                    vegtrkpnelen[b,j,:] = 0
                                    vegtrkpnedia[b,j,:] = 0

        # ======================================================================
        # RECOVERY AND MORTALITY
        # ======================================================================
        # Note that this stage is calculated prior to mangrove growth
        # Mortality events in this version of the model: 1. Inundation, 2. Desiccation, 3. Burial, 4. Senescence
        # Damage regimes for each mortality event: A. Stress, B. Mortality, C. System Removal
        # These damage regimes are defined by threshold limits, where one threshold defines stress, and a higher / more intense threshold defines mortality
        # For each damage regime, each stem require an approximate timeframe in which it can recover
        # For damage regime B. Mortality, the recovery time is 0 as the mangrove is unable to recover
        # This mangrove still remains in the system (unless system removal is activated) because uprooting or overturning is not included at this stage of model development 
        if mortactive == 'Y': # Check if user has activated the recovery and mortality stage
            for b in range(0,maxcapacity-1): # For every possible stem
                if vegtrkht[b,j] != 0: # For every non-zero value in the stem height array

                    ## Seedling phase ##
                    if vegtrkht[b,j] < sapht: # If the mangrove stem is a seedling (based on height)
                        # 1. Inundation #
                        if inundactive == 'Y': # Check if user has activated inundation
                            if dmgregin[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                if minwd[j] >= max(inundstemht * vegtrkht[b,j], inundhtthresh): # If the min water depth over the timestep is greater than threshold
                                    c_instr[b,j] += 1 # Update counter for being impacted by inundation
                                if minwd[j] < max(inundstemht * vegtrkht[b,j], inundhtthresh): # If the min water depth over the timestep is lower than threshold
                                    c_instr[b,j] = 0 # Reset counter if not impacted by inundation
                                if c_instr[b,j] == inundtimeseed1: # If the counter reaches the threshold for the inundation to induce stress
                                    dmgregin[b,j] = 1 # Set damage regime to A. Stress (1)
                                if c_instr[b,j] == inundtimeseed2: # If the counter reaches the threshold for the inundation to cause mortality
                                    dmgregin[b,j] = 2 # Set damage regime to B. Mortality (2)
                                if dmgregin[b,j] == 1 and not minwd[j] >= max(inundstemht * vegtrkht[b,j], inundhtthresh): # If in stress but now no longer impacted by inundation, set to recovery
                                    dmgregin[b,j] = 0.5 # Initialise the recovery phase
                                    c_instr[b,j] = 0 # Reset counter for being impacted by inundation
                                # If in recovery and impacted by inundation (but not for long enough to induce stress), then the recovery timeframe increases
                                if dmgregin[b,j] == 0.5 and minwd[j] >= max(inundstemht * vegtrkht[b,j], inundhtthresh):
                                    if recinseed[b,j] < recinundtimeseed + addinundtimeseed * maxinundevents: # Recovery timeframe limited to the user defined number of events
                                        recinseed[b,j] += addinundtimeseed # Replace the recovery period with a higher period due to the occurrence of an extreme event
                                # If in recovery and not impacted by inundation and the recovery timeframe has not ended, then remain in recovery
                                if dmgregin[b,j] == 0.5 and not minwd[j] >= max(inundstemht * vegtrkht[b,j], inundhtthresh) and c_inrec[b,j] < recinseed[b,j]:
                                    c_inrec[b,j] += 1 # Update recovery counter
                                # If in recovery and the recovery timeframe has ended, set to healthy
                                if dmgregin[b,j] == 0.5 and c_inrec[b,j] >= recinseed[b,j]:
                                    dmgregin[b,j] = 0 # Return damage regime to a healthy state (0)
                                    c_inrec[b,j] = 0 # Reset recovery counter

                        # 2. Desiccation #
                        if desiccactive == 'Y': # Check if user has activated desiccation
                            if dmgregde[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                if maxwd[j] < wdthresh: # If max water depth is lower than the water depth threshold for desiccation
                                    c_destr[b,j] += 1 # Update counter for being impacted by desiccation
                                if maxwd[j] > wdthresh: # If max water depth is greater than the water depth threshold for desiccation
                                    c_destr[b,j] = 0 # Reset counter if not impacted by desiccation
                                if c_destr[b,j] == desicctimeseed1: # If the counter reaches the threshold for the desiccation to induce stress
                                    dmgregde[b,j] = 1 # Set damage regime to A. Stress (1)
                                if c_destr[b,j] == desicctimeseed2: # If the counter reaches the threshold for the desiccation to cause mortality
                                    dmgregde[b,j] = 2 # Set damage regime to B. Mortality (2)
                                if dmgregde[b,j] == 1 and not maxwd[j] < wdthresh: # If in stress but now no longer being impacted by desiccation, set to recovery
                                    dmgregde[b,j] = 0.5 # Initialise the recovery phase
                                    c_destr[b,j] = 0 # Reset counter for being impacted by desiccation
                                # If in recovery and impacted by desiccation (but not for long enough to induce stress), then the recovery timeframe increases
                                if dmgregde[b,j] == 0.5 and maxwd[j] < wdthresh:
                                    if recdeseed[b,j] < recdesicctimeseed + adddesicctimeseed * maxdesiccevents: # Recovery timeframe limited to the user defined number of events
                                        recdeseed[b,j] += adddesicctimeseed # Replace the recovery period with a higher period due to the occurrence of an extreme event
                                # If in recovery and not impacted by desiccation and the recovery timeframe has not ended, then remain in recovery
                                if dmgregde[b,j] == 0.5 and not maxwd[j] < wdthresh and c_derec[b,j] < recdeseed[b,j]:
                                    c_derec[b,j] += 1 # Update recovery counter
                                # If in recovery and the recovery timeframe has ended, set to healthy
                                if dmgregde[b,j] == 0.5 and c_derec[b,j] >= recdeseed[b,j]:
                                    dmgregde[b,j] = 0 # Return damage regime to a healthy state (0)
                                    c_derec[b,j] = 0 # Reset recovery counter

                        # 3. Burial #
                        # Damage / mortality occurs when the current bed level is higher than the threshold proportion of the stem height + the bed level at establishment, for a set no. of timesteps
                        # Assumed that there is no recovery under severe burial events but recruitment still possible (Paling et al., 2008)
                        if burialactive == 'Y': # Check if user has activated burial
                            if dmgregbu[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                if bldiff[j] > 0: # If there is accretion
                                    if c_bustr[b,j] == 0: # If this is the first occasion of burial impact
                                        bltab[j] = bl[j] # Set tracker for bed level at time of burial impact 
                                        bldiffcum[j] += bldiff[j] # Add the change in bed level to a cumulative counter
                                        c_bustr[b,j] += 1 # Update counter for being impacted by burial
                                    else: # Else, the burial is continuing
                                        bldiffcum[j] += bldiff[j] # Add the change in bed level to a cumulative counter
                                        c_bustr[b,j] += 1 # Update counter for being impacted by burial
                                # If there is no further accretion this timestep, reset the burial trackers
                                # This is simplified such that the mangrove is returned to a healthy status as soon as there is a timestep where the accretion stops
                                # It is assumed that stems will be able to reach the surface should it no longer be accreting
                                if bldiff[j] <= 0:
                                    bldiffcum[j] = 0
                                    c_bustr[b,j] = 0
                                    dmgregbu[b,j] = 0 # Set damage regime to healthy (0)
                                # If cumulative burial within the event timeframe is greater than the stress proportion of stem height that was above the ground when the burial impact began, set to stress
                                if c_bustr[b,j] <= burialeventtime and bltab[j] - vegtrkbed[b,j] > 0 and bldiffcum[j] > burialseedht1 * (vegtrkht[b,j] - (bltab[j] - vegtrkbed[b,j])):
                                    dmgregbu[b,j] = 1 # Set matrix damage regime to A. Stress (1)
                                # If cumulative burial within the event timeframe is greater than the mortality proportion of stem height that was above the ground when the burial impact began, set to mortality
                                # Or if the stem has been in burial stress for longer than the mortality timeframe, set to mortality
                                if c_bustr[b,j] <= burialeventtime and bltab[j] - vegtrkbed[b,j] > 0 and bldiffcum[j] > burialseedht2 * (vegtrkht[b,j] - (bltab[j] - blini[j])) or c_bustr[b,j] > burialmorttime:
                                    dmgregbu[b,j] = 2 # Set matrix damage regime to B. Mortality (2)

                    ## Sapling / Adult phase ##
                    elif vegtrkht[b,j] >= sapht: # If stem is a sapling/adult
                        # 1. Inundation #
                        if inundactive == 'Y': # Check if user has activated inundation
                            if dmgregin[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                if minwd[j] >= max(inundhtthresh, inundpneumlen * np.average(agpnelen[b,j,:])): # If the min water depth over the timestep is greater than threshold
                                    c_instr[b,j] += 1 # Update counter for being impacted by inundation
                                if minwd[j] < max(inundhtthresh, inundpneumlen * np.average(agpnelen[b,j,:])): # If the min water depth over the timestep is lower than threshold
                                    c_instr[b,j] = 0 # Reset counter if not impacted by inundation
                                if c_instr[b,j] == inundtimesap1: # If the counter reaches the threshold for the inundation to induce stress
                                    dmgregin[b,j] = 1 # Set damage regime to A. Stress (1)
                                if c_instr[b,j] == inundtimesap2: # If the counter reaches the threshold for the inundation to cause mortality
                                    dmgregin[b,j] = 2 # Set damage regime to B. Mortality (2)
                                if dmgregin[b,j] == 1 and not minwd[j] >= max(inundhtthresh, inundpneumlen * np.average(agpnelen[b,j,:])): # If in stress but now no longer being impacted by inundation, set to recovery
                                    dmgregin[b,j] = 0.5 # Initialise the recovery stage
                                    c_instr[b,j] = 0 # Reset counter for being impacted by inundation
                                # If in recovery and impacted by inundation (but not for long enough to induce stress), then the recovery timeframe increases
                                if dmgregin[b,j] == 0.5 and minwd[j] >= max(inundhtthresh, inundpneumlen * np.average(agpnelen[b,j,:])):
                                    if recinsap[b,j] < recinundtimesap + addinundtimesap * maxinundevents: # Recovery timeframe limited to the user defined number of events
                                        recinsap[b,j] += addinundtimesap # Replace the recovery period with a higher period due to the occurrence of an extreme event
                                # If in recovery and not impacted by inundation and the recovery timeframe has not ended, then remain in recovery
                                if dmgregin[b,j] == 0.5 and not minwd[j] >= max(inundhtthresh, inundpneumlen * np.average(agpnelen[b,j,:])) and c_inrec[b,j] < recinsap[b,j]:
                                    c_inrec[b,j] += 1 # Update recovery counter
                                # If in recovery and the recovery timeframe has ended, set to healthy
                                if dmgregin[b,j] == 0.5 and c_inrec[b,j] >= recinsap[b,j]:
                                    dmgregin[b,j] = 0 # Return damage regime to a healthy state (0)
                                    c_inrec[b,j] = 0 # Reset recovery counter

                        # 2. Desiccation #
                        # If damage regime A is induced this timestep, set damage regime to A for this timestep
                        if desiccactive == 'Y': # Check if user has activated desiccation
                            if dmgregde[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                if maxwd[j] < wdthresh: # If max water depth is lower than the water depth threshold for desiccation
                                    c_destr[b,j] += 1 # Update counter for being impacted by desiccation
                                if maxwd[j] > wdthresh: # If max water depth is higher than the water depth threshold for desiccation
                                    c_destr[b,j] = 0 # Reset counter if not impacted by desiccation
                                if c_destr[b,j] == desicctimesap1: # If the counter reaches the threshold for the desiccation to induce stress
                                    dmgregde[b,j] = 1 # Set matrix damage regime to A. Stress (1)
                                if c_destr[b,j] == desicctimesap2: # If the counter reaches the threshold for the desiccation to cause mortality
                                    dmgregde[b,j] = 2 # Set matrix damage regime to B. Mortality (2)
                                if dmgregde[b,j] == 1 and not maxwd[j] < wdthresh: # If in stress but now no longer being impacted by desiccation, set to recovery
                                    dmgregde[b,j] = 0.5 # Initialise the recovery stage
                                    c_destr[b,j] = 0 # Reset counter for being impacted by desiccation
                                # If in recovery and impacted by desiccation (but not for long enough to induce stress), then the recovery timeframe increases
                                if dmgregde[b,j] == 0.5 and maxwd[j] < wdthresh:
                                    if recdesap[b,j] < recdesicctimesap + adddesicctimesap * maxdesiccevents: # Recovery timeframe limited to the user defined number of events
                                        recdesap[b,j] += adddesicctimesap # Replace the recovery period with a higher period due to the occurrence of an extreme event
                                # If in recovery and not impacted by desiccation and the recovery timeframe has not ended, then remain in recovery
                                if dmgregde[b,j] == 0.5 and not maxwd[j] < wdthresh and c_derec[b,j] < recdesap[b,j]:
                                    c_derec[b,j] += 1 # Update recovery counter
                                # If in recovery and the recovery timeframe has ended, set to healthy
                                if dmgregde[b,j] == 0.5 and c_derec[b,j] >= recdesap[b,j]:
                                    dmgregde[b,j] = 0 # Return damage regime to a healthy state (0)
                                    c_derec[b,j] = 0 # Reset recovery counter

                        # 3. Burial #
                        # Damage / mortality occurs when the current bed level is higher than the threshold proportion of the pneumatophores + the bed level at establishment, for a set no. of timesteps
                        # Assumed that there is no recovery under severe burial events but recruitment still possible (Paling et al., 2008)
                        if burialactive == 'Y': # Check if user has activated burial
                            if dmgregbu[b,j] != 2: # If the damage regime for the mangrove stem is not B (Mortality), then calculate updates to damage regime
                                nonzeropneums = agpnelen[b,j,:][agpnelen[b,j,:] != 0] # To create a parameter that only includes non-zero pneums. This is needed for using the np.average function below
                                if bldiff[j] > 0: # If there is accretion
                                    if c_bustr[b,j] == 0: # If this is the first occasion of burial impact
                                        bltab[j] = bl[j] # Set tracker for bed level at time of burial impact
                                        bldiffcum[j] += bldiff[j] # Add the change in bed level to a cumulative counter
                                        c_bustr[b,j] += 1 # Update counter for being impacted by burial
                                    else: # Else, the burial is continuing
                                        bldiffcum[j] += bldiff[j] # Add the change in bed level to a cumulative counter
                                        c_bustr[b,j] += 1 # Update counter for being impacted by burial
                                # If there is no further accretion this timestep, reset the burial trackers
                                # This is simplified such that the mangrove is returned to a healthy status as soon as there is a timestep where the accretion stops
                                # It is assumed that pneumatophores will be able to reach the surface should it no longer be accreting
                                if bldiff[j] <= 0:
                                    bldiffcum[j] = 0
                                    c_bustr[b,j] = 0
                                    dmgregbu[b,j] = 0 # Set damage regime to healthy (0)
                                # Calculate average pneumatophore height relative to recent bed level. This simplifies the computations instead of calculating this for every pneumatophore of each stem
                                if i < burialeventtime:
                                    recentpneumht = np.average(nonzeropneums)
                                else:
                                    recentpneumht = np.average(nonzeropneums) - (bltab[j] - blini[j])
                                # If cumulative burial within the event timeframe is greater than the stress proportion of stem height that was above the ground when the burial impact began, set to stress
                                if c_bustr[b,j] <= burialeventtime and bltab[j] - vegtrkbed[b,j] > 0 and bldiffcum[j] > burialpneht1 * recentpneumht:
                                    dmgregbu[b,j] = 1 # Set matrix damage regime to A. Stress (1)
                                # If cumulative burial over the event timeframe is greater than the mortality proportion of stem height that was above the ground when the burial impact began, set to mortality
                                # Or if the stem has been in burial stress for longer than the mortality timeframe, set to mortality
                                if c_bustr[b,j] <= burialeventtime and bltab[j] - vegtrkbed[b,j] > 0 and bldiffcum[j] > burialpneht2 * recentpneumht or c_bustr[b,j] > burialmorttime:
                                    dmgregbu[b,j] = 2 # Set matrix damage regime to B. Mortality (2)

                        # 4. Senescence #
                        # Mortality occurs when mangroves reach maximum age
                        if senesactive == 'Y': # Check if user has activated senescence
                            if vegtrkage[b,j] * model_years >= maxage: # If the mangrove age has reached the maximum age
                                dmgregse[b,j] = 2 # Set matrix damage regime to B. Mortality

                    # Set overall damage regime tracker
                    dmgreglist[(b,j)] = [dmgregin[b,j],dmgregde[b,j],dmgregbu[b,j],dmgregse[b,j]] # Create list of damage regimes
                    dmgregall[b,j] = max(dmgreglist[(b,j)]) # Set maximum value from the list to the overall damage regime tracker for simplicity in outputs

                    # System removal #
                    # Optional section to remove mangroves below a set threshold height if they have suffered mortality
                    if sysremactive == 'Y': # Check if user has activated system removal
                        if vegtrkht[b,j] <= sysremthresh and any(value == 2 for value in dmgreglist[(b,j)]): # If stem height is below system removal threshold and has suffered mortality from any of the events
                            for item in vegtrklist: # Set all veg properties to zero (including the tracker)
                                item[b,j] = 0 
                            vegtrkht[b,j] = 0
                            vegtrkpnechk[b,j,:] = 0
                            vegtrkpnelen[b,j,:] = 0
                            vegtrkpnedia[b,j,:] = 0

        # ======================================================================
        # GROWTH
        # ======================================================================
        # Calculate changes to physical properties
        if growthactive == 'Y': # Check if user has activated the growth stage
            for b in range(0,maxcapacity-1): # For every stem
                
                ## Growth logic ##
                # NOTE: OPTIONAL USER UPDATE - Additional conditions can be added to include sigmoidal or parabolic fitness functions to define the growth logic
                if vegtrkht[b,j] != 0:
                    if growthlogic == 'simplified':
                        # Mortality - no growth
                        if any(value == 2 for value in dmgreglist[(b,j)]):
                            growthfactor = 0
                        # Stress - proportionally reduced growth
                        elif any(value == 1 for value in dmgreglist[(b,j)]):
                            if vegtrkht[b,j] < sapht:
                                growthfactorinund = 1 - (c_instr[b,j] - inundtimeseed1)/(inundtimeseed2 - inundtimeseed1)
                                growthfactordesicc = 1 - (c_destr[b,j] - desicctimeseed1)/(desicctimeseed2 - desicctimeseed1)
                                growthfactorburial = 1 - (c_bustr[b,j] - burialeventtime)/(burialmorttime - burialeventtime)
                                growthfactor = min(growthfactorinund, growthfactordesicc, growthfactorburial, 1)
                            if vegtrkht[b,j] >= sapht:
                                growthfactorinund = 1 - (c_instr[b,j] - inundtimesap1)/(inundtimesap2 - inundtimesap1)
                                growthfactordesicc = 1 - (c_destr[b,j] - desicctimesap1)/(desicctimesap2 - desicctimesap1)
                                growthfactorburial = 1 - (c_bustr[b,j] - burialeventtime)/(burialmorttime - burialeventtime)
                                growthfactor = min(growthfactorinund, growthfactordesicc, growthfactorburial, 1)
                        # Recovery - 3/4 growth
                        elif any(value == 0.5 for value in dmgreglist[(b,j)]):
                            growthfactor = 0.75
                        # Healthy - normal growth
                        elif all(value == 0 for value in dmgreglist[(b,j)]):
                            growthfactor = 1
                        # Pneumatophores and cable roots - growth only halts if desiccation occurs (inundation and burial still enable pneums to grow)
                        # Mortality from any event - no growth
                        if any(value == 2 for value in dmgreglist[(b,j)]):
                            growthfactorrt = 0
                        # Else if there is no mortality but stressed from desiccation - half growth
                        elif dmgregde[b,j] == 1:
                            growthfactorrt = 0.5
                        # Else if there is no mortality nor stress from desiccation, but in recovery from desiccation - 3/4 growth
                        elif dmgregde[b,j] == 0.5:
                            growthfactorrt = 0.75
                        # Else if no impacts from desiccation - normal growth
                        elif dmgregde[b,j] == 0:
                            growthfactorrt = 1
                
                ## Seedling phase ##
                if vegtrkht[b,j] < sapht and vegtrk[b,j] == 1: # Stem height less than defined sapling height
                    if vegtrkht[b,j] == 0:
                        growthfactor = 1 # Set growth factor to be 1 if the seed has not yet established a shoot
                        growthfactorrt = 1 # Similarly for root growth factor
                    # Calculate cable/tap root length
                    # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for cable/tap root length growth
                    if vegtrkrt[b,j] == 0: # Determine the initial value once established
                        vegtrkrt[b,j] = startrt # Initial value
                    if vegtrkrt_seedling_growth == 'balke2015':
                        vegtrkrt[b,j] += growthfactorrt * (0.003 * mtpervt / (24*3600)) # Update root growth for these anchored seeds (gradient of Fig. 3 in Balke et al., 2015)
                    if isinstance(vegtrkrt_seedling_growth, (int, float)):
                        vegtrkrt[b,j] += growthfactorrt * (vegtrkrt_seedling_growth * mtpervt / (24*3600))
                    # Calculate cable/tap root diameter
                    # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for cable/tap root diameter growth
                    if seedrtdia2len == 'basyuni2018':
                        vegtrkrtdia[b,j] = vegtrkrt[b,j] / 43 * mtpervt / (24*3600) # Calculate root diameter
                    if isinstance(seedrtdia2len, (int, float)):
                        vegtrkrtdia[b,j] = vegtrkrt[b,j] / seedrtdia2len * mtpervt / (24 * 3600)
                    # Determine if tap root length is sufficient for stem growth
                    # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for cable/tap root length requirement for stem growth
                    if vegtrkrt_thresh == 'balke2015':
                        seedrtthresh = 0.005 # Approx. 2 days of root growth required to resist dislodgement due to inundation alone (Balke et al., 2015)
                    if isinstance(vegtrkrt_thresh, (int, float)):
                        seedrtthresh = vegtrkrt_thresh
                    if vegtrkrt[b,j] >= seedrtthresh: # If the root length is above threshold for stem growth
                        stemtrk[b,j] = 1 # Set stem tracker to 1 now that stem has begun to grow
                        # Calculate stem diameter
                        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for stem diameter growth
                        if vegtrkdia[b,j] == 0: # Determine the initial value once established
                            vegtrkdia[b,j] = startdia # Initial value (user input)
                        if vegtrkdiagrowth == 'hastuti2018':
                            vegtrkdia[b,j] += growthfactor * (0.00004 * mtpervt / (24*3600)) # Calculate shoot diameter
                        if isinstance(vegtrkdiagrowth, (int, float)):
                            vegtrkdia[b,j] += growthfactor * (vegtrkdiagrowth * mtpervt / (24*3600)) # Calculate shoot diameter    
                        # Calculate stem height
                        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for stem height to diameter relationship
                        if seedht2dia == 'jacotot2019': # This relationship is based on the basal diameter and stem height growth for A. marina under natural tidal flooding in New Caledonia (Jacotot et al., 2019)
                            vegtrkht[b,j] = vegtrkdia[b,j] * 25 # Seedling height in the order of 25 x diameter
                        if isinstance(seedht2dia, (int, float)): # If the user has defined a single numerical value instead of an equation
                            vegtrkht[b,j] = vegtrkdia[b,j] * seedht2dia
                        # Calculate pneumatophore growth
                        # Assume that pneumatophores can begin growing whenever the seedling has a stem and the cable roots reach the appropriate length
                        # If the cable root length has reached the next interval for the generation of a pneumatophore
                        if floor(vegtrkrt[b,j]/pneumint) >= 1 and vegtrkpnechk[b,j,floor(vegtrkrt[b,j]/pneumint)-1] == 0: # If the cable root reaches next interval and tracker matrix has no pneum
                            vegtrkpnechk[b,j,floor(vegtrkrt[b,j]/pneumint)-1] = 1 # This assigns a value of 1 to the tracker matrix to highlight the presence of a pneumatophore
                        # For all pneumatophores that exist and are below the maximum length and diameter, calculate the new length/height and new diameter
                        for p in range(maxpneums):
                            if vegtrkpnechk[b,j,p] == 1 and maxwd[j] >= pnelenreq * (vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j])): # If a pneum exists and the max water depth is >= the height of the pneum (since they grow to reach oxygen)
                                # Calculate pneumatophore diameter growth
                                if vegtrkpnedia[b,j,p] == 0: # Determine the initial value once established
                                    vegtrkpnedia[b,j,p] = startpnedia # Initial value
                                # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for pneumatophore diameter growth
                                if vegtrkpnedia[b,j,p] < maxpnedia: # If the base diameter is less than max base diameter
                                    if isinstance(vegtrkpnedia_adult_growth, (int, float)):
                                        vegtrkpnedia[b,j,p] += growthfactorrt * (vegtrkpnedia_adult_growth * mtpervt / (24*3600))
                                # Calculate pneumatophore height/length growth
                                # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for pneumatophore height growth
                                if pnelen2pnedia == "jereznova2022":
                                    vegtrkpnelen[b,j,p] += growthfactorrt * (10.3 * vegtrkpnedia_adult_growth * mtpervt / (24*3600)) # metres/day
                                if isinstance(pnelen2pnedia, (int, float)):
                                    vegtrkpnelen[b,j,p] += growthfactorrt * (pnelen2pnedia * vegtrkpnedia_adult_growth * mtpervt / (24*3600)) # metres/day
                                # Update aboveground pneumatophore height/length
                                if vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j]) > 0: # If the above ground pneumatophore length is positive, then update this. Otherwise leave at 0
                                    agpnelen[b,j,p] = vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j]) # Update above ground pneumatophore length
                                 
                ## Sapling / Adult phase ##
                if vegtrkht[b,j] >= sapht: # Stems that have grown beyond the seedling stage 
                    # Allometric relationships
                    # Calculate change in stem diameter
                    # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for stem diameter growth
                    if vegtrkdia[b,j] < maxstemdia:
                        if vegtrkdia_adult_growth == 'rajkumar2017':
                            vegtrkdia[b,j] += growthfactor * (0.00003 * mtpervt / (24*3600))
                        if isinstance(vegtrkdia_adult_growth, (int, float)):
                            vegtrkdia[b,j] += growthfactor * (vegtrkdia_adult_growth * mtpervt / (24*3600))
                    # Calculate change in stem height
                    # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for stem height to diameter relationship
                    if vegtrkht[b,j] < maxstemht:
                        if dbh2ht == 'thampanya2006':
                            vegtrkht[b,j] = (2.1 + 0.9 * vegtrkdia[b,j]) / 1.45 # Simultaneous equations resolved in Thampanya 2006 dissertation for A. marina
                        if isinstance(dbh2ht, (int, float)): # If user defines a constant value as the relationship between DBH and height
                            vegtrkht[b,j] = dbh2ht * vegtrkdia[b,j]
                    # Calculate change in cable root length
                    if vegtrkrt[b,j] < maxrtlen:
                        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for cable root length growth
                        if isinstance(vegtrkrt_adult_growth, (int, float)):
                            vegtrkrt[b,j] += growthfactorrt * (vegtrkrt_adult_growth * mtpervt / (24*3600))
                    # Calculate change in cable root diameter
                    if vegtrkrtdia[b,j] < maxrtdia:
                        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for cable root diameter growth
                        if vegtrkrtdia_adult_growth == 'allometric1':
                            vegtrkrtdia[b,j] = vegtrkrt[b,j] / 200 # Assumed growth rate to be 1/200 of the cable root length
                        if isinstance(vegtrkrtdia_adult_growth, (int, float)):
                            vegtrkrtdia[b,j] += growthfactorrt * (vegtrkrtdia_adult_growth * mtpervt / (24*3600))
                    # Calculate change in pneumatophore length
                    # If the cable root length has reached the next interval for the generation of a pneumatophore
                    if floor(vegtrkrt[b,j]/pneumint) >= 1 and vegtrkpnechk[b,j,floor(vegtrkrt[b,j]/pneumint)-1] == 0: # If the cable root reaches next interval and tracker matrix has no pneum
                        vegtrkpnechk[b,j,floor(vegtrkrt[b,j]/pneumint)-1] = 1 # This assigns a value of 1 to the tracker matrix to highlight the presence of a pneumatophore
                    # For all pneumatophores that exist and are below the maximum length and diameter, calculate the new length/height and new diameter
                    for p in range(maxpneums):
                        if vegtrkpnechk[b,j,p] == 1 and maxwd[j] >= pnelenreq * (vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j])): # If a pneum exists and the max water depth is >= the height of the pneum (since they grow to reach oxygen)
                            # Calculate pneumatophore diameter growth
                            if vegtrkpnedia[b,j,p] == 0: # Determine the initial value once established
                                vegtrkpnedia[b,j,p] = startpnedia # Initial value
                            # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for pneumatophore diameter growth
                            if vegtrkpnedia[b,j,p] < maxpnedia: # If the base diameter is less than max base diameter
                                if isinstance(vegtrkpnedia_adult_growth, (int, float)):
                                    vegtrkpnedia[b,j,p] += growthfactorrt * (vegtrkpnedia_adult_growth * mtpervt / (24*3600))
                            # Calculate pneumatophore height/length growth
                            # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for pneumatophore height growth
                            if pnelen2pnedia == "jereznova2022":
                                vegtrkpnelen[b,j,p] += growthfactorrt * (10.3 * vegtrkpnedia_adult_growth * mtpervt / (24*3600)) # metres/day
                            if isinstance(pnelen2pnedia, (int, float)):
                                vegtrkpnelen[b,j,p] += growthfactorrt * (pnelen2pnedia * vegtrkpnedia_adult_growth * mtpervt / (24*3600)) # metres/day
                            # Update aboveground pneumatophore height/length
                            if vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j]) > 0: # If the above ground pneumatophore length is positive, then update this. Otherwise leave at 0
                                agpnelen[b,j,p] = vegtrkpnelen[b,j,p] - cablertdepth - (blh[j] - vegtrkbed[b,j]) # Update above ground pneumatophore length
                       
                    # Autochthonous accretion
                    # Calculate autochthonous accretion over the time step
                    if autochactive == "Y": # Check if user has activated autochthonous accretion
                        vegtrkaut[b,j] = autoch * model_years # This is used to record the autochthonous accretion for each stem for the time step
            
            # The added autochthonous accretion is divided by the stem capacity for the grid cell to obtain the proportional autochthonous accretion that is added to the bed level  
            if autochactive == "Y":
                autochgrd[j] = np.sum(vegtrkaut[:,j]) / maxcapacity # This is equal to calculating the average of non-zero autoch accretion values in grid cell j, and then multiplying it by the proportion of accretion (i.e., no. of stems / maxcapacity)
            
    # ======================================================================
    # FUNCTIONALITY
    # ======================================================================
    # This stage is used to calculate the parameters that are returned to DFM, and those that are returned to the user
    # In the case of DFM, these parameters are returned for use in the Baptist roughness equations   
    ## Parameters returned to DFM ##
    if funcactive == 'Y': # Check if user activated the functionality stage
        # For the purposes of calculating projected volumes, areas, and drag coefficients, adopt a threshold for water depth
        maxwd[maxwd < epshu] = 0 # Set all water depths below epshu to 0
        
        # Stems
        # Calculate the stem volume for each grid cell (j)
        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for projected volume of stems
        if projvolstem == "cylinder":
            stem_volume_per_j = np.pi * vegtrkdia ** 2 * np.minimum(np.expand_dims(maxwd, axis=0), vegtrkht) / 4
        # Sum along the first axis for each grid cell (j)
        vmstem = np.sum(stem_volume_per_j, axis=0)
        # Calculate the stem area for each grid cell (j)
        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for projected area of stems
        if projareastem == "cylinder":
            stem_area_per_j = vegtrkdia * np.minimum(np.expand_dims(maxwd, axis=0), vegtrkht)
        amstem = np.sum(stem_area_per_j, axis=0)
        # Stem density
        stemtrk[np.where(vegtrkht != 0)[0], np.where(vegtrkht != 0)[1]] = 1
        rnveg = np.sum(stemtrk, axis = 0) / cellarea

        # Pneumatophores
        # Expand the dimensions of wd to match the shape of vegtrkpnedia
        wd_expanded = np.expand_dims(maxwd, axis=(0, 2))
        # Pneumatophore length
        # The maximum is used to make sure that negative root lengths (below-ground) are not captured in the area, volume, and drag coefficient calculations
        agpnelen = np.maximum(agpnelen, 0)
        # Calculate the volume of each pneumatophore for each grid cell (j)
        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for projected volume of pneumatophores
        if projvolpneum == 'jereznova2022':
            # Pneumatophores modelled with varying diameters in the vertical direction
            vegtrkpnedia50 = ((-0.0007 * agpnelen * 1000 + 0.9898) * vegtrkpnedia) # Calculate pneumatophore diameter at the midpoint
            vegtrkpnedia90 = ((-0.0012 * agpnelen * 1000 + 0.8954) * vegtrkpnedia) # Calculate pneumatophore diameter at 90% of its height
            # Create conditions where the water level corresponds to different heights of the pneumatophores
            # These conditions enable element-wise operations for the two arrays
            conditionvol1 = np.greater_equal(wd_expanded, agpnelen)
            conditionvol2 = np.greater(wd_expanded, 0.9 * agpnelen)
            conditionvol3 = np.greater(wd_expanded, 0.5 * agpnelen)
            conditionvol4 = np.less_equal(wd_expanded, 0.5 * agpnelen)
            if np.any(conditionvol1):
                volpneum = (np.pi * agpnelen / 12) * (vegtrkpnedia90 ** 2 * 0.1 + (vegtrkpnedia50 ** 2 - vegtrkpnedia90 ** 2) * 0.4 + (vegtrkpnedia ** 2 - vegtrkpnedia50 ** 2) * 0.5)
            elif np.any(conditionvol2):
                vegtrkpnediawl = vegtrkpnedia90 * (1 - ((0.1 * agpnelen - (agpnelen - wd_expanded)) / (0.1 * agpnelen))) # Calculate diameter of pneumatophore at water level
                volpneum = (np.pi / 12) * ((vegtrkpnedia90 ** 2 - vegtrkpnediawl ** 2) * (0.1 * agpnelen - (agpnelen - wd_expanded)) + (vegtrkpnedia50 ** 2 - vegtrkpnedia90 ** 2) * 0.4 * agpnelen + (vegtrkpnedia ** 2 - vegtrkpnedia50 ** 2) * 0.5 * agpnelen)
            elif np.any(conditionvol3):
                vegtrkpnediawl = vegtrkpnedia90 + (1 - ((0.4 * agpnelen - (0.9 * agpnelen - wd_expanded)) / (0.1 * agpnelen))) * (vegtrkpnedia50 - vegtrkpnedia90) # Calculate diameter of pneumatophore at water level
                volpneum = (np.pi / 12) * ((vegtrkpnedia50 ** 2 - vegtrkpnediawl ** 2) * (0.4 * agpnelen - (0.9 * agpnelen - wd_expanded)) + (vegtrkpnedia ** 2 - vegtrkpnedia50 ** 2) * 0.5 * agpnelen)
            elif np.any(conditionvol4):
                vegtrkpnediawl = vegtrkpnedia50 + (1 - ((0.5 * agpnelen - (0.5 * agpnelen - wd_expanded)) / (0.5 * agpnelen))) * (vegtrkpnedia - vegtrkpnedia50) # Calculate diameter of pneumatophore at water level
                volpneum = (np.pi / 12) * ((vegtrkpnedia ** 2 - vegtrkpnediawl ** 2) * (0.5 * agpnelen - (0.5 * agpnelen - wd_expanded)))
        if projvolpneum == 'du2021':
            # Pneumatophores are assumed to be conical
            volpneum = np.pi * vegtrkpnedia ** 2 * np.minimum(wd_expanded, agpnelen) / 12
        # Sum along the first and third axes to compute the total volume for each grid cell (j)
        vmpneum = np.sum(volpneum, axis=(0, 2)) * maxrts
        # Calculate the projected area of each pneum for each grid cell (j)
        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for projected area of pneumatophores
        if projareapneum == 'jereznova2022':
            # Pneumatophores modelled with varying diameters in the vertical direction
            vegtrkpnedia50 = ((-0.0007 * agpnelen * 1000 + 0.9898) * vegtrkpnedia) # Calculate pneumatophore diameter at the midpoint
            vegtrkpnedia90 = ((-0.0012 * agpnelen * 1000 + 0.8954) * vegtrkpnedia) # Calculate pneumatophore diameter at 90% of its height
            # Create conditions where the water level corresponds to different heights of the pneumatophores
            # These conditions enable element-wise operations for the two arrays
            conditionarea1 = np.greater_equal(wd_expanded, agpnelen)
            conditionarea2 = np.greater(wd_expanded, 0.9 * agpnelen)
            conditionarea3 = np.greater(wd_expanded, 0.5 * agpnelen)
            conditionarea4 = np.less_equal(wd_expanded, 0.5 * agpnelen)
            if np.any(conditionarea1):
                areapneum = (agpnelen / 2) * (vegtrkpnedia90 * 0.1 + (vegtrkpnedia50 + vegtrkpnedia90) * 0.4 + (vegtrkpnedia + vegtrkpnedia50) * 0.5)
            elif np.any(conditionarea2):
                vegtrkpnediawl = vegtrkpnedia90 * (1 - ((0.1 * agpnelen - (agpnelen - wd_expanded)) / (0.1 * agpnelen))) # Calculate diameter of pneumatophore at water level
                areapneum = (0.5) * ((vegtrkpnedia90 + vegtrkpnediawl) * (0.1 * agpnelen - (agpnelen - wd_expanded)) + (vegtrkpnedia50 + vegtrkpnedia90) * 0.4 * agpnelen + (vegtrkpnedia + vegtrkpnedia50) * 0.5 * agpnelen)
            elif np.any(conditionarea3):
                vegtrkpnediawl = vegtrkpnedia90 + (1 - ((0.4 * agpnelen - (0.9 * agpnelen - wd_expanded)) / (0.1 * agpnelen))) * (vegtrkpnedia50 - vegtrkpnedia90) # Calculate diameter of pneumatophore at water level
                areapneum = (0.5) * ((vegtrkpnedia50 + vegtrkpnediawl) * (0.4 * agpnelen - (0.9 * agpnelen - wd_expanded)) + (vegtrkpnedia + vegtrkpnedia50) * 0.5 * agpnelen)
            elif np.any(conditionarea4):
                vegtrkpnediawl = vegtrkpnedia50 + (1 - ((0.5 * agpnelen - (0.5 * agpnelen - wd_expanded)) / (0.5 * agpnelen))) * (vegtrkpnedia - vegtrkpnedia50) # Calculate diameter of pneumatophore at water level
                areapneum = (0.5) * ((vegtrkpnedia + vegtrkpnediawl) * (0.5 * agpnelen - (0.5 * agpnelen - wd_expanded)))
        if projareapneum == 'du2021':
            # Pneumatophores are assumed to be conical
            areapneum = 1/2 * vegtrkpnedia * np.minimum(wd_expanded, agpnelen)
        # Sum along the first and third axes to compute the total area for each grid cell (j)
        ampneum = np.sum(areapneum, axis=(0, 2)) * maxrts
        
        # Calculate total volume
        vm = vmstem + vmpneum # Projected volume
        am = amstem + ampneum # Projected area
        vcontrol = maxwd * cellarea # Control volume
        
        # Drag coefficient
        # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for drag coefficient
        if dragcoeff == 'vanmaanen2015':
            cd_no = 0.005 # Set to 0.005 (Beselly et al., 2023)
            am_nonzero = np.nonzero(am)
            L[am_nonzero] = (vcontrol[am_nonzero] - vm[am_nonzero])/am[am_nonzero]
            cdvegsp_updated[am_nonzero] = cd_no + 5 / L[am_nonzero]
            cdvegsp_updated[am == 0] = cd_no
        if isinstance(dragcoeff, (int, float)): # If the user has defined a single numerical value instead of an equation 
            cdvegsp_updated = dragcoeff

        # Autochthonous accretion
        if autochactive == "Y": # Check if user activated autochthonous accretion 
            for j in range(ndx): # Creating a second for-loop here for the grid cell number, so that changes to the bed level don't impact the previous for-loop
                bedlevel[j] += autochgrd[j]       

        # Return parameters to DFM
        model_dfm.set_var("bl", bedlevel)
        model_dfm.set_var('rnveg',rnveg)
        model_dfm.set_var('stemheight',hv_stem)
        model_dfm.set_var('diaveg',stemdia)
        model_dfm.set_var('Cdvegsp',cdvegsp_updated)

    # Reset statistical variables
    # Reset variable to 0 to create new statistics for the next timestep
    is_maxvalsnd.fill(0.)
    model_dfm.set_var('is_maxvalsnd',is_maxvalsnd)

    # Biomass
    # This additional stage is used to calculate the biomass of the vegetation for each grid cell
    if biomassactive == 'Y': # Check if user activated biomass 
        for j in range(ndx):
            for b in range(0,maxcapacity-1):
                # Above-ground biomass
                # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for above-ground biomass
                # Calculate above-ground biomass for stems
                if agb == 'comleymcguiness2005': # Note this excludes pneumatophores
                    if vegtrkdia[b,j] < 0.35: # Comley & McGuiness 2005 equation is applied for mangroves with a maximum diameter of 0.35m
                        agb_stem[b,j] = 0.308 * (vegtrkdia[b,j] * 100) ** 2.11 # kg. Sum of above-ground biomass for stems (not pneumatophores)
                if agb == 'komiyama2005':
                    if vegtrkdia[b,j] > 0.05 and vegtrkdia[b,j] < 0.49: # Komiyama 2005 common equation based on sample trees with diameters between 0.05m and 0.49m
                        meanwooddensity = 0.506 # t/m3. Avicennia alba. No density measurements recorded here for A. marina.
                        agb_stem[b,j] = 0.251 * meanwooddensity * (vegtrkdia[b,j] * 100) ** 2.46 # kg. Sum of above-ground biomass for stems
                # Calculate above-ground biomass for pneumatophores separately
                agb_pne_sum[b,j] = np.sum(agb_pne[b,j,:]) # Calculate the sum of above-ground biomass for all pneumatophores for each stem and store it in agb_pne_sum
                # Below-ground Biomass
                # Note that belowground biomass results from allometric relationships are often 40% larger than those from field tests (Adame et al., 2017)
                # This is true for all equations presented here, but note that practitioners should derive their own allometric relationships
                # NOTE: OPTIONAL USER UPDATE - additional if statements can be included below for below-ground biomass
                if bgb == 'comleymcguiness2005': # Note that this equation includes pneumatophores, but also likely underestimates the biomass by only including a 2m root radius
                    if vegtrkdia[b,j] < 0.35: # Comley & McGuiness 2005 equation is applied for mangroves with a maximum diameter of 0.35m.
                        bgb_stem[b,j] = 1.28 * (vegtrkdia[b,j] * 100) ** 1.17 # kg. Sum of below-ground biomass for stems
                if bgb == 'komiyama2005':
                    if vegtrkdia[b,j] > 0.05 and vegtrkdia[b,j] < 0.45: # Komiyama 2005 common equation based on sample trees with diameters between 0.05m and 0.49m
                        bgb_stem[b,j] = 0.199 * meanwooddensity ** 0.899 * (vegtrkdia[b,j] * 100) ** 2.22 # kg. Sum of below-ground biomass.        
            # Calculate total above- and below-ground biomass
            agb_total[j] = np.sum(agb_stem[:,j]) + np.sum(agb_pne_sum[:,j]) # Total above-ground biomass per grid cell by summing all stem agb and pneumatophore agb per grid cell and timestep
            bgb_total[j] = np.sum(bgb_stem[:,j])

    # ======================================================================
    # PREPARE VARIABLES FOR PLOTTING
    # ======================================================================
    # Calculate the non-zero percentiles of stem and pneumatophore heights and diameters #
    # Calculate percentile for each unique value along the second dimension

    # Stem height and bed level
    hv_stem_nzindices = np.nonzero(vegtrkht)
    unique_indices = np.unique(hv_stem_nzindices[1]) # Unique indices along the second dimension
    for idx in unique_indices:
        # Find non-zero values corresponding to current unique index
        idx_mask = hv_stem_nzindices[1] == idx
        non_zero_hvsvalues = vegtrkht[hv_stem_nzindices[0][idx_mask], hv_stem_nzindices[1][idx_mask]]
        non_zero_hvsvalues = non_zero_hvsvalues[non_zero_hvsvalues != 0] # Filter out zero values
        bedvalues = vegtrkbed[hv_stem_nzindices[0][idx_mask], hv_stem_nzindices[1][idx_mask]]
        bedvalues = bedvalues[bedvalues != 0] # Filter out zero values
        percentile_bed = np.percentile(bedvalues, percentileveg)
        # Calculate percentile for non-zero values
        percentile_value = np.percentile(non_zero_hvsvalues, percentileveg)
        # Update corresponding entry in hv_stem
        hv_stem[idx] = percentile_value
        vegbed[idx] = percentile_bed

    # Stem diameter
    stemdia_nzindices = np.nonzero(vegtrkdia)
    unique_sdindices = np.unique(stemdia_nzindices[1]) # Unique indices along the second dimension
    for idx in unique_sdindices:
        # Find non-zero values corresponding to current unique index
        idx_sdmask = stemdia_nzindices[1] == idx
        non_zero_sdvalues = vegtrkdia[stemdia_nzindices[0][idx_sdmask], stemdia_nzindices[1][idx_sdmask]]
        non_zero_sdvalues = non_zero_sdvalues[non_zero_sdvalues != 0] # Filter out zero values
        # Calculate percentile for non-zero values
        percentile_sdvalue = np.percentile(non_zero_sdvalues, percentileveg)
        # Update corresponding entry in stemdia
        stemdia[idx] = percentile_sdvalue

    # Mangrove presence
    for p in range(ndx): # For every grid cell
        vegcheck[p] = check_column(vegtrk,p) # Assign a value of 1 if there is at least one stem in the grid cell

    # Above-ground pneumatophore height/length
    agpnelen_nzindices = np.nonzero(agpnelen) # Note that this is non-zero, but there is a condition below which ensures hv_pneum only returns values >= 0 (default value is 0)
    if agpnelen_nzindices[0].size > 0:
        # Calculate percentile for each unique value along the second dimension
        unique_agplindices = np.unique(agpnelen_nzindices[1]) # Unique indices along the second dimension
        for idx in unique_agplindices:
            # Find non-zero values corresponding to current unique index
            idx_agplmask = agpnelen_nzindices[1] == idx
            non_zero_agplvalues = agpnelen[agpnelen_nzindices[0][idx_agplmask], agpnelen_nzindices[1][idx_agplmask]]
            non_zero_agplvalues = non_zero_agplvalues[non_zero_agplvalues > 0] # Filter out values at or below zero. Pneums that haven't yet surfaced might be causing negative projected area and volume results
            # Calculate percentile for non-zero values
            percentile_agplvalue = np.percentile(non_zero_agplvalues, percentileveg)
            # Update corresponding entry in hv_pneum
            if percentile_agplvalue > 0:
                hv_pneum[idx] = percentile_agplvalue

    # Pneumatophore diameter
    pneumdia_nzindices = np.nonzero(vegtrkpnedia)
    if pneumdia_nzindices[0].size > 0:
        # Calculate percentile for each unique value along the second dimension
        unique_pdindices = np.unique(pneumdia_nzindices[1]) # Unique indices along the second dimension
        for idx in unique_pdindices:
            # Find non-zero values corresponding to current unique index
            idx_pdmask = pneumdia_nzindices[1] == idx
            non_zero_pdvalues = vegtrkpnedia[pneumdia_nzindices[0][idx_pdmask], pneumdia_nzindices[1][idx_pdmask]]
            non_zero_pdvalues = non_zero_pdvalues[non_zero_pdvalues != 0] # Filter out zero values
            # Calculate percentile for non-zero values
            percentile_pdvalue = np.percentile(non_zero_pdvalues, percentileveg)
            # Update corresponding entry in pneumdia
            pneumdia[idx] = percentile_pdvalue

    # Pneumatophore density
    n_pne = (np.sum(vegtrkpnechk, axis = (0,2)) / cellarea) * maxrts

    # Set mortality variables for plotting to refer to the maximum value (worst case) in each grid cell
    dmgregingc = np.max(dmgregin, axis=0)
    dmgregdegc = np.max(dmgregde, axis=0)
    dmgregbugc = np.max(dmgregbu, axis=0)
    dmgregsegc = np.max(dmgregse, axis=0)
    dmgregallgc = np.max(dmgregall, axis=0)
    blgc = blh # Set plotting variable at each grid cell to be equal to bed level at the current timestep i.   

    # ======================================================================
    # PLOTTING
    # ======================================================================
    # Plotting key variables at the user-specified intervals
    if i % plotint == 0 or i == nt-1: # Check whether the timestep corresponds to the defined interval or the last timestep in the model run
        print('plotting for timestep: ', i)
        # This section is to save figures of variables
        # By including this section before the model dimr is finalised, the figures are saved in the dflowfm folder

        ## Cross-section figures ##
        # Note that these cross-sections are based on user defined y-coordinates for the cross-section, xsecycoords (see INPUTS)
        # Define colours to represent the bed, stems, pneumatophores, and water in the cross-section
        sand_colour = '#BA934E'
        stem_colour = '#7B9A6E' # for healthy mangroves
        stem_line_colour = '#D2B48C'
        pneum_colour = '#443D01'
        pneum_line_colour = '#443D01'
        water_colour = '#DBECEE'
        
        # Calculate the average age for each grid cell, multiply by model_months, and round to the nearest integer
        ages = np.round(np.mean(vegtrkage, axis=0) * model_months).astype(int)       

        for idx, y_coord_value in enumerate(xsecycoords):
            matchingindices = np.where(y_coord == y_coord_value) # Identify the corresponding indices in the y-coordinates
            matchxcoord = x_coord[matchingindices] # Extract corresponding x-coordinates
            matchblh = blh[matchingindices] # Extract corresponding bed levels
            matchvegbed = vegbed[matchingindices] # Extract corresponding bed level
            matchhvstem = hv_stem[matchingindices] + matchvegbed # Update stem height relative to bed level
            matchhvpneum = hv_pneum[matchingindices] + matchvegbed # Update pneumatophore height relative to bed level
            matchmaxwd = maxwd[matchingindices] + matchblh # For showing the max water level at the timestep of the plot
            matchdmgregallgc = dmgregallgc[matchingindices] # Extract corresponding mangrove health status
            # Filter for non-zero stem heights and pneum heights
            non_zero_stem_indices = (matchhvstem > 0) & (matchhvstem > matchblh)
            non_zero_pneum_indices = (matchhvpneum > 0) & (matchhvpneum > matchblh)
            non_zero_matchxcoord_stem = matchxcoord[non_zero_stem_indices]
            non_zero_matchhvstem = matchhvstem[non_zero_stem_indices]
            non_zero_matchxcoord_pneum = matchxcoord[non_zero_pneum_indices]
            non_zero_matchhvpneum = matchhvpneum[non_zero_pneum_indices]
            # Filter for water depths above the bed level
            water_above_bed_indices = matchmaxwd > matchblh
            non_zero_matchxcoord_maxwd = matchxcoord[water_above_bed_indices]
            non_zero_maxwd = matchmaxwd[water_above_bed_indices]
            matchblh_above = matchblh[water_above_bed_indices]
            # Create 2D matrix for bed levels
            xsecblh = np.column_stack((matchxcoord, matchblh))
            # Export to csv
            pandas.DataFrame(xsecblh, columns=['x-coord','bed-level']).to_csv(os.path.join(csvfolder_path, f'{today}_xsecblh_{y_coord_value}_{i}.csv'),index=False) # Output to CSV
            # Create plot
            plt.figure(figsize=(10, 4))
            # Plot bed levels
            plt.plot(matchxcoord, matchblh, linestyle='-', label=f'Bed level', color=sand_colour) # Line graph with points marked
            # Fill the area beneath the bed level to the bottom of the graph
            plt.fill_between(matchxcoord, matchblh, xplotminbl, color=sand_colour, alpha=0.5)
            # Plot water levels
            plt.plot(non_zero_matchxcoord_maxwd, non_zero_maxwd, linestyle='-', color=water_colour, label=f'Water Levels at Y = {y_coord_value}')
            # Fill the area between the water level and bed level
            plt.fill_between(non_zero_matchxcoord_maxwd, matchblh_above, non_zero_maxwd, color=water_colour, alpha=0.5)
            # Add mangroves to plot with colour coded canopies based on health status if chosen by the user. Default value is green otherwise
            if includestatus == "N":
                plt.scatter(non_zero_matchxcoord_stem, non_zero_matchhvstem, color=stem_colour, label=f'Stem Heights (Healthy)', s=100)
            if includestatus == "Y":
                # Plotting the canopy circles with colors based on mortality
                canopy_colours = [map_mortality_to_colour(mort) for mort in matchdmgregallgc[non_zero_stem_indices]]
                plt.scatter(non_zero_matchxcoord_stem, non_zero_matchhvstem, color=canopy_colours, s=100)
                # Add invisible scatter plots for legend labels
                plt.scatter([], [], color=colour_healthy, label='Healthy')
                plt.scatter([], [], color=colour_recovery, label='Recovery')
                plt.scatter([], [], color=colour_stress, label='Stress')
                plt.scatter([], [], color=colour_mortality, label='Mortality')
            # Plot pneumatophore heights as points
            plt.scatter(non_zero_matchxcoord_pneum, non_zero_matchhvpneum, color=pneum_colour, label=f'Pneumatophore Heights', marker='d')
            # Plot vertical lines between bed levels and stem heights
            for x, bl, sh in zip(non_zero_matchxcoord_stem, matchblh[non_zero_stem_indices], non_zero_matchhvstem):
                plt.plot([x,x], [bl,sh], color=stem_line_colour, linestyle='-')
            # Plot vertical lines between bed levels and pneumatophore heights
            for x, bl, ph in zip(non_zero_matchxcoord_pneum, matchblh[non_zero_pneum_indices], non_zero_matchhvpneum):
                plt.plot([x,x], [bl, ph], color=pneum_line_colour, linestyle='-')

            plt.xlabel('Cross-shore distance (m)')
            plt.ylabel('Bed level (m)')
            plt.title(f'Cross-Section at Y = {y_coord_value}')
            plt.grid(False)
            plt.legend()

            # Set the y-axis limits
            plt.ylim(-1, 5) # To be updated based on expected mangrove heights and bed levels

            # Set the x-axis limits to the edge of the data
            min_x = matchxcoord.min()
            max_x = matchxcoord.max()
            plt.xlim(min_x, max_x)

            # Annotate the stem heights with the ages
            if includeages == "Y":
                for x, y, age in zip(non_zero_matchxcoord_stem, non_zero_matchhvstem, ages[matchingindices][non_zero_stem_indices]):
                    plt.text(x, y, str(int(age)), color='white', fontsize=8, ha='center', va='center')
           
            # Save the plot
            plt.savefig(os.path.join(plotfolder_path, f'{today}_xsecblh_{y_coord_value}_{i}.png'))

        ## Plan-view figures and CSV outputs ##        
        # Set up variable list and names
        variables = [vegcheck, hv_stem, stemdia, rnveg, ages, hv_pneum, pneumdia, n_pne, \
                     dmgregingc, dmgregdegc, dmgregbugc, dmgregsegc, dmgregallgc, \
                     vm, am, cdvegsp, blgc, agb_total, bgb_total, maxwd, minwd, maxbss]
        variable_names = ['Mangrove presence (-)', 'Stem height (m)', 'Stem diameter (m)', 'Stem density (stems per m2)', 'Stem age (months)', 'Above ground pneumatophore height (m)', 'Pneumatophore diameter (m)', \
                          'Pneumatophore density (pneumatophores per m2)', 'Inundation mortality (-)', 'Desiccation mortality (-)', 'Burial mortality (-)', 'Senescence mortality (-)', \
                            'Mortality and recovery (-)', 'Projected volume (m3)', 'Projected area (m2)', 'Drag coefficient (-)', 'Bed level (m)', 'Above-ground biomass (kg)', 'Below-ground biomass (kg)', \
                                'Maximum water depth (m)', 'Minimum water depth (m)', 'Maximum bed shear stress (Nm-2)']

        for u, variable in enumerate(variables):
            # Create a figure and axis for each plot
            fig, ax = plt.subplots(figsize=(12,6))
            # Create a colour map and normalise it based on stem heights
            cmap = plt.get_cmap('YlGn') # See the codes for colour gradients at the link here: https://matplotlib.org/stable/users/explain/colors/colormaps.html
            norm = Normalize(vmin=np.min(variable), vmax=np.max(variable))
            cmap2 = plt.get_cmap('YlOrRd') # Yellow to red colour map for the damage regime plots
            norm2 = Normalize(vmin=0, vmax=2)
            cmap3 = plt.get_cmap('YlGnBu') # Yellow to green to blue colour map for the max water depth plot
            norm3 = Normalize(vmin=np.min(variable), vmax=np.max(variable))
            cmap4 = plt.get_cmap('RdPu') # Red to purple colour map for the max bed shear stress plot
            norm4 = Normalize(vmin=np.min(variable), vmax=np.max(variable))

            if u < 8 or (u > 12 and u < 19): # For plots that aren't damage regimes, set these to have a yellow to green gradient
                # Iterate through grid cell centres and plot them with colour coding
                for x_center, y_center, data in zip(x_coord, y_coord, variable):
                    # Define the coordinates of the rectangle (grid cell)
                    x_left = x_center - 2.5
                    x_right = x_center + 2.5
                    y_bottom = y_center - 2.5
                    y_top = y_center + 2.5
                    color = cmap(norm(data))
                    # Plot the rectangle with colour based on stem height
                    rect = plt.Rectangle((x_left, y_bottom), x_right - x_left, y_top - y_bottom, color=color, edgecolor='black')
                    ax.add_patch(rect)
                # Add colour bar
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)

            elif u >= 8 and u <= 12: # For damage regime plots, set these to have a green to yellow to red gradient
                # Iterate through grid cell centres and plot them with colour coding
                for x_center, y_center, data in zip(x_coord, y_coord, variable):
                    # Define the coordinates of the rectangle (grid cell)
                    x_left = x_center - 2.5
                    x_right = x_center + 2.5
                    y_bottom = y_center - 2.5
                    y_top = y_center + 2.5
                    color = cmap2(norm2(data))
                    # Plot the rectangle with colour based on stem height
                    rect = plt.Rectangle((x_left, y_bottom), x_right - x_left, y_top - y_bottom, color=color, edgecolor='black')
                    ax.add_patch(rect)
                # Add colour bar
                sm = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2)

            elif u == 19 or u == 20: # For max water depth plot
                # Iterate through grid cell centres and plot them with colour coding
                for x_center, y_center, data in zip(x_coord, y_coord, variable):
                    # Define the coordinates of the rectangle (grid cell)
                    x_left = x_center - 2.5
                    x_right = x_center + 2.5
                    y_bottom = y_center - 2.5
                    y_top = y_center + 2.5
                    color = cmap3(norm3(data))
                    # Plot the rectangle with colour based on stem height
                    rect = plt.Rectangle((x_left, y_bottom), x_right - x_left, y_top - y_bottom, color=color, edgecolor='black')
                    ax.add_patch(rect)
                # Add colour bar
                sm = plt.cm.ScalarMappable(cmap=cmap3, norm=norm3)

            elif u == 21: # For max bed shear stress plot
                # Iterate through grid cell centres and plot them with colour coding
                for x_center, y_center, data in zip(x_coord, y_coord, variable):
                    # Define the coordinates of the rectangle (grid cell)
                    x_left = x_center - 2.5
                    x_right = x_center + 2.5
                    y_bottom = y_center - 2.5
                    y_top = y_center + 2.5
                    color = cmap4(norm4(data))
                    # Plot the rectangle with colour based on stem height
                    rect = plt.Rectangle((x_left, y_bottom), x_right - x_left, y_top - y_bottom, color=color, edgecolor='black')
                    ax.add_patch(rect)
                # Add colour bar
                sm = plt.cm.ScalarMappable(cmap=cmap4, norm=norm4)
            
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, label=f'{variable_names[u]}')
            # Set the font size for axis labels
            ax.tick_params(axis='both', labelsize=12)
            # Customise the x-axis to show specific values
            x_values_to_show = [0, 50, 100, 150, 200, 250]
            ax.set_xticks(x_values_to_show, minor=False)
            ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
            x_range = max(x_coord2) - min(x_coord2)
            y_range = max(y_coord2) - min(y_coord2)
            aspect_ratio = x_range / (4*y_range)
            # Set axes scaling to be equal
            ax.set_aspect(aspect_ratio)
            y_values_to_show = [0, 25, 50]
            ax.set_yticks(y_values_to_show, minor=False)
            ax.yaxis.set_minor_locator(plt.MultipleLocator(5))
            ax.grid(True, which='both', linestyle='--', color='gray', alpha=0.7)
            # Add labels and title
            ax.set_xlabel('Cross-shore distance (m)')
            ax.set_ylabel('Long-shore distance (m)')
            ax.set_title(f'Plot {u + 1}, timestep {i}: {variable_names[u]}', fontsize = 16) # Timestep included in the axis title
            plt.xlim(min(x_coord2), max(x_coord2))
            plt.ylim(min(y_coord2), max(y_coord2))
            # Save the plot (.png)
            plt.savefig(os.path.join(plotfolder_path, f'{today}_plot_{u + 1}_{i}_{variable_names[u]}.png')) # Timestep included in the filename
            # Output a csv file
            reshaped_variable = variable.reshape(-1, variable.shape[-1]) # Need to reshape 3D arrays before converting to CSV using DataFrame from pandas
            pandas.DataFrame(reshaped_variable).to_csv(os.path.join(csvfolder_path,f'{today}_plot_{u + 1}_{i}_{variable_names[u]}.csv'),index=False)
            plt.close()

###### End of feedback loop ######
model_dimr.finalize()
print('Model run complete')
# End of script