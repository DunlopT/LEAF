When the DFM project has been created, saved and exported, this folder will be created.
This folder will contain the .dsproj file for the DFM project (PROJECTNAME.dsproj), and subfolders titled PROJECTNAME.dsproj_data, dflowfm, and wave (if D-Waves is also activated).
For the DFM example, Spencer_waves, the following files require updating to couple the DFM model with LEAF (v.1.0) and ensure that vegetation is activated in DFM.

*****
This folder:
1. dimr_config.xml

This file requires updating if the names of the FlowFM.mdu or Waves.mdw files have changed, or if the start time, stop time, or time interval for wave calculations has changed.

If the filenames have changed, update the following lines:
<inputFile>Waves.mdw</inputFile>
<inputFile>FlowFM.mdu</inputFile>

If the times of waves calculations have changed, update the following line:
<time>0 86400 31536000</time>

*****
dflowfm subfolder:
2. FlowFM.mdu

This file contains most parameters for the DFM model run.
Ensure the values of the following parameters are updated:
RefDate - this should match the start date of the model run
TStop - this should match the final stop time of the model run and water level time series
ExtForceFile - this should refer to the set_mangroves.ext file which defines the mangrove stem height
OutputDir - this should refer to the output directory

3. LEAF_inputs_example_Spencer.json

Ensure that the product of the timestep duration (mtpervt) and number of timesteps (nt) matches the TStop value in the FlowFM.mdu file.
Set the initial mangrove activation parameter (iniveg) to "N" if the project site is not to have any initial mangroves.

4. mangroves.pli

This file specifies the coordinates of the area in which mangroves are represented in DFM. The size and changes to these mangroves however, are calculated in the LEAF (v.1.0) model.

5. set_mangroves.ext

This file specifies the initial stem height of mangroves in the final row, VALUE. Note that if the LEAF inputs file has set iniveg = "N", then this file won't be activated.
The row with FILENAME, should refer to mangroves.pli.

6. WaterLevelSpencer.bc

This file should list all water level data and corresponding times. Check that the final time aligns with the product of mtpervt and nt from the LEAF inputs file.

*****
wave subfolder:
7. Waves.mdw

Update the value for the following parameter in this file:
ReferenceDate - check that this lines up with the date that waves are expected to start in the model (given by the number of seconds in the dimr_config.xml file from the start RefDate in FlowFM.mdu file)
COMFile - this file should refer to the FlowFM_com.nc file in the ...PROJECTNAME.dsproj_data/FlowFM/output folder
