When the DFM project has been created, saved and exported, this folder will be created.
This folder will contain the .dsproj file for the DFM project (PROJECTNAME.dsproj), and subfolders titled PROJECTNAME.dsproj_data, dflowfm, and wave (if D-Waves is also activated).
For the DFM example, Spencer_no_waves, the following files require updating to couple the DFM model with LEAF (v.1.0) and ensure that vegetation is activated in DFM.

*****
This folder:
1. dimr_config.xml

This file only requires updating if the name of the FlowFM.mdu file has changed.

If the filename has changed, update the following line:
<inputFile>FlowFM.mdu</inputFile>

*****
dflowfm subfolder:
2. FlowFM.mdu

This file contains most parameters for the DFM model run.
Ensure the values of the following parameters are updated:
RefDate - this should match the start date of the model run
TStop - this should match the final stop time of the model run and water level time series
ExtForceFile - this should refer to the set_mangroves.ext file which defines the mangrove stem height
OutputDir - this should refer to the output directory (e.g., C:\Users\NAME\Documents\PROJECTNAME\PROJECTNAME.dsproj_data\FlowFM\output)

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